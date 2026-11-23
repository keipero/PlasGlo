r"""CO2・CO・O2 の断面積から、電子チャンネル対応表の叩き台を作る。

LXCat 形式の断面積ファイルを Py3bolos のパーサで読み、過程の記述から
分類して、次の2つを書き出す（``SET`` は CO2 の出典 ``ist`` または ``biagi``）。

- ``co2_channel_map_SET.csv``：全過程の分類と、反応への対応案
- ``co2_electrons_profile_SET.json``：``electrons.profile`` 用の断片

CO2 の出典は ``--co2-set`` で選ぶ。``ist``（IST-Lisbon の更新版。既定）は、
解離の生成物と分岐比がデータに明記されており、基準にする。``biagi`` は、
高E/Nの解離（中性双極子解離）や解離性電離を含む、感度の比較用。CO は
IST-Lisbon、O2 は Phelps を共通に使う。

``--check`` を付けると、作った対応表で ``ElectronSwarm`` を初期化し、純CO2と
CO2/CO/O2 混合気で BOLOS が収束するかを確認する。

状態（``status``）は、``explicit``（ファイルの記述に明記）、``assumption``
（コメントからの仮定）、``review``（要確認）、``none``（化学に入れない）。
未確認の点は CSV の ``status`` 列と README に書いている。
"""

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1]
sys.path.insert(0, str(MODEL))
sys.path.insert(0, str(MODEL / "Py3bolos"))

from bolos import parser  # noqa: E402

# LXCat からダウンロードした断面積ファイルのフォルダ。環境変数 CO2_XS_DIR、
# または各スクリプトの --source-dir で指定する（既定値は持たない）。
SOURCE_DIR = (
    Path(os.environ["CO2_XS_DIR"]) if os.environ.get("CO2_XS_DIR") else None
)
# CO2 の出典ごとのファイル名。ist は LXCat（IST-Lisbon、CO2 の完全なセット）
# からダウンロードした Cross section.txt を、この名前で置く。
CO2_FILES = {
    "ist": "ISTLisbon_CO2_20261123.txt",
    "biagi": "Biagi_CO2_1027.txt",
}
OTHER_FILES = {
    "CO": "ISTLisbon_CO_TRINITI_COrot_1027.txt",
    "O2": "Phelps_O2_1027.txt",
}
CO2_SET = "ist"
# 取り込み後にケースJSONから参照する予定の相対パス（ファイルはまだ無い）。
PLANNED_DIR = "../../data/cross_sections/co2"


def sources():
    r"""選択中の出典で、標的ごとの断面積ファイル名を返す。

    Returns
    -------
    dict of str to str
        標的名（CO2、CO、O2）からファイル名への対応。
    """
    return {"CO2": CO2_FILES[CO2_SET], **OTHER_FILES}


def set_source_dir(value=None):
    r"""断面積フォルダを設定し、実在しなければ終了する。

    Parameters
    ----------
    value : str or None, optional
        フォルダのパス。None なら、環境変数 ``CO2_XS_DIR`` の値を使う。

    Returns
    -------
    pathlib.Path
        断面積ファイルのあるフォルダ。

    Raises
    ------
    SystemExit
        フォルダが未指定、または存在しない場合。
    """
    global SOURCE_DIR
    if value:
        SOURCE_DIR = Path(value)
    if SOURCE_DIR is None or not SOURCE_DIR.is_dir():
        raise SystemExit(
            "断面積のフォルダを --source-dir か環境変数 CO2_XS_DIR で指定"
            f"してください（現在の値：{SOURCE_DIR}）。"
        )
    missing = [n for n in sources().values() if not (SOURCE_DIR / n).is_file()]
    if missing:
        raise SystemExit(f"{SOURCE_DIR} に次のファイルがありません：{missing}")
    return SOURCE_DIR


def sha256(path):
    r"""ファイルのSHA-256を返す。

    Parameters
    ----------
    path : pathlib.Path
        ハッシュを計算するファイル。

    Returns
    -------
    str
        16進文字列のSHA-256。
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def comment_line(collision):
    r"""断面積コメントのうち、過程の説明にあたる最初の行を返す。

    Parameters
    ----------
    collision : dict
        Py3bolos のパーサが返す1過程分の辞書。

    Returns
    -------
    str
        ``COMMENT:`` 行の本文。無ければ空文字列。
    """
    for line in collision["comment"].splitlines():
        if line.startswith("COMMENT:"):
            return line[len("COMMENT:"):].strip()
    return ""


def classify_co2_ist(c):
    r"""IST-Lisbon のCO2の過程を、生成物の記述に基づいて分類する。

    過程の記述に生成物が明記されているため、状態は ``explicit`` とする。

    Parameters
    ----------
    c : dict
        CO2の1過程。

    Returns
    -------
    tuple
        分類名、反応式の案、状態（explicit/none）、備考。
    """
    kind, product = c["kind"], c.get("product") or ""
    if kind == "ELASTIC":
        return "elastic", "", "none", "運動量移行。EEDFのみ"
    if kind == "ATTACHMENT":
        return (
            "attachment",
            "E + CO2 -> CO + O^-",
            "explicit",
            "生成物は CO(X) + O(-,X)（Lowke, Phelps, Irwin 1973）",
        )
    if kind == "IONIZATION":
        return (
            "ion_CO2+",
            "E + CO2 -> E + E + CO2^+",
            "explicit",
            "Itikawa 2002。解離性電離は無い",
        )
    if "O(1D)" in product:
        return (
            "diss_7eV",
            "E + CO2 -> E + CO + O(1D)",
            "explicit",
            "7 eV の過程の15%（Liu 2025）。励起種 O(1D) を機構に入れる必要",
        )
    if "CO(a3P)" in product:
        return (
            "diss_10p5eV",
            "E + CO2 -> E + CO(a3Pi) + O",
            "explicit",
            "全て解離。励起種 CO(a3Pi) を機構に入れる必要。Oは基底3P",
        )
    if product == "CO2(e1)":
        return (
            "exc_7eV_e1",
            "",
            "none",
            "7 eV の過程の85%（解離しない分）。EEDFのみ",
        )
    if "v=" in product:
        return (
            "vibration",
            "",
            "none",
            "振動励起。振動種を追跡しない限りEEDFのみ",
        )
    return "excitation_other", "", "none", "EEDFのみ"


def classify_co2(c):
    r"""Biagi のCO2の過程を、コメントに基づいて分類する。

    Parameters
    ----------
    c : dict
        CO2の1過程。

    Returns
    -------
    tuple
        分類名、反応式の案、状態（assumption/review/none）、備考。
    """
    text, kind = comment_line(c), c["kind"]
    if kind == "ELASTIC":
        return "elastic", "", "none", "運動量移行。EEDFのみ"
    if kind == "ATTACHMENT":
        return (
            "attachment",
            "E + CO2 -> CO + O^-",
            "assumption",
            "Biagiは生成物を記載せず。Itikawa 2002の CO2->CO+O^- を採用",
        )
    if kind == "IONIZATION":
        table = [
            ("Ionication CO2+ from", "ion_CO2+", "E + CO2 -> E + E + CO2^+"),
            ("CO2+ (A2Piu) from", "ion_CO2+_A", "E + CO2 -> E + E + CO2^+"),
            ("CO2+ (B2Sig+u) from", "ion_CO2+_B", "E + CO2 -> E + E + CO2^+"),
            ("Ionization O+ from", "ion_O+", "E + CO2 -> E + E + O^+ + CO"),
            ("Ionization CO+ from", "ion_CO+", "E + CO2 -> E + E + CO^+ + O"),
        ]
        # 「O+」が「O++」に部分一致しないよう、直後の " from" まで含めて照合する。
        for key, name, eq in table:
            if key in text:
                note = (
                    "励起イオン状態を CO2^+ に集約する仮定"
                    if name.endswith(("_A", "_B"))
                    else "解離性電離の生成物は仮定"
                    if name in ("ion_O+", "ion_CO+")
                    else ""
                )
                return name, eq, "assumption", note
        return (
            "ion_other",
            "",
            "none",
            "C+、二価イオン、K殻。27 eV以上で寄与小と見て化学に入れない",
        )
    if "Vibrational" in text:
        return (
            "vibration",
            "",
            "none",
            "振動励起。振動種を追跡しない限りEEDFのみ",
        )
    if "Dissociation 1Delta" in text:
        return (
            "diss_1Delta",
            "E + CO2 -> E + CO + O",
            "assumption",
            "Oの電子状態（3P/1D）はファイルに記載なし。基底3Pと仮定",
        )
    if "Dissociation 1Pi" in text:
        return (
            "diss_1Pi",
            "E + CO2 -> E + CO + O",
            "assumption",
            "Oの電子状態は記載なし。基底3Pと仮定",
        )
    if "Dissociation neutral dipole" in text:
        return (
            "diss_dipole",
            "E + CO2 -> E + CO + O",
            "review",
            "生成物は記載なし。CO+Oと仮定。Magboltz資料で要確認",
        )
    if "Dissociative excitation" in text:
        return (
            "diss_excitation",
            "E + CO2 -> E + CO + O",
            "review",
            "コメントの '1Sig+1Sig+C!' の意味と生成物が不明。要確認",
        )
    return (
        "excitation_other",
        "",
        "review",
        "電子励起。解離に至る分は不明のため、いまはEEDFの損失のみ",
    )


def classify_other(target, c):
    r"""CO・O2の過程を分類する。

    Parameters
    ----------
    target : str
        標的名（CO または O2）。
    c : dict
        1過程分の辞書。

    Returns
    -------
    tuple
        分類名、反応式の案、状態、備考。
    """
    kind, product = c["kind"], c.get("product") or ""
    text = comment_line(c)
    if kind in ("ELASTIC", "EFFECTIVE"):
        return "momentum", "", "none", "運動量移行。EEDFのみ"
    if target == "CO":
        if kind == "ATTACHMENT":
            return (
                "attachment",
                "E + CO -> C + O^-",
                "assumption",
                "生成物は C(X)+O(-)。Cを追跡しない場合は除外",
            )
        if kind == "IONIZATION":
            return "ion_CO+", "E + CO -> E + E + CO^+", "assumption", ""
        if "C(X) + O(3P)" in product:
            return (
                "diss_CO",
                "E + CO -> E + C + O",
                "assumption",
                "C原子を追跡しない場合は除外",
            )
        if "v=0-v=" in product:
            return (
                "vibration",
                "",
                "none",
                "振動励起。超弾性はPlasGloで未対応（順方向のみ使用）",
            )
        return "excitation_other", "", "none", "電子励起。EEDFのみ"
    # O2
    if kind == "ATTACHMENT":
        if not product:
            return (
                "attachment_3body",
                "E + O2 + M -> O2^- + M",
                "assumption",
                "三体付着。密度規格化（空気ケースC27と同じ設定）が必要",
            )
        return "attachment_2body", "E + O2 -> O^- + O", "none", ""
    if kind == "IONIZATION":
        return "ion_O2+", "E + O2 -> E + E + O2^+", "none", ""
    if "6.0eV" in product or "8.4eV" in product:
        return (
            "excitation_dissoc_candidate",
            "",
            "review",
            "O2の解離を伴う可能性がある励起。いまはEEDFの損失のみ",
        )
    return "excitation_other", "", "none", text[:60]


def channel_id(target, category, threshold):
    r"""対応表に使うチャンネルIDを作る。

    Parameters
    ----------
    target : str
        標的名。
    category : str
        分類名。
    threshold : float
        しきい値 [eV]。

    Returns
    -------
    str
        ``標的->分類_しきい値`` の形式のID。
    """
    tag = category.upper().replace("+", "P").replace("_", "_")
    return f"{target}->{tag}_{threshold:07.3f}"


def build():
    r"""全過程を読み込み、対応表の行とプロファイル断片を作る。

    Returns
    -------
    tuple
        CSV用の行のリスト、プロファイル辞書、出典ハッシュの辞書。
    """
    rows, channels, hashes = [], [], {}
    co2_classify = classify_co2_ist if CO2_SET == "ist" else classify_co2
    for target, name in sources().items():
        path = SOURCE_DIR / name
        hashes[name] = sha256(path)
        with path.open(encoding="utf-8") as fp:
            collisions = parser.parse(fp)
        for c in collisions:
            classify = co2_classify if target == "CO2" else classify_other
            args = (c,) if target == "CO2" else (target, c)
            category, equation, status, note = classify(*args)
            threshold = float(c.get("threshold", 0.0))
            include = bool(equation) and status != "none" or category in (
                "attachment_2body",
                "ion_O2+",
            )
            include = include and bool(equation)
            cid = channel_id(target, category, threshold) if include else ""
            rows.append(
                dict(
                    file=name,
                    target=target,
                    kind=c["kind"],
                    product=c.get("product") or "",
                    threshold_eV=f"{threshold:.6g}",
                    category=category,
                    source_comment=comment_line(c),
                    in_profile="yes" if include else "no",
                    channel_id=cid,
                    proposed_equation=equation,
                    status=status,
                    note=note,
                )
            )
            if include:
                entry = dict(
                    id=cid,
                    target=target,
                    kind=c["kind"],
                    product=c.get("product"),
                    threshold_eV=threshold,
                    label=f"{target} {category}",
                )
                if category == "attachment_3body":
                    entry["density_normalization"] = dict(
                        species="O2",
                        power=1,
                        density_unit="cm^-3",
                        coefficient_scale=1e-06,
                        virtual_target="__O2_three_body__",
                    )
                channels.append(entry)
    profile = dict(
        cross_sections=[f"{PLANNED_DIR}/{n}" for n in sources().values()],
        targets={t: t for t in sources()},
        channels=channels,
    )
    return rows, profile, hashes


def check(profile):
    r"""対応表でBOLOSを初期化し、純CO2と混合気で収束を確認する。

    Parameters
    ----------
    profile : dict
        ``build`` が作ったプロファイル断片。
    """
    from plasglo.electrons.mixture import ElectronSwarm

    paths = [SOURCE_DIR / n for n in sources().values()]

    def new_swarm():
        r"""履歴のない `ElectronSwarm` を作る。

        Returns
        -------
        ElectronSwarm
            新規の電子群計算。格子の再調整を持たない旧版との比較でも
            E/Nごとの値が揃うよう、毎回作り直して使う。
        """
        return ElectronSwarm(
            paths, profile["targets"], profile["channels"], cells=400
        )

    swarm = new_swarm()
    print(f"initialised: {len(swarm.channel_ids)} channels")
    n_gas = 101325 / (1.380649e-23 * 300)
    cases = {
        "pure CO2": {"CO2": 1.0, "CO": 0.0, "O2": 0.0},
        "CO2/CO/O2 = 70/20/10": {"CO2": 0.7, "CO": 0.2, "O2": 0.1},
    }

    def total(state, *prefixes):
        r"""指定した接頭辞のチャンネルの反応係数を合計する。

        Parameters
        ----------
        state : dict
            ``ElectronSwarm.solve`` の戻り値。
        *prefixes : str
            合計対象のチャンネルIDの接頭辞。

        Returns
        -------
        float
            反応係数の合計 [m^3/s]。
        """
        return sum(
            state["rates"][i]
            for i in swarm.channel_ids
            if i.startswith(prefixes)
        )

    for label, frac in cases.items():
        print(f"== {label}")
        for td in (30.0, 100.0, 300.0):
            dens = {k: v * n_gas for k, v in frac.items()}
            swarm = new_swarm()
            state = swarm.solve(td, dens, 300.0)
            print(
                f"  E/N={td:5.0f} Td  <e>={state['energy_eV']:.3f} eV"
                f"  it={state['iterations']}  err={state['error']:.1e}"
                f"  k_diss={total(state, 'CO2->DISS', 'CO2->DISEXC'):.3e}"
                f"  k_ion={total(state, 'CO2->ION'):.3e}"
                f"  k_att={total(state, 'CO2->ATTACHMENT'):.3e}  [m3/s]"
            )


def main():
    r"""対応表を作って書き出し、``--check`` ならBOLOSの確認も行う。"""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    ap.add_argument(
        "--source-dir",
        help="断面積ファイルのフォルダ（既定：環境変数 CO2_XS_DIR）",
    )
    ap.add_argument(
        "--co2-set",
        choices=sorted(CO2_FILES),
        default="ist",
        help="CO2 の出典（既定：ist = IST-Lisbon の更新版）",
    )
    args = ap.parse_args()
    global CO2_SET
    CO2_SET = args.co2_set
    set_source_dir(args.source_dir)
    rows, profile, hashes = build()
    with (HERE / f"co2_channel_map_{CO2_SET}.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    payload = dict(source_sha256=hashes, profile=profile)
    (HERE / f"co2_electrons_profile_{CO2_SET}.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    by = {}
    for r in rows:
        by.setdefault((r["target"], r["category"], r["in_profile"]), 0)
        by[(r["target"], r["category"], r["in_profile"])] += 1
    for key, n in sorted(by.items()):
        print(f"{key[0]:4s} {key[1]:28s} profile={key[2]:3s} n={n}")
    print("channels in profile:", len(profile["channels"]))
    if args.check:
        check(profile)


if __name__ == "__main__":
    main()
