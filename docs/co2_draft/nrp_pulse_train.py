r"""CO2 中の NRP（ナノ秒繰り返しパルス）を、0次元のパルス列として計算する。

1サイクルを、パルス中（E/N あり）と、パルス間（E/N = 0）の2回の
``simulate`` に分け、状態を引き継ぐ。各サイクルで、CO2 の生成・消滅を、
``analyze_pathways`` で反応ごとに分解し、JSON Lines に記録する。

仮定（実験条件ではなく、例示）：

- E/N 一定の矩形パルス。電子密度は、各パルスの初めに ``--seed-ne`` 以上に
  保証する（不足分は、電荷が中立になるよう、CO2^+ とともに加える）。
  チャネル内のプラズマを、0次元の一様な気体として扱う。
- 閉じた定容、気体温度一定、壁・流れ・拡散なし。気体加熱なし。
- 反応機構は、``loki_to_plasglo.py`` が作ったケース（LoKI-GM 由来。GPL-3.0
  のため、リポジトリには含めない）。

使い方::

    python docs/co2_draft/nrp_pulse_train.py \
        --case outputs/co2_mechanism/case_co2.json --cycles 20
"""

import argparse
import json
import logging
import math
import sys
import time
from pathlib import Path

import numpy as np

logging.disable(logging.CRITICAL)
MODEL = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(MODEL))
from plasglo import analyze_pathways, load_case, simulate  # noqa: E402

# 射影積分でまとめて進める「遅い種」。1サイクルの相対変化が小さく、累積する。
SLOW_SPECIES = ("CO2", "CO", "O2", "O3")
ATOL_DENSITY = 1e5  # BDF の絶対許容誤差 [m^-3]（密度は 1e19〜1e25 m^-3）
RTOL = 1e-6  # BDF の相対許容誤差の既定値
# パルス間を分ける区間の境界 [ns]（パルス終了からの時間）。
GAP_EDGES_NS = (0.0, 5.0, 50.0, 500.0, 5e3, 5e4)


def final_state(result):
    r"""計算の最終時刻の密度を、種名から密度への辞書にする。

    Parameters
    ----------
    result : plasglo.Result
        ``simulate`` の結果。

    Returns
    -------
    dict of str to float
        各種の数密度 [m^-3]。丸め誤差の負値は 0 にする。
    """
    names = list(result.mechanism.species)
    return {
        n: max(float(v), 0.0) for n, v in zip(names, result.densities[:, -1])
    }


def ensure_seed(state, seed_ne):
    r"""電子密度が ``seed_ne`` 未満なら、CO2 を電離して、不足分を作る。

    CO2 を減らして、CO2^+ と電子を同じだけ増やすので、C・O の総量も
    電荷も変わらない（以前は CO2^+ だけを足しており、1サイクルあたり
    約 1e20 m^-3 の物質が増えていた）。

    Parameters
    ----------
    state : dict of str to float
        各種の数密度 [m^-3]。
    seed_ne : float
        パルスの初めに保証する電子密度 [m^-3]。

    Returns
    -------
    dict of str to float
        更新した密度（新しい辞書）。

    Examples
    --------
    >>> s = ensure_seed({"E": 1e10, "CO2^+": 1e10, "CO2": 1e13}, 1e12)
    >>> s["E"] == 1e12 and s["CO2^+"] == 1e12
    True
    >>> round(1e13 - s["CO2"])
    990000000000
    """
    out = dict(state)
    missing = seed_ne - out.get("E", 0.0)
    if missing > 0:
        out["E"] = seed_ne
        out["CO2^+"] = out.get("CO2^+", 0.0) + missing
        out["CO2"] = out["CO2"] - missing
    return out


def neutralize(state, charges):
    r"""正味の電荷を、電子の密度の調整で、ゼロに戻す。

    BDF の丸め誤差で、位相の終わりの電荷が、わずかに崩れる。ケースの検査
    （相対 1e-10）を超えると、次の位相が始められないので、電子で調整する。

    Parameters
    ----------
    state : dict of str to float
        各種の数密度 [m^-3]。
    charges : dict of str to int
        種名から電荷への辞書（電子は -1）。

    Returns
    -------
    dict of str to float
        電荷が釣り合った密度（新しい辞書）。

    Examples
    --------
    >>> ions = {"E": 1.0, "CO2^+": 1.0 + 1e-6}
    >>> s = neutralize(ions, {"E": -1, "CO2^+": 1})
    >>> round(s["E"] - s["CO2^+"], 9)
    0.0
    """
    out = dict(state)
    net = sum(charges[n] * v for n, v in out.items() if n in charges)
    if net > 0:
        out["E"] += net
    elif net < 0:
        removed = min(-net, out["E"])
        out["E"] -= removed
        out["CO2^+"] = out.get("CO2^+", 0.0) + (-net - removed)
    return out


def close_balance(state, compositions, targets):
    r"""C と O の総量が目標値になるよう、CO2 と O2 で、わずかに調整する。

    射影積分（ジャンプ）では、速い種の変化を無視するため、元素の総量が
    わずかにずれる。ジャンプのたびに、初期の総量へ戻す。

    Parameters
    ----------
    state : dict of str to float
        各種の数密度 [m^-3]。
    compositions : dict of str to dict
        種名から元素組成（``{"C": 1, "O": 2}`` など）への辞書。
    targets : dict of str to float
        元素ごとの総量の目標値 [m^-3]。

    Returns
    -------
    dict of str to float
        総量を合わせた密度（新しい辞書）。

    Examples
    --------
    >>> comp = {"CO2": {"C": 1, "O": 2}, "O2": {"O": 2}}
    >>> now, goal = {"CO2": 9.9, "O2": 0.1}, {"C": 10.0, "O": 20.0}
    >>> s = close_balance(now, comp, goal)
    >>> round(s["CO2"], 9), round(s["O2"], 9)
    (10.0, 0.0)
    """
    out = dict(state)

    def total(element):
        r"""現在の状態で、元素の総量を数える。

        Parameters
        ----------
        element : str
            元素記号（``C`` または ``O``）。

        Returns
        -------
        float
            元素の総量 [m^-3]。
        """
        return sum(
            compositions[n].get(element, 0) * v
            for n, v in out.items()
            if n in compositions
        )

    delta_c = targets["C"] - total("C")
    out["CO2"] += delta_c
    delta_o = targets["O"] - total("O")
    out["O2"] += delta_o / 2
    return out


def project(delta, before, cycle, skipped):
    r"""1サイクルの変化量が、サイクル数の冪で減るとして、進める量を求める。

    ``before`` が与えられ、``delta`` と同符号なら、2点から冪指数 alpha を
    求め、``delta * sum_j ((cycle + j) / cycle) ** (-alpha)`` を返す。
    符号が違う（0 をまたぐ）か、サイクル番号が合わなければ、直線で外挿する。
    ``before`` がなければ、一定として外挿する。

    Parameters
    ----------
    delta : float
        いまのサイクルでの、1サイクルの変化量 [m^-3]。
    before : tuple or None
        (サイクル番号, その時の変化量)。なければ None。
    cycle : int
        いまのサイクル番号（経過サイクル数）。
    skipped : int
        進めるサイクル数。

    Returns
    -------
    float
        進める量 [m^-3]。

    Examples
    --------
    >>> round(project(2.0, None, 10, 3), 6)
    6.0
    >>> round(project(1.0, (5, 2.0), 10, 3), 6) < 3.0
    True
    >>> round(project(-1.0, (5, -1.0), 10, 3), 6)
    -3.0
    """
    if before is None:
        return skipped * delta
    cycle0, delta0 = before
    if delta * delta0 > 0 and cycle > cycle0 > 0:
        alpha = -math.log(delta / delta0) / math.log(cycle / cycle0)
        return delta * sum(
            ((cycle + j) / cycle) ** (-alpha) for j in range(1, skipped + 1)
        )
    slope = (delta - delta0) / max(cycle - cycle0, 1)
    return skipped * delta + slope * skipped * (skipped + 1) / 2


def run_phase(case, state, td, end_ns, points, update_ns, rtol=RTOL):
    r"""1つの位相（パルス中、またはパルス間）を計算する。

    Parameters
    ----------
    case : plasglo.Case
        基準のケース。
    state : dict of str to float
        位相の初めの密度 [m^-3]。
    td : float
        換算電界 [Td]。0 ならパルス間。
    end_ns : float
        位相の長さ [ns]。
    points : int
        出力点数（線形の等間隔）。
    update_ns : float
        BOLOS の更新間隔 [ns]。
    rtol : float, optional
        BDF の相対許容誤差。

    Returns
    -------
    plasglo.Result
        この位相の計算結果。
    """
    # ケースの許容誤差（atol 1e-4 m^-3、rtol 1e-8）は、密度が 1e19 m^-3 を
    # 超える計算では厳しすぎ、BDF の丸め程度の負値でも結合更新が止まる。
    phase = case.with_updates(
        initial_densities_m3=state,
        field={"td": td, "duration_ns": None},
        solver={
            "end_ns": end_ns,
            "points": points,
            "rtol": rtol,
            "atol_density_m3": ATOL_DENSITY,
        },
        electrons={"update_ns": update_ns, "update_mode": "coupled"},
    )
    return simulate(phase)


def contributions(result, species="CO2"):
    r"""1つの位相での、化学種への正味の寄与を、反応ごとに求める。

    Parameters
    ----------
    result : plasglo.Result
        1つの位相の結果。
    species : str, optional
        解析する化学種。

    Returns
    -------
    tuple
        反応IDから正味の積算寄与 [m^-3]（生成を正）への辞書、反応IDから
        反応式への辞書、位相の初めと終わりの密度差 [m^-3]。
    """
    paths = analyze_pathways(result, species)
    values = dict(zip(paths.reaction_ids, paths.integrated_contributions_m3))
    equations = dict(zip(paths.reaction_ids, paths.equations))
    return values, equations, float(paths.density_change_m3)


def summarize(parts, top=6):
    r"""複数の位相の寄与を足し合わせ、生成側・消滅側の上位をまとめる。

    Parameters
    ----------
    parts : list of tuple
        ``contributions`` の戻り値のリスト。
    top : int, optional
        生成側・消滅側それぞれの上位の数。

    Returns
    -------
    dict
        密度の差 ``change``、生成側・消滅側の上位（反応ID、反応式、積算量）と
        それぞれの合計。消滅側の積算量も、正の大きさで表す。

    Examples
    --------
    >>> eqs = {"R1": "A -> B", "R2": "B -> A"}
    >>> a = ({"R1": 2.0, "R2": -1.0}, eqs, 1.0)
    >>> s = summarize([a, a])
    >>> s["change"], s["production"][0]["id"], s["loss"][0]["m3"]
    (2.0, 'R1', 2.0)
    """
    total, equations, change = {}, {}, 0.0
    for values, eqs, delta in parts:
        for rid, v in values.items():
            total[rid] = total.get(rid, 0.0) + float(v)
        equations.update(eqs)
        change += delta

    def pick(sign):
        r"""生成側（sign = 1）または消滅側（sign = -1）の反応を並べる。

        Parameters
        ----------
        sign : int
            1 なら生成側、-1 なら消滅側。

        Returns
        -------
        list of dict
            積算量の降順の反応（積算量は、どちらも正の大きさ）。
        """
        items = sorted(
            (rid for rid, v in total.items() if sign * v > 0),
            key=lambda rid: -sign * total[rid],
        )
        return [
            dict(id=rid, eq=equations[rid], m3=sign * total[rid])
            for rid in items
        ]

    production, loss = pick(1), pick(-1)
    return dict(
        change=change,
        production=production[:top],
        loss=loss[:top],
        total_production=float(sum(i["m3"] for i in production)),
        total_loss=float(sum(i["m3"] for i in loss)),
    )


# CO(a3Pi) の反応の係数を、Cenian et al. 1994 の値（Loenders らの SI Table S1。
# 単位 m^3/s）に置き換える感度計算用。LoKI の既定値は、その 0.33〜0.65 倍。
CENIAN = {
    "CO(a3Pi) + O2 -> CO + O2": "5.0e-17",
    "CO(a3Pi) + O2 -> CO + 2O": "5.0e-17",
    "CO(a3Pi) + O2 -> CO2 + O": "3.0e-17",
    "CO(a3Pi) + CO -> 2CO": "1.4e-16",
    "CO(a3Pi) + CO -> CO2 + C": "1.4e-18",
    "CO(a3Pi) + CO2 -> CO + CO2": "1.5e-17",
    "CO(a3Pi) + CO2 -> 2CO + O": "1.5e-17",
}


def apply_variant(case, variant, out):
    r"""必要なら、CO(a3Pi) の係数を置き換えた反応表を作り、ケースへ反映する。

    Parameters
    ----------
    case : plasglo.Case
        基準のケース。
    variant : str
        ``loki``（そのまま）または ``cenian``（CO(a3Pi) の7反応を置き換える）。
    out : pathlib.Path
        変更した反応表の書き出し先フォルダ。

    Returns
    -------
    plasglo.Case
        反映したケース。
    """
    if variant == "loki":
        return case
    import csv

    mechanism = dict(case.data["mechanism"])
    # 反応表のパスは、ケースJSONの場所からの相対のことがある。
    source = Path(case.source).parent / mechanism["reactions"]
    with open(source, newline="", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        fields, rows = reader.fieldnames, list(reader)
    changed = 0
    for row in rows:
        if row["equation"] in CENIAN:
            row["value"] = CENIAN[row["equation"]]
            row["kind"] = "expression"
            changed += 1
    if changed != len(CENIAN):
        raise SystemExit(f"CO(a3Pi) の反応を {changed}/{len(CENIAN)} 件しか置換できない")
    path = out / "mechanism_variant.csv"
    with path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    mechanism["reactions"] = str(path.resolve())
    return case.with_updates(mechanism=mechanism)


def main():
    r"""パルス列を計算し、サイクルごとの記録を表示・保存する。"""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--case", required=True, help="変換したケースJSON")
    ap.add_argument("--td", type=float, default=100.0, help="E/N [Td]")
    ap.add_argument("--pulse-ns", type=float, default=20.0)
    ap.add_argument("--period-us", type=float, default=100.0)
    ap.add_argument(
        "--cycles",
        type=int,
        default=10,
        help="直接計算するサイクル数（--until を指定しないとき）",
    )
    ap.add_argument(
        "--until",
        type=int,
        default=0,
        help="経過サイクル数（ジャンプで進めた分を含む）がこの値に達する"
        "まで計算する。0 なら --cycles 回の直接計算",
    )
    ap.add_argument("--seed-ne", type=float, default=1e19)
    ap.add_argument("--x0", type=float, default=0.0, help="初期の分解率")
    ap.add_argument("--out", default="outputs/co2_nrp")
    ap.add_argument(
        "--variant",
        choices=("loki", "cenian"),
        default="loki",
        help="CO(a3Pi) の係数：loki（既定）または cenian",
    )
    ap.add_argument(
        "--pulse-update-ns",
        type=float,
        default=0.5,
        help="パルス中の BOLOS の更新間隔 [ns]。EEDF は組成だけで決まり、"
        "20 ns で組成はほとんど変わらないので、5 ns でも分解率の増分の"
        "差は 0.2%% 程度（0.5 ns 比）",
    )
    ap.add_argument(
        "--bolos-cache-tol",
        type=float,
        default=0.0,
        help="0 より大きいと、BOLOS の結果を、組成が近ければ使い回す。"
        "標的の割合を量子化する幅（例 1e-3）。0 で無効",
    )
    ap.add_argument(
        "--jump",
        type=int,
        default=0,
        help="射影積分：--block 回の直接計算のあと、遅い種（CO2、CO、O2、"
        "O3）を、最後のサイクルの変化量で、この値のサイクル分まとめて"
        "進める。0 または 1 で無効",
    )
    ap.add_argument(
        "--block",
        type=int,
        default=3,
        help="ジャンプの間に直接計算するサイクル数（速い種が準定常に"
        "戻る時間）",
    )
    ap.add_argument(
        "--jump-after",
        type=int,
        default=0,
        help="経過サイクル数がこの値以上になってから、ジャンプを始める",
    )
    ap.add_argument(
        "--resume",
        action="store_true",
        help="出力フォルダの最後のサイクルの状態から再開し、記録に追記する",
    )
    ap.add_argument(
        "--stop-dx",
        type=float,
        default=0.0,
        help="1サイクルの分解率の増分（割合）がこれ未満になったら止める",
    )
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    case = apply_variant(load_case(args.case), args.variant, out)
    n0 = case.data["initial_densities_m3"]["CO2"]
    species = case.data["mechanism"]["species"]
    names = [s["name"] for s in species]
    charges = {s["name"]: s["charge"] for s in species}
    compositions = {s["name"]: s["composition"] for s in species}
    targets = {"C": n0, "O": 2 * n0}  # C・O の総量（初期は純 CO2 と同じ）
    state = {"CO2": n0 * (1 - args.x0)}
    if args.x0 > 0:
        state["CO"] = n0 * args.x0
        state["O2"] = n0 * args.x0 / 2
    state = {n: state.get(n, 0.0) for n in names}
    cache = None
    if args.bolos_cache_tol > 0:
        import bolos_cache

        cache = bolos_cache.install(args.bolos_cache_tol)
    inter_ns = args.period_us * 1e3 - args.pulse_ns
    path = out / "cycles.jsonl"
    done = 0  # 経過サイクル数（ジャンプで進めた分を含む）
    if args.resume and path.exists():
        with path.open(encoding="utf-8-sig") as fp:
            last = [json.loads(line) for line in fp if line.strip()][-1]
        state = {n: last["state_end_gap"].get(n, 0.0) for n in names}
        done = last["cycle"]
    log = path.open("a" if args.resume else "w", encoding="utf-8")
    print(f"E/N={args.td} Td, pulse={args.pulse_ns} ns, "
          f"period={args.period_us} us, seed ne={args.seed_ne:.1e}"
          f"{'  [再開：サイクル' + str(done + 1) + 'から]' if done else ''}"
          f"{'  [射影積分：' + str(args.block) + '回の直接計算ごとに '
             + str(args.jump) + ' サイクル進める]' if args.jump > 1 else ''}",
          flush=True)
    print(f"{'cycle':>6s} {'X[%]':>9s} {'dX[%]':>9s} {'ne_end':>9s} "
          f"{'O':>9s} {'CO':>9s} {'O2':>9s} {'time[s]':>8s}", flush=True)
    edges = [e for e in GAP_EDGES_NS if e < inter_ns] + [inter_ns]
    direct = in_block = 0
    previous_block = None  # 前のブロックの (サイクル番号, 変化量)
    last_delta = None  # 直前の直接計算サイクルの変化量
    while (done < args.until) if args.until > 0 else (direct < args.cycles):
        done += 1
        direct += 1
        t0 = time.time()
        start = dict(state)
        state = ensure_seed(state, args.seed_ne)
        pulse = run_phase(
            case, state, args.td, args.pulse_ns, 81, args.pulse_update_ns
        )
        # 電荷は、位相の間で中立に戻す（BDF の丸め誤差が、ケースの検査の
        # 許容 1e-10 を超えると、次の位相が始まらない）。
        mid = neutralize(final_state(pulse), charges)
        # パルス間は、失活が ns で終わる一方、全体は 100 µs 級なので、
        # 刻みを段階的に粗くした区間に分け、反応ごとの寄与を足し合わせる。
        gap_parts, state = [], mid
        for lo, hi in zip(edges[:-1], edges[1:]):
            result = run_phase(case, state, 0.0, hi - lo, 21, (hi - lo) / 4)
            gap_parts.append(contributions(result))
            state = neutralize(final_state(result), charges)
        a = summarize([contributions(pulse)])
        b = summarize(gap_parts)
        # 断面積のない中性種（O 原子など）の割合。EEDF に入らない分の目安。
        unsupported = max(
            h["unsupported_neutral_fraction"] for h in pulse.swarm_history
        )
        x = 1 - state["CO2"] / n0
        dx = x - (1 - start["CO2"] / n0)
        record = dict(cycle=done, pulse=a, gap=b, state_end_pulse=mid,
                      state_end_gap=state, unsupported_max=unsupported,
                      x_end=x, dx=dx)
        log.write(json.dumps(record, ensure_ascii=False) + "\n")
        log.flush()
        print(f"{done:6d} {100 * x:9.5f} {100 * dx:9.5f} {mid['E']:9.2e} "
              f"{state['O']:9.2e} {state['CO']:9.2e} {state['O2']:9.2e} "
              f"{time.time() - t0:8.1f}"
              + (f"  [BOLOS 使い回し {cache['hit']} / 計算 {cache['miss']}]"
                 if cache else ""), flush=True)
        if args.stop_dx > 0 and abs(dx) < args.stop_dx:
            print(f"増分 {100 * dx:.5f}% が {100 * args.stop_dx:g}% 未満になり、"
                  "停止", flush=True)
            break
        in_block += 1
        delta = {s: state[s] - start[s] for s in SLOW_SPECIES}
        previous_cycle, last_delta = last_delta, delta
        ready = in_block >= args.block and done >= args.jump_after
        if args.jump > 1 and ready:
            # 射影積分：遅い種を、最後のサイクルの変化量で、まとめて進める。
            # 速い種（O 原子、イオンなど）は、そのまま次の直接計算で、準定常
            # に戻す。C・O の総量は、初期値に戻し、電荷も中立にする。
            skipped = args.jump - 1
            # 1サイクルの変化量は、サイクルとともに下がる（分解率の増分は、
            # サイクル数の冪で減る）。一定として進めると過大になるので、
            # 前のブロックの同じ位相（速い種の緩和が同じ程度）の変化量との
            # 比から冪指数を求め、その冪で減るとして、まとめて進める。
            # 最初のブロックだけは、直前のサイクルとの差から、直線で外挿する。
            for s in SLOW_SPECIES:
                before = None if previous_block is None else (
                    previous_block[0], previous_block[1][s])
                if before is None and previous_cycle is not None:
                    before = (done - 1, previous_cycle[s])
                state[s] = max(
                    state[s] + project(delta[s], before, done, skipped), 0.0
                )
            previous_block = (done, delta)
            state = neutralize(
                close_balance(state, compositions, targets), charges
            )
            done += args.jump - 1
            x = 1 - state["CO2"] / n0
            log.write(json.dumps(
                dict(cycle=done, jump=args.jump - 1, state_end_gap=state,
                     x_end=x, dx=dx), ensure_ascii=False) + "\n")
            log.flush()
            print(f"{done:6d} {100 * x:9.5f} {100 * dx:9.5f}   "
                  f"<- {args.jump - 1} サイクル分をジャンプ", flush=True)
            in_block = 0
    log.close()


if __name__ == "__main__":
    main()
