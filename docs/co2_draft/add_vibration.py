r"""CO2 の非対称伸縮振動の梯子（v1..v21）を、変換済みの機構に追加する。

Koelman et al., *Plasma Process. Polym.* 14, 1600155 (2017) の表 A3・A5 と
スケーリング則（式14〜19）に基づく、簡略版である（Kozak & Bogaerts 2014 由来）。
入れるもの:

- 電子衝突の振動励起 ``E + CO2 -> E + CO2v1``（IST-Lisbon の v=001 断面積）。
- VT 緩和 ``CO2vn + M -> CO2v(n-1) + M``（式14、v1 は V2a-c の和）。逆反応は
  詳細釣り合いで入れる。
- VV 交換 ``CO2vn + CO2vj -> CO2v(n-1) + CO2v(j+1)``（式18、基準は V8）。
  逆反応は詳細釣り合い。
- 振動による解離の促進 ``CO2vn + M -> CO + O + M``（N1、alpha=0.8）と
  ``CO2vn + O -> CO + O2``（N2、alpha=0.5）。活性化エネルギーを
  ``alpha * E_n`` だけ下げる（0 で打ち切る）。

入れないもの（簡略化）: 対称モード（va-vd）の準位、CO・O2 の振動、電子衝突に
よる準位間の遷移（v_n -> v_m）と多量子励起、超弾性衝突、O 原子による VT 緩和。

使い方::

    python docs/co2_draft/add_vibration.py outputs/co2_mechanism \
        outputs/co2_mechanism_vib
"""

import csv
import json
import math
import sys
from pathlib import Path

N_LEVELS = 21
EPS1 = 0.291  # v1 のエネルギー [eV]
XE = 5.25e-3  # 非対称モードの非調和性
KB_EV = 8.617333262e-05  # [eV/K]
KB_J = 1.380649e-23
HBAR_J = 1.054571817e-34
EV_J = 1.602176634e-19
AMU = 1.66053907e-27
ALPHA = 17.5 / 3.94 * 1e10  # 指数型反発ポテンシャルの係数 [1/m]（CO2）
MASS = 44.0 * AMU / 2.0  # CO2-CO2 の換算質量 [kg]
V_COLLIDERS = {"CO2": 1.0, "CO": 0.3, "O2": 0.4}  # 表 A5 の注 b
TEMPERATURE_K = 300.0


def level(n):
    r"""準位の名前を返す。

    Parameters
    ----------
    n : int
        振動量子数（0 は基底状態）。

    Returns
    -------
    str
        化学種名。

    Examples
    --------
    >>> level(0), level(3)
    ('CO2', 'CO2v3')
    """
    return "CO2" if n == 0 else f"CO2v{n}"


def energy(n):
    r"""非調和振動子の準位エネルギー [eV]。

    Parameters
    ----------
    n : int
        振動量子数。

    Returns
    -------
    float
        ``EPS1 * n * (1 - XE * (n - 1))``。

    Examples
    --------
    >>> energy(1)
    0.291
    >>> 5.3 < energy(21) < 5.6
    True
    """
    return EPS1 * n * (1.0 - XE * (n - 1)) if n else 0.0


def zeta(n):
    r"""式15の係数 Z_n。

    Parameters
    ----------
    n : int
        振動量子数。

    Returns
    -------
    float
        ``n (1 - xe) / (1 - n xe)``。

    Examples
    --------
    >>> zeta(1)
    1.0
    """
    return n * (1.0 - XE) / (1.0 - n * XE)


def gamma_coefficient(delta_ev):
    r"""式17の γ を、T^(-1/2) の係数として返す。

    Parameters
    ----------
    delta_ev : float
        反応でのエネルギー差の大きさ [eV]。

    Returns
    -------
    float
        ``gamma = 係数 / sqrt(T)`` の係数 [K^0.5]。

    Examples
    --------
    >>> gamma_coefficient(0.0)
    0.0
    >>> gamma_coefficient(0.3) > gamma_coefficient(0.1)
    True
    """
    nu = abs(delta_ev) * EV_J / HBAR_J
    return math.pi * nu / ALPHA * math.sqrt(MASS / (2.0 * KB_J))


def f_factor(coefficient):
    r"""式16の F(γ) を、tgas の式（文字列）にする。

    Parameters
    ----------
    coefficient : float
        ``gamma_coefficient`` の値。

    Returns
    -------
    str
        ``F`` を表す式。``coefficient`` が 0 なら ``1``。

    Examples
    --------
    >>> f_factor(0.0)
    '1'
    """
    if coefficient == 0.0:
        return "1"
    g = f"({coefficient:.6g}/sqrt(tgas))"
    return f"(0.5*(3-exp(-2*{g}/3))*exp(-2*{g}/3))"


def vt_base():
    r"""V2a-c の和（v1 -> 基底、衝突相手 CO2 基準）を表す式を返す。

    Returns
    -------
    str
        ``tgas`` の式 [m^3/s]。
    """
    parts = ((4.25e-7, 407, 824), (8.57e-7, 404, 1096), (1.43e-11, 252, 685))
    return "(" + " + ".join(
        f"{a:g}*exp(-{b}*tgas**(-1/3)+{c}*tgas**(-2/3))"
        for a, b, c in parts
    ) + ")"


def reactions():
    r"""追加する反応を、CSV の行（辞書）のリストにして返す。

    Returns
    -------
    list of dict
        追加する反応。

    Examples
    --------
    >>> rows = reactions()
    >>> rows[0]["equation"]
    'E + CO2 -> E + CO2v1'
    >>> len({r["id"] for r in rows}) == len(rows)
    True
    """
    rows = []

    def add(rid, equation, kind, value, note):
        r"""反応を1本、``rows`` に足す。

        Parameters
        ----------
        rid : str
            反応ID。
        equation : str
            反応式。
        kind : str
            速度の種類（``electron`` または ``expression``）。
        value : str
            チャンネルIDまたは ``tgas`` の式。
        note : str
            出典の注記。
        """
        rows.append(dict(
            id=rid, equation=equation, kind=kind, value=value, parameters="",
            length_unit="m", source_file="Koelman 2017 (Table A3, A5)",
            line="", reference=note))

    add("V000", "E + CO2 -> E + CO2v1", "electron",
        "CO2->EXC_0.291_CO2(v=001)", "IST-Lisbon, Lowke-Phelps-Irwin 1973")
    f1 = f_factor(gamma_coefficient(EPS1))
    for n in range(1, N_LEVELS + 1):
        gap = energy(n) - energy(n - 1)
        scale = (f"{zeta(n):.6g}*{f1}/"
                 f"{f_factor(gamma_coefficient(gap))}") if n > 1 else "1"
        down = f"{vt_base()}*{scale}"
        add(f"V{n:03d}D", f"{level(n)} + MV -> {level(n - 1)} + MV",
            "expression", down, "VT, 式14（V2a-c の和）")
        up = f"{down}*exp(-{gap:.6f}/({KB_EV:.9g}*tgas))"
        add(f"V{n:03d}U", f"{level(n - 1)} + MV -> {level(n)} + MV",
            "expression", up, "VT 逆反応（詳細釣り合い）")
    k_vv = "1.32e-16*(tgas/300)**0.5*exp(-250/tgas)"
    for j in range(0, N_LEVELS - 1):
        for n in range(j + 2, N_LEVELS + 1):
            m = j + 1
            delta = energy(n) + energy(j) - energy(n - 1) - energy(m)
            fac = f_factor(gamma_coefficient(delta))
            forward = f"{k_vv}*{zeta(n) * zeta(m):.6g}*{fac}"
            products = (f"2{level(n - 1)}" if n - 1 == m
                        else f"{level(n - 1)} + {level(m)}")
            add(f"W{n:02d}{j:02d}F",
                f"{level(n)} + {level(j)} -> {products}", "expression",
                forward, "VV, 式18（基準 V8）")
            backward = f"{forward}*exp(-{delta:.6f}/({KB_EV:.9g}*tgas))"
            add(f"W{n:02d}{j:02d}B",
                f"{products} -> {level(n)} + {level(j)}", "expression",
                backward, "VV 逆反応（詳細釣り合い）")
    for n in range(1, N_LEVELS + 1):
        ea1 = max(49000.0 - 0.8 * energy(n) / KB_EV, 0.0)
        ea2 = max(26500.0 - 0.5 * energy(n) / KB_EV, 0.0)
        add(f"D{n:03d}A", f"{level(n)} + M -> CO + O + M", "expression",
            f"1.81e-16*exp(-{ea1:.1f}/tgas)", "N1（alpha=0.8）")
        add(f"D{n:03d}B", f"{level(n)} + O -> CO + O2", "expression",
            f"2.8e-17*exp(-{ea2:.1f}/tgas)", "N2（alpha=0.5）")
    return rows


def main():
    r"""変換済みの機構を読み、振動を追加した機構とケースを書く。"""
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    dst.mkdir(parents=True, exist_ok=True)
    case = json.loads((src / "case_co2.json").read_text("utf-8-sig"))
    mech = case["mechanism"]
    mech["species"] += [
        dict(name=level(n), composition={"C": 1, "O": 2}, charge=0,
             excited=True) for n in range(1, N_LEVELS + 1)]
    mech["third_bodies"]["MV"] = dict(V_COLLIDERS)
    mech["reactions"] = "mechanism_co2_vib.csv"
    channel = dict(
        id="CO2->EXC_0.291_CO2(v=001)", target="CO2", kind="EXCITATION",
        product="CO2(v=001)", threshold_eV=0.291,
        label="e + CO2(X) -> e + CO2(v=001)")
    electrons = case["electrons"]
    electrons["profile"]["channels"].append(channel)
    case["name"] = "co2_loki_vib_draft"
    (dst / "case_co2_vib.json").write_text(
        json.dumps(case, ensure_ascii=False, indent=1), encoding="utf-8")
    with (src / "mechanism_co2.csv").open(encoding="utf-8-sig") as fp:
        rows = list(csv.DictReader(fp))
        header = list(rows[0])
    rows += reactions()
    with (dst / "mechanism_co2_vib.csv").open(
            "w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    print(f"反応 {len(rows)} 本（振動 {len(reactions())} 本）、"
          f"化学種 {len(mech['species'])} 種 -> {dst}")


if __name__ == "__main__":
    main()
