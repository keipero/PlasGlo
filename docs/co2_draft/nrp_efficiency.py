r"""E/N と CO2 解離のエネルギー効率を、EEDF から求める（化学は入れない）。

純 CO2 と、``nrp_pulse_train.py`` の飽和までの各組成について、E/N を振って
BOLOS で EEDF を解き、次を出す。

- 電子1個あたりの入力電力 ``P = mu*N * (E/N)^2``（気体分子1個あたり）。
- CO2 の一次解離（CO+O(1D)、CO(a3Pi)+O、CO+O^- 付着）の速度。
- 解離1回あたりのエネルギーコスト ``P / 解離速度`` [eV] と、効率
  ``2.93 eV / コスト``（CO2 -> CO + 1/2 O2 の標準エンタルピー）。
- 入力電力の配分（振動、電子励起、解離、電離、CO・O2 での損失、弾性ほか）。

使い方::

    python docs/co2_draft/nrp_efficiency.py outputs/co2_nrp/loki_1e20 \
        --case outputs/co2_mechanism/case_co2.json
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

MODEL = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(MODEL))
from plasglo import load_case  # noqa: E402
from plasglo.electrons import mixture  # noqa: E402

plt.rcParams["font.family"] = ["Yu Gothic", "Meiryo", "sans-serif"]
ENTHALPY_EV = 2.93  # CO2 -> CO + 1/2 O2 [eV/分子]
GROUPS = ("振動(CO2)", "電子励起(CO2 e1)", "解離(CO2)", "電離",
          "CO・O2 の非弾性", "弾性ほか")


def classify(target, process):
    r"""非弾性過程を、電力配分の区分に分ける。

    Parameters
    ----------
    target : str
        標的の名前。
    process : object
        BOLOS の過程（``kind``、``threshold``、``product`` を持つ）。

    Returns
    -------
    tuple
        (区分名, CO2 の一次解離なら True)。

    Examples
    --------
    >>> class P:
    ...     kind, threshold, product = "EXCITATION", 0.083, "CO2(v=010)"
    >>> classify("CO2", P)
    ('振動(CO2)', False)
    """
    kind, product = str(process.kind), str(process.product)
    if target != "CO2":
        return GROUPS[4], False
    if kind == "IONIZATION":
        return GROUPS[3], False
    dissociation = product.startswith("CO(") and "O" in product[3:]
    if dissociation:
        return GROUPS[2], True
    if kind == "EXCITATION" and process.threshold < 3.0:
        return GROUPS[0], False
    return GROUPS[1], False


def efficiency_at(swarm, td, densities, gas_temperature=300.0):
    r"""1つの E/N、組成での、コストと電力配分を求める。

    Parameters
    ----------
    swarm : plasglo.electrons.mixture.ElectronSwarm
        電子群計算。
    td : float
        換算電界 [Td]。
    densities : dict of str to float
        標的の数密度 [m^-3]。
    gas_temperature : float, optional
        気体温度 [K]。

    Returns
    -------
    dict
        ``cost_eV``（解離1回あたり）、``efficiency``、``power``（電子1個、
        気体分子1個あたりの入力電力 [eV m^3/s]）、``mean_eV``、
        ``partition``（区分名から入力電力の割合）。
    """
    # E/N が大きく変わる走査では、エネルギー格子が最初の E/N に引きずられる
    # ので、毎回、格子を作り直して、初期値から解く。
    swarm.eedf, swarm.last_field_td = None, 0.0
    swarm.solve(td, densities, gas_temperature)
    b, f = swarm.bsolver, swarm.eedf
    en = td * mixture.solver.TOWNSEND
    power = b.mobility(f) * en**2  # [eV m^3 / s] / N
    parts = dict.fromkeys(GROUPS, 0.0)
    diss = 0.0
    for target, process in b.iter_inelastic():
        rate = b.rate(f, process, weighted=True)  # [m^3/s] / N
        group, primary = classify(target.name, process)
        parts[group] += process.threshold * rate
        if primary:
            diss += rate
    parts[GROUPS[5]] = max(power - sum(parts.values()), 0.0)
    cost = power / diss if diss > 0 else float("inf")
    return dict(
        cost_eV=cost, efficiency=ENTHALPY_EV / cost, power=power,
        mean_eV=float(b.mean_energy(f)),
        partition={k: v / power for k, v in parts.items()},
    )


def pick_compositions(records, targets):
    r"""分解率が目標に近い記録の状態を選ぶ。

    Parameters
    ----------
    records : list of dict
        ``cycles.jsonl`` の記録。
    targets : list of float
        目標の分解率（割合）。

    Returns
    -------
    list of tuple
        (実際の分解率, 状態) のリスト。

    Examples
    --------
    >>> r = [{"x_end": 0.1, "state_end_gap": {"CO2": 1}},
    ...      {"x_end": 0.5, "state_end_gap": {"CO2": 2}}]
    >>> pick_compositions(r, [0.4])[0][0]
    0.5
    """
    out = []
    for t in targets:
        r = min(records, key=lambda rec: abs(rec["x_end"] - t))
        out.append((r["x_end"], r["state_end_gap"]))
    return out


def main():
    r"""E/N と分解率の2次元の効率を計算し、図と JSON を書く。"""
    logging.disable(logging.WARNING)
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run")
    ap.add_argument("--case", required=True)
    ap.add_argument("--td-min", type=float, default=10.0)
    ap.add_argument("--td-max", type=float, default=500.0)
    ap.add_argument("--points", type=int, default=30)
    ap.add_argument("--x", type=float, nargs="*",
                    default=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.58])
    args = ap.parse_args()
    run = Path(args.run)
    case = load_case(args.case)
    s = case.electron_settings()
    p = s["profile_data"]
    swarm = mixture.ElectronSwarm(
        p["cross_sections"], p["targets"], p["channels"], cells=s["cells"]
    )
    n0 = case.data["initial_densities_m3"]["CO2"]
    with (run / "cycles.jsonl").open(encoding="utf-8-sig") as fp:
        records = [json.loads(line) for line in fp if line.strip()]
    for r in records:  # 直接計算の記録には x_end がないので、状態から求める
        r["x_end"] = 1 - r["state_end_gap"]["CO2"] / n0
    comps = [(0.0, {"CO2": n0, "CO": 0.0, "O2": 0.0})]
    comps += pick_compositions(records, [x for x in args.x if x > 0])
    tds = np.geomspace(args.td_min, args.td_max, args.points)
    table = []
    for x, state in comps:
        dens = {k: state.get(k, 0.0) for k in ("CO2", "CO", "O2")}
        row = [efficiency_at(swarm, float(td), dens) for td in tds]
        table.append(dict(x=x, rows=row))
        best = min(range(len(tds)), key=lambda i: row[i]["cost_eV"])
        print(f"X={100 * x:5.1f}%  最適 E/N={tds[best]:6.1f} Td  "
              f"コスト={row[best]['cost_eV']:6.2f} eV  "
              f"効率={100 * row[best]['efficiency']:5.1f}%", flush=True)
    (run / "efficiency_map.json").write_text(
        json.dumps(dict(td=tds.tolist(), table=table), ensure_ascii=False),
        encoding="utf-8")
    plot(run, tds, table)


def plot(run, tds, table):
    r"""曲線の図と、E/N × 分解率の2次元図を書く。

    Parameters
    ----------
    run : pathlib.Path
        出力フォルダ。
    tds : numpy.ndarray
        E/N の配列 [Td]。
    table : list of dict
        組成ごとの結果。
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for entry in table:
        eff = [100 * r["efficiency"] for r in entry["rows"]]
        axes[0].plot(tds, eff, label=f"X={100 * entry['x']:.0f}%")
    axes[0].set(xscale="log", xlabel="E/N [Td]", ylabel="効率 [%]",
                title="CO2 一次解離の効率（2.93 eV / コスト）")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    base = table[0]["rows"]
    bottom = np.zeros(len(tds))
    for g in GROUPS:
        v = np.array([r["partition"][g] for r in base])
        axes[1].fill_between(tds, bottom, bottom + v, label=g, alpha=0.85)
        bottom += v
    axes[1].set(xscale="log", xlabel="E/N [Td]", ylabel="入力電力の割合",
                title="純CO2の電力配分", ylim=(0, 1))
    axes[1].legend(fontsize=8, loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    fig.savefig(run / "efficiency_curves.png", dpi=150)
    plt.close(fig)

    xs = np.array([100 * t["x"] for t in table])
    z = np.array([[100 * r["efficiency"] for r in t["rows"]] for t in table])
    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(xs, tds, z.T, shading="gouraud",
                         norm=LogNorm(vmin=max(z.min(), 0.1), vmax=z.max()),
                         cmap="viridis")
    best = tds[np.argmax(z, axis=1)]
    ax.plot(xs, best, "w-o", lw=2, ms=5, label="最適 E/N")
    ax.set(yscale="log", xlabel="分解率 X [%]", ylabel="E/N [Td]",
           title="CO2 一次解離の効率 [%]（組成は飽和までのパルス列）")
    ax.legend(loc="upper right")
    fig.colorbar(mesh, label="効率 [%]")
    fig.tight_layout()
    fig.savefig(run / "efficiency_map.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
