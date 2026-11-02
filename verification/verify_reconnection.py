r"""接続テストを再実行し、サンプル履歴をZDPlasKinと比較する。
"""

from importlib import import_module

from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys

import numpy as np

MODEL = Path(__file__).resolve().parents[1]
WORKSPACE = MODEL.parent
REFERENCE = WORKSPACE.parent / "ZDplasKin_analysis"
sys.path.insert(0, str(MODEL))
VERIFICATION_OUTPUT_DIR = import_module(
    "py3zdplaschem.paths"
).VERIFICATION_OUTPUT_DIR
SAMPLE_OUTPUT_DIR = import_module("py3zdplaschem.paths").SAMPLE_OUTPUT_DIR
native_control = import_module(
    "py3zdplaschem.validation.zdplaskin"
).native_control
read_csv = import_module("py3zdplaschem.io.reference_data").read_csv

OUT = VERIFICATION_OUTPUT_DIR


def metrics(names, predicted, reference):
    r"""粒子種ごとのピーク密度で規格化した時間履歴の差を計算する。

    Parameters
    ----------
    names : sequence of str
        密度配列の行と同じ順序の粒子種名。
    predicted : numpy.ndarray
        形状が (粒子種数, 時刻数) の予測密度 [m^-3]。
    reference : numpy.ndarray
        形状が (粒子種数, 時刻数) の参照密度 [m^-3]。

    Returns
    -------
    dict
        規格化最大差、差が最大となる粒子種、種ごとの誤差。
    """
    scale = np.maximum(np.max(np.abs(reference), axis=1), 1.0)
    error = np.max(np.abs(predicted - reference), axis=1) / scale
    return dict(
        max_absolute_error_over_species_peak=float(error.max()),
        worst_species=names[int(error.argmax())],
        per_species=dict(zip(names, error.tolist())),
    )


def argon_control():
    r"""元の2反応ケースを高精度CSV出力で再ビルドする。

    Returns
    -------
    tuple
        粒子名、ネイティブ時刻 [s]、SI密度履歴。

    Raises
    ------
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。
    """
    target = OUT / "argon_native"
    target.mkdir(parents=True, exist_ok=True)
    for name in ["kinet.inp", "zdplaskin_m.F90"]:
        shutil.copy2(REFERENCE / "build/argon" / name, target / name)
    shutil.copy2(
        REFERENCE / "examples/official1/bolsigdb.dat", target / "bolsigdb.dat"
    )
    shutil.copy2(
        REFERENCE / "runs/argon/20261021_000542_588761/bolsig_x86_64_g.dll",
        target / "bolsig_x86_64_g.dll",
    )
    shutil.copy2(
        REFERENCE
        / "ZDPlasKin_2.0a_Windows/ZDPlasKin_2.0a_Windows/dvode_f90_m.F90",
        target / "dvode_f90_m.F90",
    )
    # Same original conditions; only CSV output precision is changed.
    shutil.copy2(
        MODEL / "samples/gases/argon/native_csv_driver.F90",
        target / "user_code.F90",
    )
    compiler = Path(
        "C:/rtools44/x86_64-w64-mingw32.static.posix/bin/gfortran.exe"
    )
    env = os.environ.copy()
    env["PATH"] = str(compiler.parent) + os.pathsep + env.get("PATH", "")
    temporary = target / "tmp"
    temporary.mkdir(exist_ok=True)
    env["TEMP"] = env["TMP"] = str(temporary.resolve())
    commands = [
        (
            "compile.log",
            [
                str(compiler),
                "-std=legacy",
                "-O2",
                "-ffree-line-length-none",
                "dvode_f90_m.F90",
                "zdplaskin_m.F90",
                "user_code.F90",
                "bolsig_x86_64_g.dll",
                "-o",
                "run_control.exe",
            ],
        ),
        ("stdout.log", [str((target / "run_control.exe").resolve())]),
    ]
    for log, command in commands:
        completed = subprocess.run(
            command,
            cwd=target,
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        (target / log).write_text(
            completed.stdout + completed.stderr, encoding="utf-8"
        )
        if completed.returncode:
            raise RuntimeError(
                completed.stdout[-1000:] + completed.stderr[-2000:]
            )
    return read_csv(target / "densities.csv")


def main():
    r"""モデルのテストとサンプル履歴をネイティブ計算で検証する。

    Raises
    ------
    RuntimeError
        ソルバーが収束しないか、計算した物理量が不正な場合。
    ValueError
        入力が要求する形状、単位、モデルの制約を満たさない場合。

    Examples
    --------
    プロジェクトのルートから実行します::

        python verification/verify_reconnection.py
    """
    OUT.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [
            sys.executable,
            "-X",
            "utf8",
            "-m",
            "unittest",
            "discover",
            "-s",
            str(MODEL),
            "-p",
            "test_*.py",
            "-v",
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        timeout=120,
    )
    (OUT / "unittest.log").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(completed.stdout + completed.stderr)
    summary = {"unittests_passed": True, "cases": {}}
    for name in ["nitrogen", "argon"]:
        print(f"Checking {name} against a fresh native control...", flush=True)
        folder = SAMPLE_OUTPUT_DIR / name
        with np.load(folder / "densities.npz") as data:
            names, time, predicted = (
                data["species"].tolist(),
                data["time_s"],
                data["densities_m3"],
            )
        reference_dir = REFERENCE / (
            "runs/nitrogen_1000Td_10ns/20261021_234859_884"
            if name == "nitrogen"
            else "runs/argon/20261021_000542_588761"
        )
        input_files = list((MODEL / "samples/gases" / name).glob("*"))
        input_files += list(reference_dir.glob("*.F90")) + list(
            reference_dir.glob("*.inp")
        )
        input_files += [
            reference_dir / "densities.csv",
            reference_dir / "bolsigdb.dat",
        ]
        input_files = [path for path in input_files if path.is_file()]
        before = {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in input_files
        }
        if name == "nitrogen":
            control_names, control_time, control = native_control(
                reference_dir, OUT, disable_ee=True
            )
        else:
            control_names, control_time, control = argon_control()
        if [n.upper() for n in control_names] != names:
            raise ValueError("Species order mismatch")
        np.testing.assert_allclose(time, control_time, rtol=1e-12, atol=1e-25)
        result = metrics(names, predicted, control)
        after = {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in input_files
        }
        if before != after:
            raise RuntimeError("An original or staged sample input changed")
        result["inputs_unchanged"] = True
        result["native_electron_electron_collisions"] = (
            False if name == "nitrogen" else "original configuration"
        )
        # Also isolate chemistry integration from Boltzmann coefficient
        # differences.
        if name == "argon":
            from scipy.integrate import solve_ivp
            from samples.gases.run_samples import run_case

            _, _, _, mechanism, _ = run_case("argon")
            coefficients = np.loadtxt(
                OUT / "argon_native/coefficients.csv",
                delimiter=",",
                skiprows=1,
            )
            mechanism.k[:] = coefficients[0, 1:3] * [1e-6, 1e-12]
            mechanism.scaled = mechanism.k * mechanism.n0 ** (
                mechanism.orders - 1
            )
            sol = solve_ivp(
                mechanism.rhs,
                (0, time[-1]),
                control[:, 0] / mechanism.n0,
                t_eval=time,
                method="BDF",
                jac=mechanism.jac,
                rtol=1e-9,
                atol=1e-4 / mechanism.n0,
            )
            if not sol.success:
                raise RuntimeError(sol.message)
            result["same_native_coefficients_chemistry_error"] = metrics(
                names, sol.y * mechanism.n0, control
            )
            result["bolos_ionization_coefficient_m3_s"] = float(
                np.loadtxt(
                    folder / "coefficients.csv",
                    delimiter=",",
                    skiprows=1,
                    usecols=2,
                )[0]
            )
            result["native_ionization_coefficient_m3_s"] = float(
                coefficients[0, 1] * 1e-6
            )
        summary["cases"][name] = result
        print(json.dumps({name: result}, indent=2), flush=True)
    status = subprocess.run(
        ["git", "-C", str(MODEL / "Py3bolos"), "status", "--porcelain"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    summary["clone_tracked_files_unchanged"] = not bool(status.strip())
    if status.strip():
        raise RuntimeError(f"Unexpected clone changes: {status}")
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
