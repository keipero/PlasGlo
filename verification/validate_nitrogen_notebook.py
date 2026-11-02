r"""新しいローカルPythonカーネルで全セルを実行し、出力を検証する。
"""

from pathlib import Path
import argparse
import json
import os
import sys

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager

cli = argparse.ArgumentParser(description=__doc__)
cli.add_argument(
    "--python",
    type=Path,
    default=Path(sys.executable),
    help="Python interpreter used by the notebook kernel",
)
args = cli.parse_args()
kernel_python = args.python.resolve()
if not kernel_python.is_file():
    raise FileNotFoundError(kernel_python)

MODEL = Path(__file__).resolve().parents[1]
NOTEBOOK = MODEL / "notebooks/run_nitrogen.ipynb"
runtime = MODEL / "outputs/verification/jupyter_runtime"
runtime.mkdir(parents=True, exist_ok=True)
os.environ["JUPYTER_RUNTIME_DIR"] = str(runtime)
os.environ["IPYTHONDIR"] = str(runtime / "ipython")
os.environ["MPLCONFIGDIR"] = str(runtime / "matplotlib")

notebook = nbformat.read(NOTEBOOK, as_version=4)
print("Starting a fresh Python kernel...", flush=True)
manager = KernelManager(kernel_name="python3")
manager.kernel_spec.argv[0] = str(kernel_python)
client = NotebookClient(
    notebook,
    timeout=180,
    kernel_name="python3",
    km=manager,
    resources={"metadata": {"path": str(NOTEBOOK.parent)}},
    allow_errors=False,
    record_timing=False,
    startup_timeout=30,
    on_cell_start=lambda cell, cell_index: print(
        f"Cell {cell_index + 1}/{len(notebook.cells)}", flush=True
    ),
)
client.execute()
nbformat.validate(notebook)
code_cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
assert all(cell.execution_count is not None for cell in code_cells)
assert not [
    output
    for cell in code_cells
    for output in cell.outputs
    if output.output_type == "error"
]
assert any(
    "image/png" in output.get("data", {})
    for cell in code_cells
    for output in cell.outputs
)
nbformat.write(notebook, NOTEBOOK)
output = MODEL / "outputs/notebooks/nitrogen"
for filename in [
    "densities.csv",
    "densities.npz",
    "coefficients.csv",
    "diagnostics.json",
    "densities.png",
    "pathways/N/reaction_rates.npz",
    "pathways/N/reaction_ranking.csv",
    "pathways/N/pathway_graph.json",
    "pathways/N/contributions.png",
]:
    assert (output / filename).is_file(), filename
report = json.loads((output / "diagnostics.json").read_text(encoding="utf-8"))
assert report["species_count"] == 21
assert report["reaction_count"] == 96
assert Path(report["python_executable"]).resolve() == kernel_python
graph = json.loads(
    (output / "pathways/N/pathway_graph.json").read_text(encoding="utf-8")
)
assert graph["species"] == "N"
assert graph["nodes"] and graph["edges"]
assert graph["weight_unit"] == "m^-3"
result = dict(
    notebook=str(NOTEBOOK),
    executed_code_cells=len(code_cells),
    kernel_python=report["python_executable"],
    inline_plot_verified=True,
    pathway_outputs_verified=True,
    species=21,
    reactions=96,
    max_atom_relative_error=report["max_atom_relative_error"],
    output_directory=str(output),
)
(MODEL / "outputs/verification/nitrogen_notebook_validation.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(result, ensure_ascii=False, indent=2))
