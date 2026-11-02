r"""プロジェクトのデータ、ローカル依存クローン、出力先の定義。
"""

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
WORKSPACE_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUTS_DIR = PROJECT_DIR / "outputs"
MECHANISMS_DIR = DATA_DIR / "mechanisms"
CROSS_SECTIONS_DIR = DATA_DIR / "cross_sections"
CHEMISTRY = MECHANISMS_DIR / "air_chemistry_Komuro_reactions.txt"
BOLSIG = CROSS_SECTIONS_DIR / "N2-O2_80-20_2000Td.dat"
CROSS_SECTIONS = CROSS_SECTIONS_DIR / "Phelps_N2-O2_20261026.txt"
BOLOS_DIR = PROJECT_DIR / "Py3bolos"
SAMPLES_DIR = PROJECT_DIR / "samples/gases"
NOTEBOOK_OUTPUT_DIR = OUTPUTS_DIR / "notebooks/nitrogen"
SAMPLE_OUTPUT_DIR = OUTPUTS_DIR / "samples"
VERIFICATION_OUTPUT_DIR = OUTPUTS_DIR / "verification"
