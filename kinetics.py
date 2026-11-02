r"""互換API・CLI。実装本体はpy3zdplaschem配下にある。
"""

from py3zdplaschem import *
from py3zdplaschem.constants import KB, KB_EV
from py3zdplaschem.paths import WORKSPACE_DIR as ROOT, CHEMISTRY, BOLSIG
from py3zdplaschem.chemistry.reactions import parse_side, composition
from py3zdplaschem.chemistry.rates import analytic_rate
from py3zdplaschem.solvers.simulation import _integrate_coupled
from py3zdplaschem.cli import main

if __name__ == "__main__":
    main()
