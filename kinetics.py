r"""互換API・CLI。実装本体はplasglo配下にある。
"""

from plasglo import *
from plasglo.constants import KB, KB_EV
from plasglo.paths import WORKSPACE_DIR as ROOT, CHEMISTRY, BOLSIG
from plasglo.chemistry.reactions import parse_side, composition
from plasglo.chemistry.rates import analytic_rate
from plasglo.solvers.simulation import _integrate_coupled
from plasglo.cli import main

if __name__ == "__main__":
    main()
