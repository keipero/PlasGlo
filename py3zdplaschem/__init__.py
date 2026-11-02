r"""モジュール化した0次元化学反応モデルの公開API。
"""

from .chemistry.reactions import Reaction, read_reactions
from .chemistry.mechanism import Mechanism
from .electrons.swarm_table import SwarmData, read_swarm
from .solvers.simulation import Result, simulate
from .diagnostics.balance import species_balance
from .visualization.plots import plot_result, plot_species_balance
from .io.results import save_result
from .io.cases import Case, load_case

__all__ = [
    "Reaction",
    "read_reactions",
    "Mechanism",
    "SwarmData",
    "read_swarm",
    "Result",
    "simulate",
    "species_balance",
    "plot_result",
    "plot_species_balance",
    "save_result",
    "Case",
    "load_case",
]
