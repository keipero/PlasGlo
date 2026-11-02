r"""互換用の公開名。BOLOS接続本体はplasglo/electrons配下にある。
"""

from plasglo.electrons.bolos_backend import (
    BOLOS_DIR,
    CROSS_SECTIONS,
    THREE_BODY_TARGET,
    CHANNELS,
    grid,
    parser,
    solver,
    normalize_eedf,
    BolosSwarm,
)
