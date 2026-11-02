r"""ZDPlasKin比較用の互換CLI・API。
"""

from plasglo.chemistry.sample_mechanism import (
    SampleMechanism,
    Expression,
    side,
    number,
)
from plasglo.electrons.sample_swarm import solve_sample_swarm
from plasglo.io.reference_data import read_csv, driver_conditions
from plasglo.validation.zdplaskin import (
    replay_archive,
    native_control,
    build_report,
    add_all_species_verification,
    write_markdown,
    main,
)

if __name__ == "__main__":
    main()
