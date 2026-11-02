r"""ZDPlasKin比較用の互換CLI・API。
"""

from py3zdplaschem.chemistry.sample_mechanism import (
    SampleMechanism,
    Expression,
    side,
    number,
)
from py3zdplaschem.electrons.sample_swarm import solve_sample_swarm
from py3zdplaschem.io.reference_data import read_csv, driver_conditions
from py3zdplaschem.validation.zdplaskin import (
    replay_archive,
    native_control,
    build_report,
    add_all_species_verification,
    write_markdown,
    main,
)

if __name__ == "__main__":
    main()
