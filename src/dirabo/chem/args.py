from dataclasses import dataclass

from dirabo.args.check_types import non_negative_float, positive_float, positive_int


@dataclass
class ChemRenderArgs:
    stop_at_start: bool
    nucleus_repulsion_factor: float
    nucleus_repulsion_radius: float


CHEM_RENDER_ARGS = ChemRenderArgs(
    stop_at_start=False,
    nucleus_repulsion_factor=5.0,
    nucleus_repulsion_radius=5.0,
)


def add_chem_render_args(parser):
    parser.add_argument(
        "--stop-at-start",
        action="store_true",
        help="Simulation zunächst pausiert starten",
    )
    parser.add_argument(
        "--nucleus-repulsion-factor",
        type=non_negative_float,
        default=CHEM_RENDER_ARGS.nucleus_repulsion_factor,
        metavar="FAKTOR",
        help=(
            "Kraftmultiplikator bei Kernabstand 0; nähert sich mit Abstand "
            f"dem Basiswert 1 (Standard: {CHEM_RENDER_ARGS.nucleus_repulsion_factor})"
        ),
    )
    parser.add_argument(
        "--nucleus-repulsion-radius",
        type=positive_float,
        default=CHEM_RENDER_ARGS.nucleus_repulsion_radius,
        metavar="RADIUS",
        help=(
            "Optionaler Radius-Override; bei diesem Abstand ist der Exzess "
            "über Faktor 1 auf 10 % abgeklungen (Standard: Atomradien)"
        ),
    )


RESERVED_OPTIONS = {"nucleus-repulsion-factor", "nucleus-repulsion-radius", "stop-at-start", "help"}
