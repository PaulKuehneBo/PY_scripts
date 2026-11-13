"""CLI zum Simulieren eines Clusters aus Wassermolekülen."""

import argparse
from math import ceil, isfinite, sqrt
from typing import Sequence

from dirabo.chem.particles import CoulombForceModel, MoleculeCluster
from dirabo.chem.render import show_molecules
from dirabo.chem.water import create_water_molecule

WATER_MOLECULE_SPACING = 3.2
DEFAULT_NUCLEUS_REPULSION_FACTOR = 0.2
DEFAULT_NUCLEUS_REPULSION_RADIUS = 5.0


def create_water_cluster(
    lake_size: int = 2,
    spacing: float = WATER_MOLECULE_SPACING,
    nucleus_repulsion_factor: float = DEFAULT_NUCLEUS_REPULSION_FACTOR,
    nucleus_repulsion_radius: float | None = DEFAULT_NUCLEUS_REPULSION_RADIUS,
) -> MoleculeCluster:
    """Erzeugt einen zentrierten, rasterförmigen Cluster aus Wassermolekülen."""
    if lake_size < 1:
        raise ValueError("Die Seegröße muss mindestens ein Wassermolekül betragen")
    if spacing <= 0:
        raise ValueError("Der Abstand zwischen Wassermolekülen muss positiv sein")

    column_count = ceil(sqrt(lake_size))
    grid_positions = [
        (
            (index % column_count) * spacing,
            (index // column_count) * spacing,
            0.0,
        )
        for index in range(lake_size)
    ]
    center_x = sum(position[0] for position in grid_positions) / lake_size
    center_y = sum(position[1] for position in grid_positions) / lake_size

    molecules = tuple(
        create_water_molecule(
            f"H2O-{index + 1}",
            (position[0] - center_x, position[1] - center_y, position[2]),
        )
        for index, position in enumerate(grid_positions)
    )
    return MoleculeCluster(
        molecules,
        force_model=CoulombForceModel(
            inter_molecular_nucleus_repulsion_factor=nucleus_repulsion_factor,
            inter_molecular_nucleus_repulsion_radius=nucleus_repulsion_radius,
        ),
    )


def _positive_int(value: str) -> int:
    try:
        size = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("muss eine ganze Zahl sein") from error
    if size < 1:
        raise argparse.ArgumentTypeError("muss mindestens 1 sein")
    return size


def _non_negative_float(value: str) -> float:
    try:
        factor = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("muss eine nicht negative Zahl sein") from error
    if not isfinite(factor) or factor < 0:
        raise argparse.ArgumentTypeError("muss endlich und mindestens 0 sein")
    return factor


def _positive_float(value: str) -> float:
    try:
        radius = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("muss eine positive Zahl sein") from error
    if not isfinite(radius) or radius <= 0:
        raise argparse.ArgumentTypeError("muss endlich und größer als 0 sein")
    return radius


def _parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Ein See aus interagierenden Wassermolekülen."
    )
    parser.add_argument(
        "--lake-size",
        type=_positive_int,
        default=2,
        metavar="N",
        help="Anzahl der H2O-Moleküle (Standard: 2)",
    )
    parser.add_argument(
        "--nucleus-repulsion-factor",
        type=_non_negative_float,
        default=DEFAULT_NUCLEUS_REPULSION_FACTOR,
        metavar="FAKTOR",
        help=(
            "Skalierung der Kern-Kern-Abstoßung zwischen verschiedenen "
            f"Molekülen, 0=aus (Standard: {DEFAULT_NUCLEUS_REPULSION_FACTOR})"
        ),
    )
    parser.add_argument(
        "--nucleus-repulsion-radius",
        type=_positive_float,
        default=DEFAULT_NUCLEUS_REPULSION_RADIUS,
        metavar="RADIUS",
        help=(
            "Maximaler Abstand zwischen Atomkernen für die zusätzliche "
            f"Abstoßung (Standard: {DEFAULT_NUCLEUS_REPULSION_RADIUS})"
        ),
    )
    return parser.parse_args(args)


def main(args: Sequence[str] | None = None) -> None:
    options = _parse_args(args)
    show_molecules(
        create_water_cluster(
            options.lake_size,
            nucleus_repulsion_factor=options.nucleus_repulsion_factor,
            nucleus_repulsion_radius=options.nucleus_repulsion_radius,
        )
    )


if __name__ == "__main__":
    main()
