"""CLI zum Simulieren eines Clusters aus Wassermolekülen."""

import argparse
from math import ceil, sqrt
from typing import Sequence

from dirabo.chem.particles import MoleculeCluster
from dirabo.chem.render import show_molecules
from dirabo.chem.water import create_water_molecule

WATER_MOLECULE_SPACING = 3.2


def create_water_cluster(
    lake_size: int = 2,
    spacing: float = WATER_MOLECULE_SPACING,
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
    return MoleculeCluster(molecules)


def _positive_int(value: str) -> int:
    try:
        size = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("muss eine ganze Zahl sein") from error
    if size < 1:
        raise argparse.ArgumentTypeError("muss mindestens 1 sein")
    return size


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
    return parser.parse_args(args)


def main(args: Sequence[str] | None = None) -> None:
    options = _parse_args(args)
    show_molecules(create_water_cluster(options.lake_size))


if __name__ == "__main__":
    main()
