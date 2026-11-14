"""CLI zum Simulieren eines Clusters aus Wassermolekülen."""

import argparse
from math import isfinite
from typing import Sequence

from dirabo.chem.molecule_cluster import MoleculeCluster
from dirabo.chem.particles import CoulombForceModel
from dirabo.chem.render import show_molecules
from dirabo.chem.water import create_water_molecule

DEFAULT_NUCLEUS_REPULSION_FACTOR = 0.2
DEFAULT_NUCLEUS_REPULSION_RADIUS: float | None = None


def create_water_cluster(
    lake_size: int = 2,
    *,
    placement_clearance: float = 0.0,
    nucleus_repulsion_factor: float = DEFAULT_NUCLEUS_REPULSION_FACTOR,
    nucleus_repulsion_radius: float | None = DEFAULT_NUCLEUS_REPULSION_RADIUS,
) -> MoleculeCluster:
    """Erzeugt einen Cluster, indem jedes Wassermolekül einzeln hinzugefügt wird."""
    if lake_size < 1:
        raise ValueError("Die Seegröße muss mindestens ein Wassermolekül betragen")
    if not isfinite(placement_clearance) or placement_clearance < 0:
        raise ValueError("Der Platzierungsabstand muss endlich und nicht negativ sein")

    cluster = MoleculeCluster(
        create_water_molecule("H2O-1"),
        force_model=CoulombForceModel(
            inter_molecular_nucleus_repulsion_factor=nucleus_repulsion_factor,
            inter_molecular_nucleus_repulsion_radius=nucleus_repulsion_radius,
        ),
    )
    for index in range(1, lake_size):
        cluster.add(
            create_water_molecule(f"H2O-{index + 1}"),
            clearance=placement_clearance,
        )
    return cluster


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
        "--stop-at-start",
        action="store_true",
        help="Simulation zunächst pausiert starten",
    )
    parser.add_argument(
        "--nucleus-repulsion-factor",
        type=_non_negative_float,
        default=DEFAULT_NUCLEUS_REPULSION_FACTOR,
        metavar="FAKTOR",
        help=(
            "Kraftmultiplikator bei Kernabstand 0; nähert sich mit Abstand "
            f"dem Basiswert 1 (Standard: {DEFAULT_NUCLEUS_REPULSION_FACTOR})"
        ),
    )
    parser.add_argument(
        "--nucleus-repulsion-radius",
        type=_positive_float,
        default=DEFAULT_NUCLEUS_REPULSION_RADIUS,
        metavar="RADIUS",
        help=(
            "Optionaler Radius-Override; bei diesem Abstand ist der Exzess "
            "über Faktor 1 auf 10 % abgeklungen (Standard: Atomradien)"
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
        ),
        start_paused=options.stop_at_start,
    )


if __name__ == "__main__":
    main()
