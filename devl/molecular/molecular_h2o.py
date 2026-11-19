"""CLI zum Simulieren eines Clusters aus Wassermolekülen."""

import argparse
from math import isfinite
from typing import Sequence

from dirabo.args.check_types import non_negative_float, positive_float, positive_int
from dirabo.chem.args import add_chem_render_args, CHEM_RENDER_ARGS
from dirabo.chem.molecule_cluster import MoleculeCluster
from dirabo.chem.particles import CoulombForceModel
from dirabo.chem.render import show_molecules
from dirabo.chem.water import create_water_molecule


def create_water_cluster(
    lake_size: int = 2,
    *,
    placement_clearance: float = 0.0,
) -> MoleculeCluster:
    """Erzeugt einen Cluster, indem jedes Wassermolekül einzeln hinzugefügt wird."""
    if lake_size < 1:
        raise ValueError("Die Seegröße muss mindestens ein Wassermolekül betragen")
    if not isfinite(placement_clearance) or placement_clearance < 0:
        raise ValueError("Der Platzierungsabstand muss endlich und nicht negativ sein")

    cluster = MoleculeCluster(create_water_molecule("H2O-1"))
    for index in range(1, lake_size):
        cluster.add(
            create_water_molecule(f"H2O-{index + 1}"),
            clearance=placement_clearance,
        )
    return cluster


def _parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI-Argumente für die Wassermolekül-Simulation."""

    parser = argparse.ArgumentParser(description="Ein See aus interagierenden Wassermolekülen.")

    # da wir die Moleküle einzeln hinzufügen, ist die Seegröße die Anzahl der Moleküle
    parser.add_argument(
        "--lake-size",
        type=positive_int,
        default=2,
        metavar="N",
        help="Anzahl der H2O-Moleküle (Standard: 2)",
    )

    # weitere
    add_chem_render_args(parser)

    return parser.parse_args(args)


def main(args: Sequence[str] | None = None) -> None:
    options = _parse_args(args)
    show_molecules(
        create_water_cluster(
            options.lake_size,
            nucleus_repulsion_factor=CHEM_RENDER_ARGS.nucleus_repulsion_factor,
            nucleus_repulsion_radius=CHEM_RENDER_ARGS.nucleus_repulsion_radius,
        ),
        start_paused=CHEM_RENDER_ARGS.stop_at_start,
    )


if __name__ == "__main__":
    main()
