"""Visualize configured numbers of molecules from a JSON group library."""

import argparse
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Sequence

from dirabo.args.check_types import non_negative_int
import dirabo.chem.args as chem_args
from dirabo.chem.generic_molecular_groups import GenericMolecularGroups
from dirabo.chem.molecule_cluster import MoleculeCluster
from dirabo.chem.particles import Molecule
from dirabo.chem.render import show_molecules

DEFAULT_CONFIG = Path(__file__).with_suffix(".json")
RESERVED_OPTIONS = chem_args.RESERVED_OPTIONS.copy()
RESERVED_OPTIONS.update({"config", "list-molecules"})


def _config_path(args: Sequence[str] | None) -> Path:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    options, _ = parser.parse_known_args(args)
    return options.config


def _parse_args(
    args: Sequence[str] | None = None,
    molecule_names: tuple[str, ...] = (),
) -> argparse.Namespace:
    """Parse CLI-Argumente für die Molekül-Simulation."""
    parser = argparse.ArgumentParser(
        description=("Molekülanzahlen aus der JSON-Gruppenbibliothek angeben, " "z. B. --ethanol 2 --acetone 1.")
    )

    # config file für die Definition der Molekülgruppen
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"JSON-Gruppenbibliothek (Standard: {DEFAULT_CONFIG.name})",
    )

    # optionale Auflistung der verfügbaren Moleküle
    parser.add_argument(
        "--list-molecules",
        action="store_true",
        help="Verfügbare Moleküloptionen auflisten und beenden",
    )

    # weitere Parameter für die Chemie-Visualisierung hinzufügen
    chem_args.add_chem_render_args(parser)

    destinations: dict[str, str] = {}
    for index, molecule_name in enumerate(molecule_names):
        if molecule_name in RESERVED_OPTIONS or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", molecule_name) is None:
            parser.error(f"Molekülname {molecule_name!r} kann nicht als CLI-Option verwendet werden")
        destination = f"molecule_count_{index}"
        destinations[molecule_name] = destination
        parser.add_argument(
            f"--{molecule_name}",
            dest=destination,
            type=non_negative_int,
            default=0,
            metavar="ANZAHL",
            help=f"Anzahl der Instanzen von {molecule_name}",
        )

    options = parser.parse_args(args)
    options.molecule_counts = {
        molecule_name: getattr(options, destination) for molecule_name, destination in destinations.items()
    }
    if not options.list_molecules and not any(options.molecule_counts.values()):
        parser.error("Mindestens eine Molekülanzahl muss größer als 0 sein")
    return options


def create_cluster(
    library: GenericMolecularGroups,
    molecule_counts: Mapping[str, int],
) -> MoleculeCluster:
    unknown_names = molecule_counts.keys() - set(library.molecule_names)
    if unknown_names:
        raise ValueError(f"Unbekannte Moleküle: {', '.join(sorted(unknown_names))}")
    if any(not isinstance(count, int) or isinstance(count, bool) or count < 0 for count in molecule_counts.values()):
        raise ValueError("Molekülanzahlen müssen nichtnegative ganze Zahlen sein")
    if not any(molecule_counts.values()):
        raise ValueError("Mindestens eine Molekülanzahl muss größer als 0 sein")

    cluster: MoleculeCluster | None = None
    for molecule_name in library.molecule_names:
        for index in range(molecule_counts.get(molecule_name, 0)):
            molecule: Molecule = library.create_molecule(
                molecule_name,
                instance_id=f"{molecule_name}-{index + 1}",
            )
            if cluster is None:
                cluster = MoleculeCluster(molecule)
            else:
                cluster.add(molecule)

    if cluster is None:
        raise ValueError("Die angegebenen Molekülanzahlen erzeugen keinen Cluster")
    return cluster


def main(args: Sequence[str] | None = None) -> None:
    config_path = _config_path(args)
    try:
        library = GenericMolecularGroups.load(config_path)
    except (OSError, ValueError) as error:
        raise SystemExit(f"{Path(__file__).name}: error: {error}") from error

    options = _parse_args(args, library.molecule_names)
    if options.list_molecules:
        print("\n".join(library.molecule_names))
        return

    try:
        cluster = create_cluster(library, options.molecule_counts)
    except (ValueError, RuntimeError) as error:
        raise SystemExit(f"{Path(__file__).name}: error: {error}") from error

    show_molecules(cluster, start_paused=options.stop_at_start)


if __name__ == "__main__":
    main()
