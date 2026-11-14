"""CLI zum Visualisieren einfacher aromatischer Ringmuster."""

import argparse
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Sequence

from dirabo.chem.aromatic_functional_groups import create_aromatic_ring
from dirabo.chem.math import Vec3
from dirabo.chem.organic_functional_groups import (
    FunctionalGroupType,
    OrganicFunctionalGroup,
)
from dirabo.chem.particles import Molecule, MoleculeCluster
from dirabo.chem.render import show_molecules


@dataclass(frozen=True)
class RingTemplate:
    name: str
    elements: tuple[str, ...]
    lone_pair_donor_indices: tuple[int, ...] = ()


RING_TEMPLATES = {
    "benzol": RingTemplate("Benzol", ("C",) * 6),
    "pyridin": RingTemplate("Pyridin", ("N", "C", "C", "C", "C", "C")),
    "pyrrol": RingTemplate("Pyrrol", ("N", "C", "C", "C", "C"), (0,)),
    "furan": RingTemplate("Furan", ("O", "C", "C", "C", "C"), (0,)),
    "imidazol": RingTemplate("Imidazol", ("N", "C", "N", "C", "C"), (0,)),
    "pyrimidin": RingTemplate("Pyrimidin", ("N", "C", "N", "C", "C", "C")),
    "oxazol": RingTemplate("Oxazol", ("O", "C", "N", "C", "C"), (0,)),
}
MAX_RING_CARBONS = max(template.elements.count("C") for template in RING_TEMPLATES.values())


def create_ring_molecule(
    ring_name: str,
    instance_id: str | None = None,
    position: Vec3 = (0.0, 0.0, 0.0),
    *,
    substituents: Mapping[int, FunctionalGroupType] | None = None,
) -> Molecule:
    """Erzeugt eines der vordefinierten Ringmuster als Molekül."""
    try:
        template = RING_TEMPLATES[ring_name]
    except KeyError:
        choices = ", ".join(RING_TEMPLATES)
        raise ValueError(f"Unbekannter Ring {ring_name!r}; verfügbar: {choices}") from None

    molecule_id = instance_id or ring_name
    ring: OrganicFunctionalGroup = create_aromatic_ring(
        molecule_id,
        template.elements,
        position=position,
        lone_pair_donor_indices=template.lone_pair_donor_indices,
        substituents=substituents,
    )
    return ring.to_molecule(name=f"{template.name} {molecule_id}")


def create_ring_cluster(
    ring_name: str,
    instance_id: str | None = None,
    position: Vec3 = (0.0, 0.0, 0.0),
    *,
    substituents: Mapping[int, FunctionalGroupType] | None = None,
) -> MoleculeCluster:
    """Erzeugt einen simulierbaren Cluster für ein vordefiniertes Ringmuster."""
    return MoleculeCluster(
        create_ring_molecule(
            ring_name,
            instance_id,
            position,
            substituents=substituents,
        )
    )


def _parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aromatische Ringmuster als vereinfachte Molekülmodelle visualisieren."
    )
    parser.add_argument(
        "--ring",
        choices=tuple(RING_TEMPLATES),
        default="benzol",
        help="Ringmuster (Standard: benzol)",
    )
    for carbon_number in range(1, MAX_RING_CARBONS + 1):
        parser.add_argument(
            f"--C{carbon_number}",
            dest=f"substituent_C{carbon_number}",
            choices=tuple(group_type.value for group_type in FunctionalGroupType),
            metavar="GRUPPE",
            help=f"Substituent am {carbon_number}. Kohlenstoffatom",
        )

    options = parser.parse_args(args)
    carbon_indices = [index for index, element in enumerate(RING_TEMPLATES[options.ring].elements) if element == "C"]
    substituents: dict[int, FunctionalGroupType] = {}
    for carbon_number in range(1, MAX_RING_CARBONS + 1):
        group_name = getattr(options, f"substituent_C{carbon_number}")
        if group_name is None:
            continue
        if carbon_number > len(carbon_indices):
            parser.error(
                f"--C{carbon_number} ist für den Ring {options.ring!r} "
                f"nicht verfügbar (nur {len(carbon_indices)} C-Atome)"
            )
        substituents[carbon_indices[carbon_number - 1]] = FunctionalGroupType(group_name)
    options.substituents = substituents
    return options


def main(args: Sequence[str] | None = None) -> None:
    options = _parse_args(args)
    show_molecules(create_ring_cluster(options.ring, substituents=options.substituents))


if __name__ == "__main__":
    main()
