"""Erzeuger für aromatische Ringfragmente mit impliziten Wasserstoffatomen."""

from collections.abc import Mapping, Sequence
from math import cos, pi, sin

from dirabo.chem.math import Vec3, add, scale
from dirabo.chem.organic_functional_groups import (
    FunctionalGroupConnection,
    FunctionalGroupType,
    OrganicFunctionalGroup,
    combine_functional_groups,
    create_functional_group,
    functional_group_attachment_distance,
)
from dirabo.chem.particles import Atom, Bond, create_bond
from dirabo.chem.psystems import ElementParticles, carbon, hydrogen, nitrogen, oxygen

RING_BOND_ELECTRON_OFFSET = 0.12
PI_BOND_INWARD_OFFSET = 0.22
CARBON_HYDROGEN_BOND_LENGTH = 1.09
NITROGEN_HYDROGEN_BOND_LENGTH = 1.01

_ELEMENT_VALENCES = {"C": 4, "N": 3, "O": 2}


def _default_pi_bond_indices(
    ring_size: int,
    lone_pair_donor_indices: set[int],
) -> tuple[int, ...]:
    candidates: list[tuple[int, ...]] = []
    for parity in (0, 1):
        used_vertices: set[int] = set()
        selected: list[int] = []
        for edge_index in range(parity, ring_size, 2):
            next_index = (edge_index + 1) % ring_size
            if (
                edge_index in lone_pair_donor_indices
                or next_index in lone_pair_donor_indices
                or edge_index in used_vertices
                or next_index in used_vertices
            ):
                continue
            selected.append(edge_index)
            used_vertices.update((edge_index, next_index))
        candidates.append(tuple(selected))
    return max(candidates, key=len)


def _make_ring_atom(element: str, atom_id: str, position: Vec3) -> ElementParticles:
    if element == "C":
        return carbon(atom_id, position, lone_pair_count=0)
    if element == "N":
        return nitrogen(atom_id, position, lone_pair_count=1)
    if element == "O":
        return oxygen(atom_id, position, lone_pair_count=2)
    raise ValueError(f"Element {element!r} wird für aromatische Ringatome nicht unterstützt")


def _hydrogen_direction(
    radial_direction: Vec3,
    normal: Vec3,
    hydrogen_index: int,
    hydrogen_count: int,
) -> Vec3:
    if hydrogen_count == 1:
        return radial_direction
    side_direction = (
        -radial_direction[1],
        radial_direction[0],
        0.0,
    )
    if hydrogen_count == 2:
        side_scale = 0.5
        normal_scale = 0.5
        side_sign = 1.0 if hydrogen_index == 0 else -1.0
        return add(
            scale(radial_direction, 0.5),
            add(
                scale(side_direction, side_scale * side_sign),
                scale(normal, normal_scale),
            ),
        )
    return radial_direction


def _ring_bond(
    bond_id: str,
    first_atom: Atom,
    second_atom: Atom,
    normal_offset: Vec3,
) -> Bond:
    return create_bond(
        bond_id,
        first_atom,
        second_atom,
        electron_z_offsets=(0.0, 0.0),
        electron_direction_offsets=(normal_offset, normal_offset),
    )


def create_aromatic_ring(
    instance_id: str,
    ring_elements: Sequence[str],
    *,
    substituents: Mapping[int, FunctionalGroupType] | None = None,
    position: Vec3 = (0.0, 0.0, 0.0),
    bond_length: float = 1.39,
    pi_bond_indices: Sequence[int] | None = None,
    lone_pair_donor_indices: Sequence[int] = (),
) -> OrganicFunctionalGroup:
    """Erzeugt einen planaren aromatischen Ring und füllt freie Valenzen mit H.

    ``ring_elements`` gibt die Ringatome in Reihenfolge an. Bindungsindex ``i``
    bezeichnet die Kante von Atom ``i`` zu Atom ``(i + 1) % n``. Ohne
    ``pi_bond_indices`` wird eine nicht überlappende Kekulé-Belegung erzeugt.
    ``lone_pair_donor_indices`` nimmt die Ringatome auf, deren freies Paar zum
    aromatischen System beiträgt, etwa das N in Pyrrol. Substituenten-Indizes
    sind nullbasiert und ersetzen an dieser Position ein implizites H.
    """
    if not instance_id:
        raise ValueError("Ein aromatischer Ring benötigt eine eindeutige ID")
    elements = tuple(ring_elements)
    ring_size = len(elements)
    if ring_size < 3:
        raise ValueError("Ein Ring benötigt mindestens drei Ringatome")
    if any(element not in _ELEMENT_VALENCES for element in elements):
        raise ValueError("Aromatische Ringe unterstützen derzeit nur C, N und O")
    if bond_length <= 0:
        raise ValueError("Die Ring-Bindungslänge muss positiv sein")

    donor_indices = set(lone_pair_donor_indices)
    if any(
        not isinstance(index, int)
        or isinstance(index, bool)
        or not 0 <= index < ring_size
        for index in donor_indices
    ):
        raise ValueError("Lone-Pair-Donor-Indizes müssen gültige nullbasierte Ringpositionen sein")
    if any(elements[index] not in ("N", "O") for index in donor_indices):
        raise ValueError("Nur Stickstoff oder Sauerstoff können als Lone-Pair-Donor dienen")

    substituent_types = dict(substituents or {})
    if any(
        not isinstance(index, int)
        or isinstance(index, bool)
        or not 0 <= index < ring_size
        for index in substituent_types
    ):
        raise ValueError("Substituenten-Indizes müssen gültige nullbasierte Ringpositionen sein")
    if any(not isinstance(group_type, FunctionalGroupType) for group_type in substituent_types.values()):
        raise TypeError("Substituenten müssen Werte aus FunctionalGroupType sein")

    if pi_bond_indices is None:
        selected_pi_bonds = _default_pi_bond_indices(ring_size, donor_indices)
    else:
        selected_pi_bonds = tuple(pi_bond_indices)
    if any(
        not isinstance(index, int)
        or isinstance(index, bool)
        or not 0 <= index < ring_size
        for index in selected_pi_bonds
    ):
        raise ValueError("Pi-Bindungs-Indizes müssen gültige nullbasierte Kanten-Indizes sein")
    if len(set(selected_pi_bonds)) != len(selected_pi_bonds):
        raise ValueError("Pi-Bindungs-Indizes dürfen nicht doppelt vorkommen")

    pi_bond_degrees = [0] * ring_size
    for edge_index in selected_pi_bonds:
        next_index = (edge_index + 1) % ring_size
        pi_bond_degrees[edge_index] += 1
        pi_bond_degrees[next_index] += 1
        if edge_index in donor_indices or next_index in donor_indices:
            raise ValueError("Eine Pi-Bindung darf keinen Lone-Pair-Donor enthalten")
    if any(degree > 1 for degree in pi_bond_degrees):
        raise ValueError("Kekulé-Pi-Bindungen dürfen sich an Ringatomen nicht überlappen")
    if any(
        degree == 0 and index not in donor_indices
        for index, degree in enumerate(pi_bond_degrees)
    ):
        raise ValueError("Jedes Ringatom muss eine Pi-Bindung oder ein beitragendes Lone-Pair besitzen")
    pi_electron_count = 2 * len(selected_pi_bonds) + 2 * len(donor_indices)
    if pi_electron_count < 2 or (pi_electron_count - 2) % 4:
        raise ValueError("Die Pi-Elektronenzahl des Rings muss der 4n+2-Regel entsprechen")

    radius = bond_length / (2 * sin(pi / ring_size))
    radial_directions = tuple(
        (cos(2 * pi * index / ring_size), sin(2 * pi * index / ring_size), 0.0)
        for index in range(ring_size)
    )
    ring_particles = tuple(
        _make_ring_atom(
            element,
            f"{instance_id}:{element}{index + 1}",
            add(position, scale(radial_directions[index], radius)),
        )
        for index, element in enumerate(elements)
    )
    ring_atoms = tuple(particles.atom for particles in ring_particles)

    sigma_normal_offset = (0.0, 0.0, RING_BOND_ELECTRON_OFFSET)
    bonds: list[Bond] = []
    for edge_index, first_atom in enumerate(ring_atoms):
        next_index = (edge_index + 1) % ring_size
        second_atom = ring_atoms[next_index]
        bonds.append(
            _ring_bond(
                f"{instance_id}:{first_atom.atom_id.rsplit(':', 1)[-1]}-"
                f"{second_atom.atom_id.rsplit(':', 1)[-1]}:sigma",
                first_atom,
                second_atom,
                sigma_normal_offset,
            )
        )
        if edge_index in selected_pi_bonds:
            midpoint = scale(add(first_atom.position, second_atom.position), 0.5)
            inward_offset = scale(
                add(position, scale(midpoint, -1.0)),
                PI_BOND_INWARD_OFFSET / bond_length,
            )
            bonds.append(
                create_bond(
                    f"{instance_id}:{first_atom.atom_id.rsplit(':', 1)[-1]}-"
                    f"{second_atom.atom_id.rsplit(':', 1)[-1]}:pi",
                    first_atom,
                    second_atom,
                    electron_z_offsets=(0.0, 0.0),
                    electron_direction_offsets=(inward_offset, inward_offset),
                )
            )

    hydrogen_atoms: list[Atom] = []
    for index, atom in enumerate(ring_atoms):
        substituent_count = int(index in substituent_types)
        implicit_hydrogen_count = (
            _ELEMENT_VALENCES[elements[index]]
            - 2
            - pi_bond_degrees[index]
            - substituent_count
        )
        if implicit_hydrogen_count < 0:
            raise ValueError(
                f"Ringatom {index} ({elements[index]}) überschreitet seine Valenz"
            )

        hydrogen_bond_length = (
            NITROGEN_HYDROGEN_BOND_LENGTH if elements[index] == "N" else CARBON_HYDROGEN_BOND_LENGTH
        )
        for hydrogen_index in range(implicit_hydrogen_count):
            hydrogen_direction = _hydrogen_direction(
                radial_directions[index],
                (0.0, 0.0, 1.0),
                hydrogen_index,
                implicit_hydrogen_count,
            )
            hydrogen_atom = hydrogen(
                f"{instance_id}:{atom.atom_id.rsplit(':', 1)[-1]}:H{hydrogen_index + 1}",
                add(atom.position, scale(hydrogen_direction, hydrogen_bond_length)),
            ).atom
            hydrogen_atoms.append(hydrogen_atom)
            bonds.append(
                _ring_bond(
                    f"{instance_id}:{atom.atom_id.rsplit(':', 1)[-1]}-"
                    f"H{hydrogen_index + 1}",
                    atom,
                    hydrogen_atom,
                    sigma_normal_offset,
                )
            )

    ring_lone_pairs = tuple(pair for particles in ring_particles for pair in particles.lone_pairs)
    ring_attachment_atoms = tuple(ring_atoms[index] for index in sorted(substituent_types))
    ring_fragment = OrganicFunctionalGroup(
        name=f"Aromatischer Ring {instance_id}",
        atoms=(*ring_atoms, *hydrogen_atoms),
        bonds=tuple(bonds),
        lone_pairs=ring_lone_pairs,
        attachment_atom_ids=tuple(atom.atom_id for atom in ring_attachment_atoms),
    )
    if not substituent_types:
        return ring_fragment

    substituent_groups: list[OrganicFunctionalGroup] = []
    connections: list[FunctionalGroupConnection] = []
    for attachment_index, ring_index in enumerate(sorted(substituent_types)):
        group_type = substituent_types[ring_index]
        direction = radial_directions[ring_index]
        group_position = add(
            ring_atoms[ring_index].position,
            scale(direction, functional_group_attachment_distance(group_type)),
        )
        substituent = create_functional_group(
            group_type,
            f"{instance_id}:ring-{ring_index + 1}:{group_type.value}",
            group_position,
            direction,
        )
        if not substituent.attachment_atom_ids:
            raise ValueError(f"Die Gruppe {group_type.value} hat keinen offenen Bindungsplatz")
        substituent_groups.append(substituent)
        connections.append(
            FunctionalGroupConnection(
                first_group=0,
                first_attachment=attachment_index,
                second_group=attachment_index + 1,
                second_attachment=0,
                bond_id=f"{instance_id}:ring-{ring_index + 1}-{group_type.value}",
            )
        )

    return combine_functional_groups(
        name=f"Aromatischer Ring {instance_id}",
        groups=(ring_fragment, *substituent_groups),
        connections=tuple(connections),
    )
