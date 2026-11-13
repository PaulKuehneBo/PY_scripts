"""Erzeuger und Komposition organischer funktioneller Gruppen."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from math import acos, cos, radians, sin, sqrt

from dirabo.chem.math import Vec3, add, cross, length, scale
from dirabo.chem.particles import (
    Atom,
    Bond,
    Electron,
    ElectronPair,
    Molecule,
    create_bond,
)
from dirabo.chem.psystems import carbon, hydrogen, nitrogen, oxygen

SINGLE_BOND_ELECTRON_OFFSET = 0.12
PI_BOND_OFFSET = 0.22


class FunctionalGroupType(str, Enum):
    HYDROXYL = "hydroxyl"
    CARBONYL = "carbonyl"
    CARBOXYL = "carboxyl"
    ETHER = "ether"
    ESTER = "ester"
    AMINO = "amino"
    AMIDE = "amide"
    NITRILE = "nitrile"


@dataclass(frozen=True)
class OrganicFunctionalGroup:
    """Ein internes Molekülfragment mit offenen Bindungsplätzen."""

    name: str
    atoms: tuple[Atom, ...]
    bonds: tuple[Bond, ...]
    lone_pairs: tuple[ElectronPair, ...]
    attachment_atom_ids: tuple[str, ...]
    free_electrons: tuple[Electron, ...] = ()

    def __post_init__(self) -> None:
        known_atom_ids = {atom.atom_id for atom in self.atoms}
        if any(atom_id not in known_atom_ids for atom_id in self.attachment_atom_ids):
            raise ValueError("Alle Bindungsplätze müssen auf Atome des Fragments verweisen")
        self.to_molecule()

    def to_molecule(self, name: str | None = None) -> Molecule:
        """Erzeugt ein Molecule für Simulation und Rendering des Fragments."""
        return Molecule(
            name=name or self.name,
            atoms=self.atoms,
            bonds=self.bonds,
            lone_pairs=self.lone_pairs,
            free_electrons=self.free_electrons,
        )


@dataclass(frozen=True)
class FunctionalGroupConnection:
    """Verknüpft zwei offene Bindungsplätze aus einer Fragmentliste."""

    first_group: int
    first_attachment: int
    second_group: int
    second_attachment: int
    bond_id: str


@dataclass(frozen=True)
class _Frame:
    origin: Vec3
    x_axis: Vec3
    y_axis: Vec3
    z_axis: Vec3

    def point(self, x: float, y: float = 0.0, z: float = 0.0) -> Vec3:
        return add(
            self.origin,
            add(
                scale(self.x_axis, x),
                add(scale(self.y_axis, y), scale(self.z_axis, z)),
            ),
        )

    def vector(self, x: float, y: float = 0.0, z: float = 0.0) -> Vec3:
        return add(
            scale(self.x_axis, x),
            add(scale(self.y_axis, y), scale(self.z_axis, z)),
        )


def _frame(position: Vec3, direction: Vec3) -> _Frame:
    direction_length = length(direction)
    if direction_length == 0:
        raise ValueError("Die Ausrichtung einer funktionellen Gruppe darf nicht null sein")

    x_axis = scale(direction, 1.0 / direction_length)
    reference = (0.0, 0.0, 1.0) if abs(x_axis[2]) < 0.9 else (0.0, 1.0, 0.0)
    y_axis = cross(reference, x_axis)
    y_axis = scale(y_axis, 1.0 / length(y_axis))
    z_axis = cross(x_axis, y_axis)
    return _Frame(position, x_axis, y_axis, z_axis)


def _check_instance_id(instance_id: str) -> None:
    if not instance_id:
        raise ValueError("Eine funktionelle Gruppe benötigt eine eindeutige ID")


def _carbon_chain_positions(
    count: int,
    position: Vec3,
    bond_length: float,
) -> tuple[Vec3, ...]:
    positions = [position]
    heading = 0.0
    turn_angle = acos(1.0 / 3.0)
    for link_index in range(1, count):
        heading += turn_angle if link_index % 2 else -turn_angle
        step = (cos(heading) * bond_length, sin(heading) * bond_length, 0.0)
        positions.append(add(positions[-1], step))
    return tuple(positions)


def _carbon_valence_directions(
    neighbor_directions: tuple[Vec3, ...],
) -> tuple[Vec3, ...]:
    if not neighbor_directions:
        tetrahedral_scale = 1.0 / sqrt(3.0)
        return tuple(
            scale(direction, tetrahedral_scale)
            for direction in (
                (1.0, 1.0, 1.0),
                (1.0, -1.0, -1.0),
                (-1.0, 1.0, -1.0),
                (-1.0, -1.0, 1.0),
            )
        )

    frame = _frame((0.0, 0.0, 0.0), neighbor_directions[0])
    if len(neighbor_directions) == 1:
        tetrahedral_axis = scale(neighbor_directions[0], -1.0 / 3.0)
        radial_scale = sqrt(8.0 / 9.0)
        return tuple(
            add(
                tetrahedral_axis,
                scale(
                    add(
                        scale(frame.y_axis, cos(radians(angle))),
                        scale(frame.z_axis, sin(radians(angle))),
                    ),
                    radial_scale,
                ),
            )
            for angle in (0.0, 120.0, 240.0)
        )

    if len(neighbor_directions) != 2:
        raise ValueError("Eine lineare Kohlenstoffkette darf höchstens zwei Nachbarn je C haben")

    first, second = neighbor_directions
    normal = cross(first, second)
    normal_length = length(normal)
    if normal_length == 0:
        raise ValueError("Die Richtungen der C-C-Bindungen dürfen nicht kollinear sein")
    normal = scale(normal, 1.0 / normal_length)
    inward_bisector = scale(add(first, second), -0.5)
    out_of_plane_component = scale(normal, sqrt(2.0 / 3.0))
    return (
        add(inward_bisector, out_of_plane_component),
        add(inward_bisector, scale(out_of_plane_component, -1.0)),
    )


def create_functional_group(
    group_type: FunctionalGroupType,
    instance_id: str,
    position: Vec3,
    direction: Vec3,
) -> OrganicFunctionalGroup:
    """Erzeugt einen funktionellen Gruppenbaustein anhand seines Enum-Werts."""
    if group_type is FunctionalGroupType.HYDROXYL:
        return create_hydroxyl_group(instance_id, position, direction)
    if group_type is FunctionalGroupType.CARBONYL:
        return create_carbonyl_group(instance_id, position, direction)
    if group_type is FunctionalGroupType.CARBOXYL:
        return create_carboxyl_group(instance_id, position, direction)
    if group_type is FunctionalGroupType.ETHER:
        return create_ether_group(instance_id, position)
    if group_type is FunctionalGroupType.ESTER:
        return create_ester_group(instance_id, position, direction)
    if group_type is FunctionalGroupType.AMINO:
        return create_amino_group(instance_id, position, direction)
    if group_type is FunctionalGroupType.AMIDE:
        return create_amide_group(instance_id, position, direction)
    if group_type is FunctionalGroupType.NITRILE:
        return create_nitrile_group(instance_id, position, direction)
    raise ValueError(f"Unbekannter funktioneller Gruppentyp: {group_type!r}")


def functional_group_attachment_distance(group_type: FunctionalGroupType) -> float:
    """Gibt den Modellabstand vom Ankeratom zum Bindungspartner zurück."""
    if group_type in (
        FunctionalGroupType.HYDROXYL,
        FunctionalGroupType.ETHER,
    ):
        return 1.43
    if group_type is FunctionalGroupType.AMINO:
        return 1.47
    if group_type is FunctionalGroupType.NITRILE:
        return 1.46
    return 1.54


def create_carbon_chain(
    instance_id: str,
    hydrogen_counts: Sequence[int],
    *,
    substituents: Mapping[int, FunctionalGroupType] | None = None,
    position: Vec3 = (0.0, 0.0, 0.0),
    bond_length: float = 1.54,
) -> OrganicFunctionalGroup:
    """Erzeugt eine gesättigte C-Kette; Substituenten-Indizes sind nullbasiert.

    ``hydrogen_counts`` gibt die Anzahl der H-Atome je C vor. Ein
    ``substituents``-Eintrag bindet genau eine funktionelle Gruppe an das
    entsprechende C. Nicht verwendete C-Valenzen bleiben als offene
    Bindungsplätze im zurückgegebenen Fragment erhalten. Beispiel:
    ``create_carbon_chain("chain", [3, 2], substituents={1: FunctionalGroupType.HYDROXYL})``.
    """
    _check_instance_id(instance_id)
    counts = tuple(hydrogen_counts)
    if not counts:
        raise ValueError("Eine Kohlenstoffkette benötigt mindestens ein C-Atom")
    if any(not isinstance(count, int) or isinstance(count, bool) or count < 0 for count in counts):
        raise ValueError("Die Zahl der H-Atome je C muss eine nichtnegative ganze Zahl sein")
    if bond_length <= 0:
        raise ValueError("Die C-C-Bindungslänge muss positiv sein")

    substituent_types = dict(substituents or {})
    if any(
        not isinstance(index, int) or isinstance(index, bool) or not 0 <= index < len(counts)
        for index in substituent_types
    ):
        raise ValueError("Substituenten-Indizes müssen gültige nullbasierte C-Indizes sein")
    if any(not isinstance(group_type, FunctionalGroupType) for group_type in substituent_types.values()):
        raise TypeError("Substituenten müssen Werte aus FunctionalGroupType sein")

    positions = _carbon_chain_positions(len(counts), position, bond_length)
    carbon_atoms = tuple(
        carbon(f"{instance_id}:C{index + 1}", atom_position, lone_pair_count=0).atom
        for index, atom_position in enumerate(positions)
    )
    chain_bonds = tuple(
        create_bond(
            f"{instance_id}:C{index + 1}-C{index + 2}",
            carbon_atoms[index],
            carbon_atoms[index + 1],
            electron_z_offsets=(-SINGLE_BOND_ELECTRON_OFFSET, SINGLE_BOND_ELECTRON_OFFSET),
        )
        for index in range(len(carbon_atoms) - 1)
    )

    hydrogen_atoms: list[Atom] = []
    hydrogen_bonds: list[Bond] = []
    substituent_groups: list[OrganicFunctionalGroup] = []
    connections: list[FunctionalGroupConnection] = []
    open_carbon_sites: list[str] = []

    for index, atom in enumerate(carbon_atoms):
        neighbor_directions: list[Vec3] = []
        if index > 0:
            neighbor_directions.append(
                scale(
                    add(atom.position, scale(carbon_atoms[index - 1].position, -1.0)),
                    1.0 / bond_length,
                )
            )
        if index + 1 < len(carbon_atoms):
            neighbor_directions.append(
                scale(
                    add(atom.position, scale(carbon_atoms[index + 1].position, -1.0)),
                    1.0 / bond_length,
                )
            )

        free_directions = _carbon_valence_directions(tuple(neighbor_directions))
        substituent_count = int(index in substituent_types)
        occupied_valences = len(neighbor_directions) + counts[index] + substituent_count
        if occupied_valences > 4:
            raise ValueError(
                f"C{index + 1} hat {occupied_valences} Bindungen; Kohlenstoff darf höchstens vier haben"
            )

        carbon_attachment_index = len(open_carbon_sites)
        open_carbon_sites.extend(
            [atom.atom_id] * (4 - len(neighbor_directions) - counts[index])
        )

        for hydrogen_index, direction in enumerate(free_directions[: counts[index]], start=1):
            hydrogen_atom = hydrogen(
                f"{instance_id}:C{index + 1}:H{hydrogen_index}",
                add(atom.position, scale(direction, 1.09)),
            ).atom
            hydrogen_atoms.append(hydrogen_atom)
            hydrogen_bonds.append(
                create_bond(
                    f"{instance_id}:C{index + 1}-H{hydrogen_index}",
                    atom,
                    hydrogen_atom,
                    electron_z_offsets=(-SINGLE_BOND_ELECTRON_OFFSET, SINGLE_BOND_ELECTRON_OFFSET),
                )
            )

        if substituent_count:
            direction = free_directions[counts[index]]
            group_type = substituent_types[index]
            group_position = add(
                atom.position,
                scale(direction, functional_group_attachment_distance(group_type)),
            )
            group_id = f"{instance_id}:C{index + 1}:{group_type.value}"
            group = create_functional_group(group_type, group_id, group_position, direction)
            if not group.attachment_atom_ids:
                raise ValueError(f"Die Gruppe {group_type.value} hat keinen offenen Bindungsplatz")
            substituent_groups.append(group)
            connections.append(
                FunctionalGroupConnection(
                    first_group=0,
                    first_attachment=carbon_attachment_index,
                    second_group=len(substituent_groups),
                    second_attachment=0,
                    bond_id=f"{instance_id}:C{index + 1}-{group_type.value}",
                )
            )

    chain_fragment = OrganicFunctionalGroup(
        name=f"Kohlenstoffkette {instance_id}",
        atoms=(*carbon_atoms, *hydrogen_atoms),
        bonds=(*chain_bonds, *hydrogen_bonds),
        lone_pairs=(),
        attachment_atom_ids=tuple(open_carbon_sites),
    )
    if not substituent_groups:
        return chain_fragment

    return combine_functional_groups(
        name=f"Kohlenstoffkette {instance_id}",
        groups=(chain_fragment, *substituent_groups),
        connections=tuple(connections),
    )


def _bond(
    instance_id: str,
    bond_name: str,
    first_atom: Atom,
    second_atom: Atom,
    frame: _Frame,
    *,
    pi_offset: Vec3 | None = None,
) -> Bond:
    bond_id = f"{instance_id}:{bond_name}"
    if pi_offset is None:
        offset = frame.vector(0.0, 0.0, SINGLE_BOND_ELECTRON_OFFSET)
        electron_direction_offsets = (offset, scale(offset, -1.0))
    else:
        offset = frame.vector(*pi_offset)
        electron_direction_offsets = (offset, offset)

    return create_bond(
        bond_id,
        first_atom,
        second_atom,
        electron_z_offsets=(0.0, 0.0),
        electron_direction_offsets=electron_direction_offsets,
    )


def _group(
    name: str,
    atoms: tuple[Atom, ...],
    bonds: tuple[Bond, ...],
    lone_pairs: tuple[ElectronPair, ...],
    attachment_atoms: tuple[Atom, ...],
) -> OrganicFunctionalGroup:
    return OrganicFunctionalGroup(
        name=name,
        atoms=atoms,
        bonds=bonds,
        lone_pairs=lone_pairs,
        attachment_atom_ids=tuple(atom.atom_id for atom in attachment_atoms),
    )


def create_hydroxyl_group(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
    direction: Vec3 = (1.0, 0.0, 0.0),
) -> OrganicFunctionalGroup:
    """Erzeugt –OH; der offene Bindungsplatz liegt am Sauerstoff."""
    _check_instance_id(instance_id)
    frame = _frame(position, direction)
    oxygen_particles = oxygen(f"{instance_id}:O", position, lone_pair_count=2)
    hydrogen_particles = hydrogen(f"{instance_id}:H", frame.point(0.96))
    bond = _bond(instance_id, "O-H", oxygen_particles.atom, hydrogen_particles.atom, frame)
    return _group(
        f"Hydroxyl {instance_id}",
        (oxygen_particles.atom, hydrogen_particles.atom),
        (bond,),
        oxygen_particles.lone_pairs,
        (oxygen_particles.atom,),
    )


def create_carbonyl_group(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
    direction: Vec3 = (1.0, 0.0, 0.0),
) -> OrganicFunctionalGroup:
    """Erzeugt >C=O mit zwei offenen Bindungsplätzen am Carbonyl-Kohlenstoff."""
    _check_instance_id(instance_id)
    frame = _frame(position, direction)
    carbon_particles = carbon(f"{instance_id}:C", position, lone_pair_count=0)
    oxygen_particles = oxygen(f"{instance_id}:O", frame.point(1.23), lone_pair_count=2)
    bonds = (
        _bond(instance_id, "C=O:sigma", carbon_particles.atom, oxygen_particles.atom, frame),
        _bond(
            instance_id,
            "C=O:pi",
            carbon_particles.atom,
            oxygen_particles.atom,
            frame,
            pi_offset=(0.0, 0.0, PI_BOND_OFFSET),
        ),
    )
    return _group(
        f"Carbonyl {instance_id}",
        (carbon_particles.atom, oxygen_particles.atom),
        bonds,
        oxygen_particles.lone_pairs,
        (carbon_particles.atom, carbon_particles.atom),
    )


def create_carboxyl_group(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
    direction: Vec3 = (1.0, 0.0, 0.0),
) -> OrganicFunctionalGroup:
    """Erzeugt –COOH mit einem offenen Bindungsplatz am Carboxyl-Kohlenstoff."""
    _check_instance_id(instance_id)
    frame = _frame(position, direction)
    carbon_particles = carbon(f"{instance_id}:C", position, lone_pair_count=0)
    carbonyl_oxygen = oxygen(f"{instance_id}:O1", frame.point(1.23), lone_pair_count=2)
    hydroxyl_oxygen = oxygen(f"{instance_id}:O2", frame.point(-0.67, 1.12), lone_pair_count=2)
    hydroxyl_hydrogen = hydrogen(f"{instance_id}:H", frame.point(-1.48, 0.63))
    atoms = (
        carbon_particles.atom,
        carbonyl_oxygen.atom,
        hydroxyl_oxygen.atom,
        hydroxyl_hydrogen.atom,
    )
    bonds = (
        _bond(instance_id, "C=O:sigma", carbon_particles.atom, carbonyl_oxygen.atom, frame),
        _bond(
            instance_id,
            "C=O:pi",
            carbon_particles.atom,
            carbonyl_oxygen.atom,
            frame,
            pi_offset=(0.0, 0.0, PI_BOND_OFFSET),
        ),
        _bond(instance_id, "C-OH", carbon_particles.atom, hydroxyl_oxygen.atom, frame),
        _bond(instance_id, "O-H", hydroxyl_oxygen.atom, hydroxyl_hydrogen.atom, frame),
    )
    return _group(
        f"Carboxyl {instance_id}",
        atoms,
        bonds,
        (*carbonyl_oxygen.lone_pairs, *hydroxyl_oxygen.lone_pairs),
        (carbon_particles.atom,),
    )


def create_ether_group(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
) -> OrganicFunctionalGroup:
    """Erzeugt –O– mit zwei offenen Bindungsplätzen am Sauerstoff."""
    _check_instance_id(instance_id)
    oxygen_particles = oxygen(f"{instance_id}:O", position, lone_pair_count=2)
    return _group(
        f"Ether {instance_id}",
        (oxygen_particles.atom,),
        (),
        oxygen_particles.lone_pairs,
        (oxygen_particles.atom, oxygen_particles.atom),
    )


def create_ester_group(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
    direction: Vec3 = (1.0, 0.0, 0.0),
) -> OrganicFunctionalGroup:
    """Erzeugt –COO– mit offenen Plätzen am Carbonyl-Kohlenstoff und O."""
    _check_instance_id(instance_id)
    frame = _frame(position, direction)
    carbon_particles = carbon(f"{instance_id}:C", position, lone_pair_count=0)
    carbonyl_oxygen = oxygen(f"{instance_id}:O1", frame.point(1.23), lone_pair_count=2)
    ester_oxygen = oxygen(f"{instance_id}:O2", frame.point(-0.67, 1.12), lone_pair_count=2)
    bonds = (
        _bond(instance_id, "C=O:sigma", carbon_particles.atom, carbonyl_oxygen.atom, frame),
        _bond(
            instance_id,
            "C=O:pi",
            carbon_particles.atom,
            carbonyl_oxygen.atom,
            frame,
            pi_offset=(0.0, 0.0, PI_BOND_OFFSET),
        ),
        _bond(instance_id, "C-O", carbon_particles.atom, ester_oxygen.atom, frame),
    )
    return _group(
        f"Ester {instance_id}",
        (carbon_particles.atom, carbonyl_oxygen.atom, ester_oxygen.atom),
        bonds,
        (*carbonyl_oxygen.lone_pairs, *ester_oxygen.lone_pairs),
        (carbon_particles.atom, ester_oxygen.atom),
    )


def create_amino_group(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
    direction: Vec3 = (1.0, 0.0, 0.0),
) -> OrganicFunctionalGroup:
    """Erzeugt –NH2 mit einem offenen Bindungsplatz am Stickstoff."""
    _check_instance_id(instance_id)
    frame = _frame(position, direction)
    nitrogen_particles = nitrogen(f"{instance_id}:N", position, lone_pair_count=1)
    hydrogen1 = hydrogen(f"{instance_id}:H1", frame.point(-0.34, 0.91, 0.43))
    hydrogen2 = hydrogen(f"{instance_id}:H2", frame.point(-0.34, 0.91, -0.43))
    bonds = (
        _bond(instance_id, "N-H1", nitrogen_particles.atom, hydrogen1.atom, frame),
        _bond(instance_id, "N-H2", nitrogen_particles.atom, hydrogen2.atom, frame),
    )
    return _group(
        f"Amino {instance_id}",
        (nitrogen_particles.atom, hydrogen1.atom, hydrogen2.atom),
        bonds,
        nitrogen_particles.lone_pairs,
        (nitrogen_particles.atom,),
    )


def create_amide_group(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
    direction: Vec3 = (1.0, 0.0, 0.0),
) -> OrganicFunctionalGroup:
    """Erzeugt –CONH2 mit einem offenen Bindungsplatz am Carbonyl-Kohlenstoff."""
    _check_instance_id(instance_id)
    frame = _frame(position, direction)
    carbon_particles = carbon(f"{instance_id}:C", position, lone_pair_count=0)
    oxygen_particles = oxygen(f"{instance_id}:O", frame.point(1.23), lone_pair_count=2)
    nitrogen_particles = nitrogen(f"{instance_id}:N", frame.point(-0.67, 1.12), lone_pair_count=1)
    hydrogen1 = hydrogen(f"{instance_id}:H1", frame.point(-1.04, 2.00, 0.40))
    hydrogen2 = hydrogen(f"{instance_id}:H2", frame.point(-1.04, 2.00, -0.40))
    bonds = (
        _bond(instance_id, "C=O:sigma", carbon_particles.atom, oxygen_particles.atom, frame),
        _bond(
            instance_id,
            "C=O:pi",
            carbon_particles.atom,
            oxygen_particles.atom,
            frame,
            pi_offset=(0.0, 0.0, PI_BOND_OFFSET),
        ),
        _bond(instance_id, "C-N", carbon_particles.atom, nitrogen_particles.atom, frame),
        _bond(instance_id, "N-H1", nitrogen_particles.atom, hydrogen1.atom, frame),
        _bond(instance_id, "N-H2", nitrogen_particles.atom, hydrogen2.atom, frame),
    )
    return _group(
        f"Amide {instance_id}",
        (carbon_particles.atom, oxygen_particles.atom, nitrogen_particles.atom, hydrogen1.atom, hydrogen2.atom),
        bonds,
        (*oxygen_particles.lone_pairs, *nitrogen_particles.lone_pairs),
        (carbon_particles.atom,),
    )


def create_nitrile_group(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
    direction: Vec3 = (1.0, 0.0, 0.0),
) -> OrganicFunctionalGroup:
    """Erzeugt –C≡N mit einem offenen Bindungsplatz am Kohlenstoff."""
    _check_instance_id(instance_id)
    frame = _frame(position, direction)
    carbon_particles = carbon(f"{instance_id}:C", position, lone_pair_count=0)
    nitrogen_particles = nitrogen(f"{instance_id}:N", frame.point(1.16), lone_pair_count=1)
    bonds = (
        _bond(instance_id, "C#N:sigma", carbon_particles.atom, nitrogen_particles.atom, frame),
        _bond(
            instance_id,
            "C#N:pi-1",
            carbon_particles.atom,
            nitrogen_particles.atom,
            frame,
            pi_offset=(0.0, 0.0, PI_BOND_OFFSET),
        ),
        _bond(
            instance_id,
            "C#N:pi-2",
            carbon_particles.atom,
            nitrogen_particles.atom,
            frame,
            pi_offset=(0.0, PI_BOND_OFFSET, 0.0),
        ),
    )
    return _group(
        f"Nitril {instance_id}",
        (carbon_particles.atom, nitrogen_particles.atom),
        bonds,
        nitrogen_particles.lone_pairs,
        (carbon_particles.atom,),
    )


def combine_functional_groups(
    name: str,
    groups: tuple[OrganicFunctionalGroup, ...],
    connections: tuple[FunctionalGroupConnection, ...] = (),
) -> OrganicFunctionalGroup:
    """Kombiniert Fragmente und verbindet angegebene offene Bindungsplätze."""
    if not name:
        raise ValueError("Das kombinierte Fragment benötigt einen Namen")
    if not groups:
        raise ValueError("Mindestens eine funktionelle Gruppe muss angegeben werden")

    atoms = tuple(atom for group in groups for atom in group.atoms)
    bonds = tuple(bond for group in groups for bond in group.bonds)
    lone_pairs = tuple(pair for group in groups for pair in group.lone_pairs)
    free_electrons = tuple(electron for group in groups for electron in group.free_electrons)

    atom_by_id = {atom.atom_id: atom for atom in atoms}
    bond_ids = {bond.bond_id for bond in bonds}
    electron_ids = {
        electron.electron_id
        for group in groups
        for bond in group.bonds
        for electron in bond.electrons
    }
    attachment_slots = [list(group.attachment_atom_ids) for group in groups]
    used_slots: set[tuple[int, int]] = set()

    for connection in connections:
        first_slot = (connection.first_group, connection.first_attachment)
        second_slot = (connection.second_group, connection.second_attachment)
        for group_index, attachment_index in (first_slot, second_slot):
            if not 0 <= group_index < len(groups):
                raise ValueError(f"Unbekannter Gruppenindex: {group_index}")
            if not 0 <= attachment_index < len(attachment_slots[group_index]):
                raise ValueError(f"Unbekannter Bindungsplatz: {attachment_index}")
        if first_slot == second_slot:
            raise ValueError("Ein Bindungsplatz kann nicht mit sich selbst verbunden werden")
        if first_slot in used_slots or second_slot in used_slots:
            raise ValueError("Jeder Bindungsplatz darf höchstens einmal verbunden werden")
        if not connection.bond_id or connection.bond_id in bond_ids:
            raise ValueError(f"Die Bindungs-ID {connection.bond_id!r} ist leer oder bereits vergeben")

        first_atom_id = attachment_slots[first_slot[0]][first_slot[1]]
        second_atom_id = attachment_slots[second_slot[0]][second_slot[1]]
        connection_bond = create_bond(
            connection.bond_id,
            atom_by_id[first_atom_id],
            atom_by_id[second_atom_id],
        )
        new_electron_ids = {electron.electron_id for electron in connection_bond.electrons}
        if electron_ids.intersection(new_electron_ids):
            raise ValueError("Die erzeugte Bindung verwendet bereits vergebene Elektronen-IDs")
        bonds += (connection_bond,)
        bond_ids.add(connection.bond_id)
        electron_ids.update(new_electron_ids)
        used_slots.update((first_slot, second_slot))

    remaining_attachments = tuple(
        atom_id
        for group_index, group_slots in enumerate(attachment_slots)
        for attachment_index, atom_id in enumerate(group_slots)
        if (group_index, attachment_index) not in used_slots
    )
    return OrganicFunctionalGroup(
        name=name,
        atoms=atoms,
        bonds=bonds,
        lone_pairs=lone_pairs,
        attachment_atom_ids=remaining_attachments,
        free_electrons=free_electrons,
    )
