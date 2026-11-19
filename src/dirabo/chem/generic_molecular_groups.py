"""Load and assemble molecule templates declared as JSON group graphs."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from math import cos, isfinite, radians, sin
from pathlib import Path
from typing import Any

from dirabo.chem.math import Vec3
from dirabo.chem.particles import Atom, Bond, Electron, ElectronPair, Molecule, create_bond
from dirabo.chem.psystems import create_element


@dataclass(frozen=True)
class _Component:
    component_id: str
    group_name: str
    position: Vec3
    rotation_degrees: Vec3


def _mapping(value: Any, context: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{context} muss ein JSON-Objekt sein")
    return value


def _check_keys(
    value: Mapping[str, Any],
    required: set[str],
    optional: set[str],
    context: str,
) -> None:
    missing = required - value.keys()
    extra = value.keys() - required - optional
    if missing:
        raise ValueError(f"{context} fehlt erforderliche Felder: {', '.join(sorted(missing))}")
    if extra:
        raise ValueError(f"{context} enthält unbekannte Felder: {', '.join(sorted(extra))}")


def _identifier(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} muss eine nichtleere Zeichenkette sein")
    return value


def _vector(value: Any, context: str) -> Vec3:
    if (
        not isinstance(value, list)
        or len(value) != 3
        or any(
            isinstance(component, bool)
            or not isinstance(component, (int, float))
            or not isfinite(component)
            for component in value
        )
    ):
        raise ValueError(f"{context} muss ein Zahlenvektor mit drei Komponenten sein")
    return float(value[0]), float(value[1]), float(value[2])


def _unique_ids(values: Sequence[str], context: str) -> None:
    if len(values) != len(set(values)):
        raise ValueError(f"{context} müssen eindeutig sein")


def _rotate_vector(vector: Vec3, rotation_degrees: Vec3) -> Vec3:
    """Apply X, then Y, then Z Euler rotations to a local vector."""
    x_angle, y_angle, z_angle = (radians(angle) for angle in rotation_degrees)
    x_cos, x_sin = cos(x_angle), sin(x_angle)
    rotated = (
        vector[0],
        x_cos * vector[1] - x_sin * vector[2],
        x_sin * vector[1] + x_cos * vector[2],
    )
    y_cos, y_sin = cos(y_angle), sin(y_angle)
    rotated = (
        y_cos * rotated[0] + y_sin * rotated[2],
        rotated[1],
        -y_sin * rotated[0] + y_cos * rotated[2],
    )
    z_cos, z_sin = cos(z_angle), sin(z_angle)
    return (
        z_cos * rotated[0] - z_sin * rotated[1],
        z_sin * rotated[0] + z_cos * rotated[1],
        rotated[2],
    )


class GenericMolecularGroups:
    """Validated JSON definitions for reusable groups and assembled molecules."""

    def __init__(self, definition: Mapping[str, Any]) -> None:
        self._definition = dict(definition)
        _check_keys(
            self._definition,
            {"schema_version", "groups", "molecules"},
            {"$schema"},
            "Konfiguration",
        )
        if (
            isinstance(self._definition.get("schema_version"), bool)
            or self._definition.get("schema_version") != 1
        ):
            raise ValueError("Nicht unterstützte oder fehlende schema_version; erwartet wird 1")
        self._groups = _mapping(self._definition.get("groups"), "groups")
        self._molecules = _mapping(self._definition.get("molecules"), "molecules")
        if not self._groups:
            raise ValueError("groups muss mindestens eine Gruppe enthalten")
        if not self._molecules:
            raise ValueError("molecules muss mindestens ein Molekül enthalten")
        self._validate_groups()
        self._validate_molecules()

    @classmethod
    def load(cls, path: str | Path) -> GenericMolecularGroups:
        config_path = Path(path)
        try:
            definition = json.loads(config_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as error:
            raise ValueError(
                f"Ungültiges JSON in {config_path} (Zeile {error.lineno}, Spalte {error.colno}): {error.msg}"
            ) from error
        except OSError as error:
            raise OSError(f"Konfigurationsdatei kann nicht gelesen werden: {config_path}") from error
        return cls(_mapping(definition, "Konfiguration"))

    @property
    def molecule_names(self) -> tuple[str, ...]:
        return tuple(self._molecules)

    def _validate_groups(self) -> None:
        for group_name, raw_group in self._groups.items():
            group = _mapping(raw_group, f"Gruppe {group_name!r}")
            _check_keys(group, {"atoms"}, {"bonds", "ports"}, f"Gruppe {group_name!r}")
            atoms = group.get("atoms")
            bonds = group.get("bonds", [])
            ports = group.get("ports", {})
            if not isinstance(atoms, list) or not atoms:
                raise ValueError(f"Gruppe {group_name!r} benötigt mindestens ein Atom")
            if not isinstance(bonds, list):
                raise ValueError(f"bonds der Gruppe {group_name!r} muss eine Liste sein")
            ports = _mapping(ports, f"ports der Gruppe {group_name!r}")

            atom_ids: list[str] = []
            for index, raw_atom in enumerate(atoms):
                atom = _mapping(raw_atom, f"Atom {index} in Gruppe {group_name!r}")
                _check_keys(
                    atom,
                    {"id", "element", "position"},
                    {"lone_pairs", "free_electrons"},
                    f"Atom {index} in Gruppe {group_name!r}",
                )
                atom_id = _identifier(atom.get("id"), f"id von Atom {index} in Gruppe {group_name!r}")
                atom_ids.append(atom_id)
                _identifier(atom.get("element"), f"element von Atom {atom_id!r}")
                _vector(atom.get("position"), f"position von Atom {atom_id!r}")
                if "lone_pairs" in atom and (
                    not isinstance(atom["lone_pairs"], int)
                    or isinstance(atom["lone_pairs"], bool)
                    or atom["lone_pairs"] < 0
                ):
                    raise ValueError(f"lone_pairs von Atom {atom_id!r} muss eine nichtnegative Ganzzahl sein")
                free_electrons = atom.get("free_electrons", [])
                if not isinstance(free_electrons, list):
                    raise ValueError(f"free_electrons von Atom {atom_id!r} muss eine Liste sein")
                for electron_index, direction in enumerate(free_electrons):
                    _vector(direction, f"free_electrons[{electron_index}] von Atom {atom_id!r}")
            _unique_ids(atom_ids, f"Atom-IDs in Gruppe {group_name!r}")
            known_atoms = set(atom_ids)

            bond_ids: list[str] = []
            for index, raw_bond in enumerate(bonds):
                bond = _mapping(raw_bond, f"Bindung {index} in Gruppe {group_name!r}")
                _check_keys(
                    bond,
                    {"id", "from", "to"},
                    {"electron_direction_offsets"},
                    f"Bindung {index} in Gruppe {group_name!r}",
                )
                bond_id = _identifier(bond.get("id"), f"id von Bindung {index} in Gruppe {group_name!r}")
                bond_ids.append(bond_id)
                first = _identifier(bond.get("from"), f"from von Bindung {bond_id!r}")
                second = _identifier(bond.get("to"), f"to von Bindung {bond_id!r}")
                if first not in known_atoms or second not in known_atoms or first == second:
                    raise ValueError(f"Bindung {bond_id!r} verweist auf ungültige Atome")
                offsets = bond.get("electron_direction_offsets", [[0.0, 0.0, 0.0]] * 2)
                if (
                    not isinstance(offsets, list)
                    or len(offsets) != 2
                ):
                    raise ValueError(
                        f"electron_direction_offsets von Bindung {bond_id!r} muss zwei Vektoren enthalten"
                    )
                for electron_index, offset in enumerate(offsets):
                    _vector(offset, f"electron_direction_offsets[{electron_index}] von Bindung {bond_id!r}")
            _unique_ids(bond_ids, f"Bindungs-IDs in Gruppe {group_name!r}")

            for port_name, raw_atom_id in ports.items():
                _identifier(port_name, f"Portname in Gruppe {group_name!r}")
                atom_id = _identifier(raw_atom_id, f"Port {port_name!r} in Gruppe {group_name!r}")
                if atom_id not in known_atoms:
                    raise ValueError(f"Port {port_name!r} verweist auf unbekanntes Atom {atom_id!r}")

    def _validate_molecules(self) -> None:
        for molecule_name, raw_molecule in self._molecules.items():
            molecule = _mapping(raw_molecule, f"Molekül {molecule_name!r}")
            _check_keys(molecule, {"components"}, {"name", "connections"}, f"Molekül {molecule_name!r}")
            components = molecule.get("components")
            connections = molecule.get("connections", [])
            if not isinstance(components, list) or not components:
                raise ValueError(f"Molekül {molecule_name!r} benötigt mindestens eine Komponente")
            if not isinstance(connections, list):
                raise ValueError(f"connections von Molekül {molecule_name!r} muss eine Liste sein")
            _identifier(molecule.get("name", molecule_name), f"name von Molekül {molecule_name!r}")

            component_ids: list[str] = []
            for index, raw_component in enumerate(components):
                component = _mapping(raw_component, f"Komponente {index} in Molekül {molecule_name!r}")
                _check_keys(
                    component,
                    {"id", "group"},
                    {"position", "rotation_degrees"},
                    f"Komponente {index} in Molekül {molecule_name!r}",
                )
                component_id = _identifier(
                    component.get("id"),
                    f"id von Komponente {index} in Molekül {molecule_name!r}",
                )
                component_ids.append(component_id)
                group_name = _identifier(component.get("group"), f"group von Komponente {component_id!r}")
                if group_name not in self._groups:
                    raise ValueError(f"Komponente {component_id!r} nutzt unbekannte Gruppe {group_name!r}")
                _vector(component.get("position", [0.0, 0.0, 0.0]), f"position von Komponente {component_id!r}")
                _vector(
                    component.get("rotation_degrees", [0.0, 0.0, 0.0]),
                    f"rotation_degrees von Komponente {component_id!r}",
                )
            _unique_ids(component_ids, f"Komponenten-IDs in Molekül {molecule_name!r}")
            known_components = set(component_ids)
            connection_ids: list[str] = []
            used_ports: set[tuple[str, str]] = set()
            for index, raw_connection in enumerate(connections):
                connection = _mapping(raw_connection, f"Verbindung {index} in Molekül {molecule_name!r}")
                _check_keys(
                    connection,
                    {"id", "from", "to"},
                    set(),
                    f"Verbindung {index} in Molekül {molecule_name!r}",
                )
                connection_id = _identifier(
                    connection.get("id"),
                    f"id von Verbindung {index} in Molekül {molecule_name!r}",
                )
                connection_ids.append(connection_id)
                endpoints = []
                for endpoint_name in ("from", "to"):
                    endpoint = _mapping(
                        connection.get(endpoint_name),
                        f"{endpoint_name} von Verbindung {connection_id!r}",
                    )
                    _check_keys(
                        endpoint,
                        {"component", "port"},
                        set(),
                        f"{endpoint_name} von Verbindung {connection_id!r}",
                    )
                    component_id = _identifier(
                        endpoint.get("component"),
                        f"component von {endpoint_name} in Verbindung {connection_id!r}",
                    )
                    port_name = _identifier(
                        endpoint.get("port"),
                        f"port von {endpoint_name} in Verbindung {connection_id!r}",
                    )
                    if component_id not in known_components:
                        raise ValueError(
                            f"Verbindung {connection_id!r} verweist auf unbekannte Komponente {component_id!r}"
                        )
                    component_definition = next(
                        raw_component
                        for raw_component in components
                        if raw_component["id"] == component_id
                    )
                    group = self._groups[component_definition["group"]]
                    if port_name not in group.get("ports", {}):
                        raise ValueError(
                            f"Gruppe {component_definition['group']!r} hat keinen Port {port_name!r}"
                        )
                    port_key = component_id, port_name
                    if port_key in used_ports:
                        raise ValueError(f"Port {port_name!r} der Komponente {component_id!r} wird mehrfach genutzt")
                    used_ports.add(port_key)
                    endpoints.append(port_key)
                if endpoints[0] == endpoints[1]:
                    raise ValueError(f"Verbindung {connection_id!r} darf keinen Port mit sich selbst verbinden")
            _unique_ids(connection_ids, f"Verbindungs-IDs in Molekül {molecule_name!r}")

    def create_molecule(self, molecule_name: str, instance_id: str | None = None) -> Molecule:
        try:
            raw_molecule = self._molecules[molecule_name]
        except KeyError:
            choices = ", ".join(self.molecule_names)
            raise ValueError(f"Unbekanntes Molekül {molecule_name!r}; verfügbar: {choices}") from None

        molecule_id = molecule_name if instance_id is None else instance_id
        _identifier(molecule_id, "instance_id")
        components = tuple(self._parse_component(item) for item in raw_molecule["components"])
        atoms_by_component: dict[str, dict[str, Atom]] = {}
        atoms: list[Atom] = []
        bonds: list[Bond] = []
        lone_pairs: list[ElectronPair] = []
        free_electrons: list[Electron] = []

        for component in components:
            group = self._groups[component.group_name]
            component_atoms: dict[str, Atom] = {}
            for atom_definition in group["atoms"]:
                local_position = _vector(atom_definition["position"], "Atomposition")
                position = _rotate_vector(local_position, component.rotation_degrees)
                position = (
                    position[0] + component.position[0],
                    position[1] + component.position[1],
                    position[2] + component.position[2],
                )
                atom_id = f"{molecule_id}:{component.component_id}:{atom_definition['id']}"
                particles = create_element(
                    atom_definition["element"],
                    atom_id,
                    position,
                    lone_pair_count=atom_definition.get("lone_pairs"),
                    free_electron_directions=tuple(
                        _rotate_vector(_vector(direction, "Elektronenrichtung"), component.rotation_degrees)
                        for direction in atom_definition.get("free_electrons", [])
                    ),
                )
                component_atoms[atom_definition["id"]] = particles.atom
                atoms.append(particles.atom)
                lone_pairs.extend(particles.lone_pairs)
                free_electrons.extend(particles.free_electrons)
                for pair in particles.lone_pairs:
                    for electron in pair.electrons:
                        electron.direction = _rotate_vector(electron.direction, component.rotation_degrees)
                        electron.orbit_normal = _rotate_vector(electron.orbit_normal, component.rotation_degrees)
                for electron in particles.free_electrons:
                    electron.direction = _rotate_vector(electron.direction, component.rotation_degrees)
                    electron.orbit_normal = _rotate_vector(electron.orbit_normal, component.rotation_degrees)
            atoms_by_component[component.component_id] = component_atoms
            for bond_definition in group.get("bonds", []):
                bond_id = f"{molecule_id}:{component.component_id}:{bond_definition['id']}"
                direction_offsets = tuple(
                    _rotate_vector(_vector(offset, "Elektronenoffset"), component.rotation_degrees)
                    for offset in bond_definition.get(
                        "electron_direction_offsets",
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    )
                )
                bond = create_bond(
                    bond_id,
                    component_atoms[bond_definition["from"]],
                    component_atoms[bond_definition["to"]],
                    electron_z_offsets=(0.0, 0.0),
                    electron_direction_offsets=direction_offsets,
                )
                for electron in bond.electrons:
                    electron.orbit_normal = _rotate_vector(
                        electron.orbit_normal,
                        component.rotation_degrees,
                    )
                bonds.append(bond)

        for connection in raw_molecule.get("connections", []):
            first_component = next(
                component for component in components if component.component_id == connection["from"]["component"]
            )
            second_component = next(
                component for component in components if component.component_id == connection["to"]["component"]
            )
            first_group = self._groups[first_component.group_name]
            second_group = self._groups[second_component.group_name]
            first_atom_id = first_group["ports"][connection["from"]["port"]]
            second_atom_id = second_group["ports"][connection["to"]["port"]]
            bonds.append(
                create_bond(
                    f"{molecule_id}:{connection['id']}",
                    atoms_by_component[first_component.component_id][first_atom_id],
                    atoms_by_component[second_component.component_id][second_atom_id],
                )
            )

        return Molecule(
            name=f"{raw_molecule.get('name', molecule_name)} {molecule_id}",
            atoms=tuple(atoms),
            bonds=tuple(bonds),
            lone_pairs=tuple(lone_pairs),
            free_electrons=tuple(free_electrons),
        )

    @staticmethod
    def _parse_component(raw_component: Mapping[str, Any]) -> _Component:
        return _Component(
            component_id=raw_component["id"],
            group_name=raw_component["group"],
            position=_vector(raw_component.get("position", [0.0, 0.0, 0.0]), "Komponentenposition"),
            rotation_degrees=_vector(
                raw_component.get("rotation_degrees", [0.0, 0.0, 0.0]),
                "Komponentenrotation",
            ),
        )
