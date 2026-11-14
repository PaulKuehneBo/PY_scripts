from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from functools import lru_cache
from math import isfinite, sqrt
from random import gauss
from typing import Any

import numpy as np

from .math import (
    Vec3,
    add,
    cross,
    direction_towards,
    dot,
    length,
    remove_radial_velocity,
    scale,
    subtract,
)

ELECTRON_DISPLAY_RADIUS = 0.055
NUCLEUS_REPULSION_FACTOR = 1.0
BOND_MINIMUM_LENGTH = 0.22


@lru_cache(maxsize=32)
def _upper_triangle_indices(particle_count: int) -> tuple[np.ndarray, np.ndarray]:
    return np.triu_indices(particle_count, 1)


@dataclass
class Atom:
    """Ein Atom mit eindeutiger ID und einem Radius für seine Elektronenbahnen."""

    atom_id: str
    element: str
    position: Vec3
    nuclear_charge: float
    core_mass: float
    electron_radius: float
    display_radius: float
    velocity: Vec3 = field(default=(0.0, 0.0, 0.0), init=False)
    visual: Any = field(default=None, init=False, repr=False, compare=False)
    normal_display_radius: float = field(default=0.0, init=False, repr=False)


@dataclass
class NuclearTrace:
    positions: deque[Vec3]
    visual: Any
    marker_visuals: list[tuple[Any, Any]] = field(default_factory=list)
    next_marker_slot: int = 0


@dataclass
class Electron:
    electron_id: str
    atom_id: str
    direction: Vec3
    orbit_direction: float = 1.0
    initial_orbit_speed: float = 0.48
    charge: float = -1.0
    mass: float = 1.0
    bond_id: str | None = None
    orbit_normal: Vec3 = (0.0, 0.0, 1.0)
    display_radius: float = ELECTRON_DISPLAY_RADIUS
    position: Vec3 = field(default=(0.0, 0.0, 0.0), init=False)
    velocity: Vec3 = field(default=(0.0, 0.0, 0.0), init=False)
    visual: Any = field(default=None, init=False, repr=False, compare=False)

    def initialize(self, atom: Atom) -> None:
        if atom.atom_id != self.atom_id:
            raise ValueError(f"Elektron {self.electron_id} gehört nicht zu Atom {atom.atom_id}")

        direction_length = length(self.direction)
        normal_length = length(self.orbit_normal)
        if direction_length == 0 or normal_length == 0:
            raise ValueError(f"Elektron {self.electron_id} hat keine gültige Richtung")
        u = scale(self.direction, 1 / direction_length)
        normal = scale(self.orbit_normal, 1 / normal_length)
        normal_component = dot(normal, u)
        normal = subtract(normal, scale(u, normal_component))
        normal_length = length(normal)
        if normal_length == 0:
            raise ValueError(f"Elektron {self.electron_id} hat eine ungültige Bahnebene")
        normal = scale(normal, 1 / normal_length)
        tangent = scale(cross(normal, u), self.orbit_direction)
        self.position = add(atom.position, scale(u, atom.electron_radius))
        self.velocity = scale(tangent, self.initial_orbit_speed)

    def advance(self, atom: Atom, force: Vec3, dt: float) -> None:
        acceleration = scale(force, 1 / self.mass)
        half_step_velocity = add(self.velocity, scale(acceleration, dt / 2))
        unconstrained_position = add(self.position, scale(half_step_velocity, dt))
        relative_position = subtract(unconstrained_position, atom.position)
        relative_length = length(relative_position)
        if relative_length == 0:
            relative_position = subtract(self.position, atom.position)
            relative_length = length(relative_position)
        if relative_length == 0:
            raise ValueError(f"Elektron {self.electron_id} liegt im Atomkern")

        orbit_direction = scale(relative_position, 1 / relative_length)
        self.position = add(atom.position, scale(orbit_direction, atom.electron_radius))
        relative_velocity = subtract(half_step_velocity, atom.velocity)
        self.velocity = add(
            atom.velocity,
            remove_radial_velocity(relative_velocity, orbit_direction),
        )

    def complete_step(self, atom: Atom, force: Vec3, dt: float) -> None:
        acceleration = scale(force, 1 / self.mass)
        relative_velocity = subtract(self.velocity, atom.velocity)
        velocity = add(relative_velocity, scale(acceleration, dt / 2))
        orbit_direction = subtract(self.position, atom.position)
        orbit_direction = scale(orbit_direction, 1 / length(orbit_direction))
        self.velocity = add(
            atom.velocity,
            remove_radial_velocity(velocity, orbit_direction),
        )


@dataclass
class ElectronPair:
    pair_id: str
    electrons: tuple[Electron, Electron]
    visual: Any = field(default=None, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if len(self.electrons) != 2:
            raise ValueError(f"Elektronenpaar {self.pair_id} muss genau zwei Elektronen enthalten")
        if self.electrons[0].electron_id == self.electrons[1].electron_id:
            raise ValueError(f"Elektronenpaar {self.pair_id} enthält ein Elektron doppelt")


@dataclass
class Bond:
    bond_id: str
    electrons: tuple[Electron, Electron]
    minimum_length: float = BOND_MINIMUM_LENGTH
    visual: Any = field(default=None, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        ElectronPair(self.bond_id, self.electrons)
        if self.minimum_length <= 0:
            raise ValueError("Die Mindestlänge einer Bindung muss positiv sein")
        if self.electrons[0].atom_id == self.electrons[1].atom_id:
            raise ValueError(f"Bindung {self.bond_id} muss Elektronen aus zwei Atomen verbinden")
        if any(electron.bond_id != self.bond_id for electron in self.electrons):
            raise ValueError(f"Beide Elektronen müssen auf Bindung {self.bond_id} verweisen")


def create_bond(
    bond_id: str,
    first_atom: Atom,
    second_atom: Atom,
    *,
    electron_ids: tuple[str, str] | None = None,
    electron_z_offsets: tuple[float, float] = (-0.12, 0.12),
    orbit_directions: tuple[float, float] = (1.0, -1.0),
    electron_direction_offsets: tuple[Vec3, Vec3] = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
) -> Bond:
    """Erzeugt eine Bindung mit je einem Elektron an beiden Atomkernen."""
    first_electron_id, second_electron_id = (
        (f"{bond_id}:e1", f"{bond_id}:e2")
        if electron_ids is None
        else electron_ids
    )
    first_z_offset, second_z_offset = electron_z_offsets
    first_orbit_direction, second_orbit_direction = orbit_directions
    first_direction_offset, second_direction_offset = electron_direction_offsets

    return Bond(
        bond_id,
        (
            Electron(
                first_electron_id,
                first_atom.atom_id,
                add(
                    direction_towards(
                        first_atom.position,
                        second_atom.position,
                        first_z_offset,
                    ),
                    first_direction_offset,
                ),
                bond_id=bond_id,
                orbit_direction=first_orbit_direction,
            ),
            Electron(
                second_electron_id,
                second_atom.atom_id,
                add(
                    direction_towards(
                        second_atom.position,
                        first_atom.position,
                        second_z_offset,
                    ),
                    second_direction_offset,
                ),
                bond_id=bond_id,
                orbit_direction=second_orbit_direction,
            ),
        ),
    )


@dataclass
class ParticleState:
    particle_id: str
    atom_id: str
    position: Vec3
    charge: float
    particle_type: str
    mass: float
    bond_id: str | None = None
    pair_id: str | None = None
    molecule_id: int | None = None
    nuclear_radius: float | None = None
    parent_nucleus_position: Vec3 | None = None
    parent_nuclear_radius: float | None = None
    electron_radius: float | None = None
    bond_minimum_length: float = BOND_MINIMUM_LENGTH


@dataclass(frozen=True)
class ForceTotals:
    particle_forces: dict[str, Vec3]
    nucleus_repulsions: dict[str, Vec3]
    nucleus_electron_forces: dict[str, Vec3]


def _nuclear_radius_force_factor(
    distance: float,
    first_radius: float | None,
    second_radius: float | None = None,
) -> float:
    if first_radius is None:
        return 1.0
    if first_radius <= 0 or (second_radius is not None and second_radius <= 0):
        raise ValueError("Kernradien müssen positiv sein")
    combined_radius = first_radius + (second_radius or 0.0)
    return min(1.0, distance / (2 * combined_radius))


def _inter_molecular_nuclear_force_factor(
    distance: float,
    bond_minimum_length: float,
    target_electron_radius: float | None,
) -> float:
    if bond_minimum_length <= 0:
        raise ValueError("Die Mindestlänge der Bindung muss positiv sein")
    if target_electron_radius is None:
        return 1.0
    if target_electron_radius <= 0:
        raise ValueError("Elektronenbahnradius muss positiv sein")
    cutoff_radius = bond_minimum_length + target_electron_radius
    if distance >= cutoff_radius:
        return 1.0
    return (distance / cutoff_radius) ** 3


def _nuclear_electron_coulomb_force(
    nucleus: ParticleState,
    electron: ParticleState,
    strength: float,
    softening_squared: float,
    *,
    inter_molecular: bool = False,
) -> Vec3:
    displacement = subtract(nucleus.position, electron.position)
    distance_squared = dot(displacement, displacement)
    distance = sqrt(distance_squared)
    if electron.parent_nucleus_position is None:
        radius_distance = distance
    else:
        radius_distance = length(subtract(nucleus.position, electron.parent_nucleus_position))
    softened_distance_squared = distance_squared + softening_squared
    radius_factor = (
        _inter_molecular_nuclear_force_factor(
            distance,
            electron.bond_minimum_length,
            nucleus.electron_radius,
        )
        if inter_molecular
        else _nuclear_radius_force_factor(
            radius_distance,
            nucleus.nuclear_radius,
            electron.parent_nuclear_radius,
        )
    )
    scale_value = (
        strength
        * nucleus.charge
        * electron.charge
        / (softened_distance_squared * sqrt(softened_distance_squared))
        * radius_factor
    )
    return scale(displacement, scale_value)


@dataclass(frozen=True)
class CoulombForceModel:
    """Coulomb-Wechselwirkungen in skalierten Visualisierungseinheiten.

    Der Faktor für intermolekulare Kernabstoßung skaliert nur Kern-Kern-Kräfte
    zwischen verschiedenen Molekülen; intramolekulare Kräfte bleiben unverändert.
    Ein optionaler Radius lässt diese zusätzliche Kraft weich auf null abklingen.
    """

    strength: float = 0.02
    softening: float = 0.12
    pair_target_length: float = 0.22
    pair_length_stiffness: float = 1.0
    inter_molecular_nucleus_repulsion_factor: float = 1.0
    inter_molecular_nucleus_repulsion_radius: float | None = None

    def __post_init__(self) -> None:
        if self.strength < 0:
            raise ValueError("Die Coulomb-Stärke darf nicht negativ sein")
        if self.softening <= 0:
            raise ValueError("Der Softening-Abstand muss positiv sein")
        if self.pair_target_length <= 0:
            raise ValueError("Die Ziel-Paarlänge muss positiv sein")
        if self.pair_length_stiffness < 0:
            raise ValueError("Die Steifigkeit der Paarlänge darf nicht negativ sein")
        if (
            not isfinite(self.inter_molecular_nucleus_repulsion_factor)
            or self.inter_molecular_nucleus_repulsion_factor < 0
        ):
            raise ValueError("Der Faktor der intermolekularen Kernabstoßung muss endlich und nicht negativ sein")
        if (
            self.inter_molecular_nucleus_repulsion_radius is not None
            and (
                not isfinite(self.inter_molecular_nucleus_repulsion_radius)
                or self.inter_molecular_nucleus_repulsion_radius <= 0
            )
        ):
            raise ValueError("Der Radius der intermolekularen Kernabstoßung muss positiv sein")

    def calculate(
        self,
        particles: tuple[ParticleState, ...],
        *,
        inter_molecular_bond_forces_enabled: bool = True,
    ) -> ForceTotals:
        particle_count = len(particles)
        particle_ids = tuple(particle.particle_id for particle in particles)
        nucleus_indices = np.fromiter(
            (index for index, particle in enumerate(particles) if particle.particle_type == "nucleus"),
            dtype=np.intp,
        )
        nucleus_atom_ids = tuple(particles[index].atom_id for index in nucleus_indices)
        if particle_count < 2:
            zero_forces = {particle_id: (0.0, 0.0, 0.0) for particle_id in particle_ids}
            zero_nucleus_forces = {atom_id: (0.0, 0.0, 0.0) for atom_id in nucleus_atom_ids}
            return ForceTotals(zero_forces, zero_nucleus_forces, zero_nucleus_forces.copy())

        positions = np.asarray([particle.position for particle in particles], dtype=np.float64)
        charges = np.asarray([particle.charge for particle in particles], dtype=np.float64)
        is_electron = np.fromiter(
            (particle.particle_type == "electron" for particle in particles),
            dtype=np.bool_,
            count=particle_count,
        )
        atom_ids = tuple(particle.atom_id for particle in particles)
        molecule_ids = tuple(particle.molecule_id for particle in particles)
        bond_ids = tuple(particle.bond_id for particle in particles)
        pair_ids = tuple(particle.pair_id for particle in particles)
        atom_index_by_id = {atom_id: index for index, atom_id in enumerate(atom_ids)}
        atom_indices = np.fromiter((atom_index_by_id[atom_id] for atom_id in atom_ids), dtype=np.intp)
        molecule_indices = np.fromiter(
            (-1 if molecule_id is None else molecule_id for molecule_id in molecule_ids),
            dtype=np.int64,
        )
        nuclear_radii = np.asarray(
            [np.nan if particle.nuclear_radius is None else particle.nuclear_radius for particle in particles],
            dtype=np.float64,
        )
        electron_radii = np.asarray(
            [np.nan if particle.electron_radius is None else particle.electron_radius for particle in particles],
            dtype=np.float64,
        )
        parent_nuclear_radii = np.asarray(
            [
                np.nan if particle.parent_nuclear_radius is None else particle.parent_nuclear_radius
                for particle in particles
            ],
            dtype=np.float64,
        )
        parent_nucleus_positions = np.asarray(
            [
                particle.position
                if particle.parent_nucleus_position is None
                else particle.parent_nucleus_position
                for particle in particles
            ],
            dtype=np.float64,
        )
        bond_minimum_lengths = np.fromiter(
            (particle.bond_minimum_length for particle in particles),
            dtype=np.float64,
            count=particle_count,
        )

        first_indices, second_indices = _upper_triangle_indices(particle_count)
        first_is_electron = is_electron[first_indices]
        second_is_electron = is_electron[second_indices]
        is_electron_nucleus = first_is_electron != second_is_electron
        is_nucleus_nucleus = ~first_is_electron & ~second_is_electron
        first_molecule_indices = molecule_indices[first_indices]
        second_molecule_indices = molecule_indices[second_indices]
        is_inter_molecular = (
            (first_molecule_indices >= 0)
            & (second_molecule_indices >= 0)
            & (first_molecule_indices != second_molecule_indices)
        )
        electron_indices = np.where(first_is_electron, first_indices, second_indices)
        electron_is_bonded = np.fromiter(
            (bond_ids[index] is not None for index in electron_indices),
            dtype=np.bool_,
            count=len(first_indices),
        )
        same_atom = atom_indices[first_indices] == atom_indices[second_indices]
        skip_electron_nucleus = is_electron_nucleus & (
            (~is_inter_molecular & (electron_is_bonded | same_atom))
            | (is_inter_molecular & electron_is_bonded & (not inter_molecular_bond_forces_enabled))
        )
        active_pairs = ~skip_electron_nucleus

        displacements = positions[first_indices] - positions[second_indices]
        distances_squared = np.einsum("ij,ij->i", displacements, displacements)
        softening_squared = self.softening * self.softening
        softened_distances_squared = distances_squared + softening_squared
        charge_products = charges[first_indices] * charges[second_indices]
        force_scales = self.strength * charge_products / (
            softened_distances_squared * np.sqrt(softened_distances_squared)
        )
        inter_molecular_nucleus_pairs = is_nucleus_nucleus & is_inter_molecular
        force_scales[inter_molecular_nucleus_pairs] *= self.inter_molecular_nucleus_repulsion_factor
        if self.inter_molecular_nucleus_repulsion_radius is not None:
            pair_distances = np.sqrt(distances_squared[inter_molecular_nucleus_pairs])
            normalized_distances = np.clip(
                pair_distances / self.inter_molecular_nucleus_repulsion_radius,
                0.0,
                1.0,
            )
            attenuation = (1.0 - normalized_distances) ** 2 * (
                1.0 + 2.0 * normalized_distances
            )
            force_scales[inter_molecular_nucleus_pairs] *= attenuation

        electron_nucleus_pair_indices = np.flatnonzero(is_electron_nucleus & active_pairs)
        radius_factors = np.ones(len(first_indices), dtype=np.float64)
        if electron_nucleus_pair_indices.size:
            target_indices = np.where(
                first_is_electron[electron_nucleus_pair_indices],
                second_indices[electron_nucleus_pair_indices],
                first_indices[electron_nucleus_pair_indices],
            )
            source_electron_indices = electron_indices[electron_nucleus_pair_indices]
            inter_pair_mask = is_inter_molecular[electron_nucleus_pair_indices]

            if np.any(inter_pair_mask):
                pair_indices = electron_nucleus_pair_indices[inter_pair_mask]
                target_indices_for_pairs = target_indices[inter_pair_mask]
                source_indices_for_pairs = source_electron_indices[inter_pair_mask]
                bond_minimum_lengths_for_pairs = bond_minimum_lengths[source_indices_for_pairs]
                target_electron_radii = electron_radii[target_indices_for_pairs]
                if np.any(bond_minimum_lengths_for_pairs <= 0):
                    raise ValueError("Die Mindestlänge der Bindung muss positiv sein")
                has_target_radius = ~np.isnan(target_electron_radii)
                if np.any(target_electron_radii[has_target_radius] <= 0):
                    raise ValueError("Elektronenbahnradius muss positiv sein")

                cutoff_radii = bond_minimum_lengths_for_pairs + np.where(
                    has_target_radius,
                    target_electron_radii,
                    0.0,
                )
                pair_distances = np.sqrt(distances_squared[pair_indices])
                taper_mask = has_target_radius & (pair_distances < cutoff_radii)
                factors = np.ones(len(pair_indices), dtype=np.float64)
                factors[taper_mask] = (
                    pair_distances[taper_mask] / cutoff_radii[taper_mask]
                ) ** 3
                radius_factors[pair_indices] = factors

            if np.any(~inter_pair_mask):
                pair_indices = electron_nucleus_pair_indices[~inter_pair_mask]
                target_indices_for_pairs = target_indices[~inter_pair_mask]
                source_indices_for_pairs = source_electron_indices[~inter_pair_mask]
                target_nuclear_radii_for_pairs = nuclear_radii[target_indices_for_pairs]
                parent_nuclear_radii_for_pairs = parent_nuclear_radii[source_indices_for_pairs]
                has_target_radius = ~np.isnan(target_nuclear_radii_for_pairs)
                if np.any(target_nuclear_radii_for_pairs[has_target_radius] <= 0) or np.any(
                    parent_nuclear_radii_for_pairs[
                        has_target_radius & ~np.isnan(parent_nuclear_radii_for_pairs)
                    ]
                    <= 0
                ):
                    raise ValueError("Kernradien müssen positiv sein")

                source_parent_positions = parent_nucleus_positions[source_indices_for_pairs]
                radius_displacements = positions[target_indices_for_pairs] - source_parent_positions
                radius_distances = np.sqrt(
                    np.einsum("ij,ij->i", radius_displacements, radius_displacements)
                )
                combined_radii = target_nuclear_radii_for_pairs + np.where(
                    np.isnan(parent_nuclear_radii_for_pairs),
                    0.0,
                    parent_nuclear_radii_for_pairs,
                )
                factors = np.ones(len(pair_indices), dtype=np.float64)
                factors[has_target_radius] = np.minimum(
                    1.0,
                    radius_distances[has_target_radius] / (2 * combined_radii[has_target_radius]),
                )
                radius_factors[pair_indices] = factors

        force_scales[is_electron_nucleus] *= radius_factors[is_electron_nucleus]
        force_scales[~active_pairs] = 0.0
        pair_forces = displacements * force_scales[:, np.newaxis]

        particle_force_values = np.zeros((particle_count, 3), dtype=np.float64)
        np.add.at(particle_force_values, first_indices, pair_forces)
        np.add.at(particle_force_values, second_indices, -pair_forces)

        nucleus_repulsion_values = np.zeros((particle_count, 3), dtype=np.float64)
        nucleus_nucleus_pair_indices = np.flatnonzero(is_nucleus_nucleus)
        np.add.at(
            nucleus_repulsion_values,
            first_indices[nucleus_nucleus_pair_indices],
            pair_forces[nucleus_nucleus_pair_indices],
        )
        np.add.at(
            nucleus_repulsion_values,
            second_indices[nucleus_nucleus_pair_indices],
            -pair_forces[nucleus_nucleus_pair_indices],
        )

        nucleus_electron_values = np.zeros((particle_count, 3), dtype=np.float64)
        if electron_nucleus_pair_indices.size:
            target_nucleus_indices = target_indices
            target_forces = np.where(
                first_is_electron[electron_nucleus_pair_indices, np.newaxis],
                -pair_forces[electron_nucleus_pair_indices],
                pair_forces[electron_nucleus_pair_indices],
            )
            np.add.at(nucleus_electron_values, target_nucleus_indices, target_forces)
            nucleus_indices_by_atom = {
                particles[index].atom_id: index for index in nucleus_indices
            }
            owner_nucleus_indices = np.fromiter(
                (nucleus_indices_by_atom[atom_ids[index]] for index in source_electron_indices),
                dtype=np.intp,
                count=len(source_electron_indices),
            )
            np.add.at(nucleus_electron_values, owner_nucleus_indices, -target_forces)

        paired_electrons: dict[str, list[int]] = {}
        for index, pair_id in enumerate(pair_ids):
            if is_electron[index] and pair_id is not None:
                paired_electrons.setdefault(pair_id, []).append(index)
        for pair_indices in paired_electrons.values():
            if len(pair_indices) != 2:
                raise ValueError("Eine Bindung oder ein Elektronenpaar muss genau zwei Elektronen enthalten")

        if paired_electrons:
            pair_first_indices = np.fromiter(
                (indices[0] for indices in paired_electrons.values()),
                dtype=np.intp,
            )
            pair_second_indices = np.fromiter(
                (indices[1] for indices in paired_electrons.values()),
                dtype=np.intp,
            )
            pair_displacements = positions[pair_first_indices] - positions[pair_second_indices]
            pair_distances = np.linalg.norm(pair_displacements, axis=1)
            nonzero_pairs = pair_distances > 0
            pair_first_indices = pair_first_indices[nonzero_pairs]
            pair_second_indices = pair_second_indices[nonzero_pairs]
            pair_displacements = pair_displacements[nonzero_pairs]
            pair_distances = pair_distances[nonzero_pairs]
            deviations = pair_distances - self.pair_target_length
            pair_force_values = (
                -self.pair_length_stiffness * deviations / pair_distances
            )[:, np.newaxis] * pair_displacements
            np.add.at(particle_force_values, pair_first_indices, pair_force_values)
            np.add.at(particle_force_values, pair_second_indices, -pair_force_values)

        particle_forces = {
            particle_id: tuple(float(component) for component in particle_force_values[index])
            for index, particle_id in enumerate(particle_ids)
        }
        nucleus_repulsions = {
            particles[index].atom_id: tuple(float(component) for component in nucleus_repulsion_values[index])
            for index in nucleus_indices
        }
        nucleus_electron_forces = {
            particles[index].atom_id: tuple(float(component) for component in nucleus_electron_values[index])
            for index in nucleus_indices
        }
        return ForceTotals(particle_forces, nucleus_repulsions, nucleus_electron_forces)


@dataclass
class Molecule:
    name: str
    atoms: tuple[Atom, ...]
    bonds: tuple[Bond, ...]
    lone_pairs: tuple[ElectronPair, ...]
    free_electrons: tuple[Electron, ...] = ()

    def __post_init__(self) -> None:
        atom_ids = [atom.atom_id for atom in self.atoms]
        if len(atom_ids) != len(set(atom_ids)):
            raise ValueError("Jedes Atom benötigt eine eindeutige ID")

        bond_ids = [bond.bond_id for bond in self.bonds]
        if len(bond_ids) != len(set(bond_ids)):
            raise ValueError("Jede Bindung benötigt eine eindeutige ID")

        known_atom_ids = set(atom_ids)
        known_bond_ids = set(bond_ids)
        bonded_electrons = [electron for bond in self.bonds for electron in bond.electrons]
        lone_pair_electrons = [electron for pair in self.lone_pairs for electron in pair.electrons]
        all_electrons = (*bonded_electrons, *lone_pair_electrons, *self.free_electrons)

        electron_ids = [electron.electron_id for electron in all_electrons]
        if len(electron_ids) != len(set(electron_ids)):
            raise ValueError("Ein Elektron darf nur einmal im Molekül vorkommen")

        for electron in all_electrons:
            if electron.atom_id not in known_atom_ids:
                raise ValueError(
                    f"Elektron {electron.electron_id} verweist auf unbekanntes Atom " f"{electron.atom_id}"
                )
            if electron.bond_id is not None and electron.bond_id not in known_bond_ids:
                raise ValueError(
                    f"Elektron {electron.electron_id} verweist auf unbekannte Bindung " f"{electron.bond_id}"
                )

        for pair in self.lone_pairs:
            if any(electron.bond_id is not None for electron in pair.electrons):
                raise ValueError(f"Freies Elektronenpaar {pair.pair_id} darf nicht gebunden sein")
            if pair.electrons[0].atom_id != pair.electrons[1].atom_id:
                raise ValueError(f"Freies Elektronenpaar {pair.pair_id} muss zu einem Atom gehören")

        if any(electron.bond_id is not None for electron in self.free_electrons):
            raise ValueError("Ein freies Elektron darf keine Bindungs-ID haben")

        if len(bonded_electrons) != sum(electron.bond_id is not None for electron in all_electrons):
            raise ValueError("Bindungs-IDs und Elektronen der Bindungen stimmen nicht überein")


def _interaction_type(first: ParticleState, second: ParticleState) -> str:
    first_is_electron = first.particle_type == "electron"
    second_is_electron = second.particle_type == "electron"
    if first_is_electron and second_is_electron:
        return "electron-electron"
    if not first_is_electron and not second_is_electron:
        return "nucleus-nucleus"
    return "electron-nucleus"


# Geschwindigkeitserhalt pro Referenzintervall: kleiner = mehr Dämpfung,
# größer = längere Bewegung; 1.0 bedeutet keine Dämpfung.
DAMPING_FACTOR = 0.99
# Zeitspanne, über die DAMPING_FACTOR angewendet wird (1/60 s = 60 FPS).
# DAMPING_REFERENCE_DT = 1 / 60
DAMPING_REFERENCE_DT = 1 / 20
# Dauer eines Physik-Teilschritts: kleiner = genauer/stabiler, aber rechenintensiver;
# größer = schnellere, aber gröbere und potenziell instabilere Bewegung.
# SIMULATION_DT = 1 / 240
SIMULATION_DT = 1 / 20
# Anzahl Physik-Teilschritte pro Bild: größer = mehr Simulation pro Bild
# und höhere Rechenlast. Zusammen mit SIMULATION_DT bestimmt dies die Simulationszeit pro Bild.
SIMULATION_STEPS_PER_FRAME = 4
# Standardabweichung eines zufälligen Geschwindigkeitsstoßes pro Raumachse
# und DAMPING_REFERENCE_DT. Kleiner = ruhigere Bewegung, größer = stärkere
# Brownsche Bewegung. 0.0 deaktiviert die zufälligen Stöße.
# BROWNIAN_VELOCITY_AMPLITUDE = 0.02
BROWNIAN_VELOCITY_AMPLITUDE = 0.0015
# Separater Faktor für die Brownsche Geschwindigkeitsstreuung der Elektronen
# in freien Elektronenpaaren. 2.5 bedeutet 2.5-mal die Standardstreuung.
LONE_PAIR_BROWNIAN_VELOCITY_FACTOR = 2.5
# Anzahl der zuletzt gespeicherten Kernpositionen je Atom. Bei 60 Bildern/s
# entsprechen 1000 Punkte etwa 16,7 Sekunden.
TRACE_POINT_LIMIT = 1000
# Markiere die Kernposition einmal pro simulierter Sekunde mit einem Kreuz.
TRACE_MARKER_INTERVAL_SECONDS = 1.0

force_model = CoulombForceModel()


def _default_force_model() -> CoulombForceModel:
    return force_model


class MoleculeCluster:
    """Molekülcluster mit konfigurierbarem Kraftmodell."""

    def __init__(
        self,
        molecules: Molecule | tuple[Molecule, ...],
        *,
        force_model: CoulombForceModel | None = None,
    ) -> None:
        self.molecules = (molecules,) if isinstance(molecules, Molecule) else tuple(molecules)
        if not all(isinstance(molecule, Molecule) for molecule in self.molecules):
            raise ValueError("Mindestens ein Molekül muss simuliert werden")
        self.force_model = _default_force_model() if force_model is None else force_model

        self._electrons = tuple(
            electron
            for item in self.molecules
            for electron in (
                *(electron for bond in item.bonds for electron in bond.electrons),
                *(electron for pair in item.lone_pairs for electron in pair.electrons),
                *item.free_electrons,
            )
        )
        self._atoms = tuple(atom for molecule in self.molecules for atom in molecule.atoms)
        self._atoms_by_id = {atom.atom_id: atom for atom in self._atoms}

        self._bonds = tuple(bond for molecule in self.molecules for bond in molecule.bonds)

        self._lone_pairs = tuple(pair for molecule in self.molecules for pair in molecule.lone_pairs)

        self._electron_parent_atoms = tuple(self._atoms_by_id[electron.atom_id] for electron in self._electrons)
        self._particle_states = self._calculate_particle_states()
        self._nucleus_particle_states = self._particle_states[: len(self._atoms)]
        self._electron_particle_states = self._particle_states[len(self._atoms) :]

    def _calculate_particle_states(self) -> tuple[ParticleState, ...]:
        molecule_ids_by_atom = {atom.atom_id: id(item) for item in self.molecules for atom in item.atoms}
        pair_ids_by_electron = {
            electron.electron_id: bond.bond_id
            for item in self.molecules
            for bond in item.bonds
            for electron in bond.electrons
        }
        pair_ids_by_electron.update(
            {
                electron.electron_id: pair.pair_id
                for item in self.molecules
                for pair in item.lone_pairs
                for electron in pair.electrons
            }
        )
        particles = [
            ParticleState(
                f"{atom.atom_id}:nucleus",
                atom.atom_id,
                atom.position,
                atom.nuclear_charge,
                "nucleus",
                atom.core_mass,
                None,
                None,
                id(item),
                atom.display_radius,
                None,
                None,
                atom.electron_radius,
            )
            for item in self.molecules
            for atom in item.atoms
        ]
        bond_minimum_lengths_by_electron_id: dict[str, float] = {}
        for item in self.molecules:
            for bond in item.bonds:
                for electron in bond.electrons:
                    bond_minimum_lengths_by_electron_id[electron.electron_id] = bond.minimum_length

        particles.extend(
            ParticleState(
                electron.electron_id,
                electron.atom_id,
                electron.position,
                electron.charge,
                "electron",
                electron.mass,
                electron.bond_id,
                pair_ids_by_electron.get(electron.electron_id),
                molecule_ids_by_atom[electron.atom_id],
                None,
                self._atoms_by_id[electron.atom_id].position,
                self._atoms_by_id[electron.atom_id].display_radius,
                None,
                bond_minimum_lengths_by_electron_id.get(
                    electron.electron_id,
                    BOND_MINIMUM_LENGTH,
                ),
            )
            for electron in self.electrons
        )
        particle_ids = [particle.particle_id for particle in particles]
        if len(particle_ids) != len(set(particle_ids)):
            raise ValueError("Partikel-IDs müssen über alle Moleküle hinweg eindeutig sein")
        return tuple(particles)

    @property
    def electrons(self) -> tuple[Electron, ...]:
        return self._electrons

    @property
    def atoms(self) -> tuple[Atom, ...]:
        return self._atoms

    @property
    def bonds(self) -> tuple[Bond, ...]:
        return self._bonds

    @property
    def lone_pairs(self) -> tuple[ElectronPair, ...]:
        return self._lone_pairs

    @property
    def particle_states(self) -> tuple[ParticleState, ...]:
        for state, atom in zip(self._nucleus_particle_states, self._atoms):
            state.position = atom.position

        for state, electron, parent_atom in zip(
            self._electron_particle_states,
            self._electrons,
            self._electron_parent_atoms,
        ):
            state.position = electron.position
            state.parent_nucleus_position = parent_atom.position
        return self._particle_states


def _brownian_velocity_kick(dt: float, velocity_factor: float = 1.0) -> Vec3:
    if BROWNIAN_VELOCITY_AMPLITUDE < 0:
        raise ValueError("Die Brownsche Bewegungsamplitude darf nicht negativ sein")
    if DAMPING_REFERENCE_DT <= 0:
        raise ValueError("Das Dämpfungs-Referenzintervall muss positiv sein")
    if dt < 0:
        raise ValueError("Der Simulationszeitschritt darf nicht negativ sein")
    if velocity_factor < 0:
        raise ValueError("Der Brownsche Geschwindigkeitsfaktor darf nicht negativ sein")

    scale = BROWNIAN_VELOCITY_AMPLITUDE * velocity_factor * sqrt(dt / DAMPING_REFERENCE_DT)
    return gauss(0.0, scale), gauss(0.0, scale), gauss(0.0, scale)


def recenter_system(
    molecules: MoleculeCluster,
    remove_com_velocity: bool,
) -> None:
    atoms = molecules.atoms
    electrons = molecules.electrons
    total_mass = sum(atom.core_mass for atom in atoms) + sum(electron.mass for electron in electrons)
    if total_mass <= 0:
        raise ValueError("Die Gesamtmasse des Molekülsystems muss positiv sein")

    center_of_mass = (0.0, 0.0, 0.0)
    center_velocity = (0.0, 0.0, 0.0)
    for atom in atoms:
        mass_fraction = atom.core_mass / total_mass
        center_of_mass = add(center_of_mass, scale(atom.position, mass_fraction))
        center_velocity = add(center_velocity, scale(atom.velocity, mass_fraction))
    for electron in electrons:
        mass_fraction = electron.mass / total_mass
        center_of_mass = add(center_of_mass, scale(electron.position, mass_fraction))
        center_velocity = add(center_velocity, scale(electron.velocity, mass_fraction))

    for atom in atoms:
        atom.position = subtract(atom.position, center_of_mass)
        if remove_com_velocity:
            atom.velocity = subtract(atom.velocity, center_velocity)
    for electron in electrons:
        electron.position = subtract(electron.position, center_of_mass)
        if remove_com_velocity:
            electron.velocity = subtract(electron.velocity, center_velocity)


def integrate_electrons(
    molecules: MoleculeCluster,
    dt: float,
    *,
    inter_molecular_bond_forces_enabled: bool = True,
) -> None:
    atoms = molecules.atoms
    electrons = molecules.electrons
    lone_pair_electron_ids = {
        electron.electron_id for item in molecules.molecules for pair in item.lone_pairs for electron in pair.electrons
    }
    atoms_by_id = {atom.atom_id: atom for atom in atoms}
    reactions_by_atom = {atom.atom_id: (0.0, 0.0, 0.0) for atom in atoms}
    current_totals = molecules.force_model.calculate(
        molecules.particle_states,
        inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
    )
    current_forces = current_totals.particle_forces

    for electron in electrons:
        atom = atoms_by_id[electron.atom_id]
        if electron.bond_id is not None:
            reactions_by_atom[atom.atom_id] = add(
                reactions_by_atom[atom.atom_id],
                current_forces[electron.electron_id],
            )

    old_atom_positions = {atom.atom_id: atom.position for atom in atoms}
    old_atom_velocities = {atom.atom_id: atom.velocity for atom in atoms}
    for atom in atoms:
        applied_force = add(
            add(
                reactions_by_atom[atom.atom_id],
                scale(
                    current_totals.nucleus_repulsions[atom.atom_id],
                    NUCLEUS_REPULSION_FACTOR,
                ),
            ),
            current_totals.nucleus_electron_forces[atom.atom_id],
        )
        atom.velocity = scale(
            add(atom.velocity, scale(applied_force, dt / atom.core_mass)),
            DAMPING_FACTOR ** (dt / DAMPING_REFERENCE_DT),
        )
        atom.velocity = add(atom.velocity, _brownian_velocity_kick(dt))
        atom.position = add(atom.position, scale(atom.velocity, dt))

    for electron in electrons:
        atom = atoms_by_id[electron.atom_id]
        atom_displacement = subtract(atom.position, old_atom_positions[atom.atom_id])
        atom_velocity_change = subtract(atom.velocity, old_atom_velocities[atom.atom_id])
        electron.position = add(electron.position, atom_displacement)
        electron.velocity = add(electron.velocity, atom_velocity_change)

    updated_forces = molecules.force_model.calculate(
        molecules.particle_states,
        inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
    ).particle_forces
    for electron in electrons:
        atom = atoms_by_id[electron.atom_id]
        velocity_factor = LONE_PAIR_BROWNIAN_VELOCITY_FACTOR if electron.electron_id in lone_pair_electron_ids else 1.0
        electron.velocity = add(
            electron.velocity,
            _brownian_velocity_kick(dt, velocity_factor),
        )
        electron.advance(atom, updated_forces[electron.electron_id], dt)

    next_forces = molecules.force_model.calculate(
        molecules.particle_states,
        inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
    ).particle_forces
    for electron in electrons:
        atom = atoms_by_id[electron.atom_id]
        electron.complete_step(atom, next_forces[electron.electron_id], dt)
        electron.velocity = scale(
            electron.velocity,
            DAMPING_FACTOR ** (dt / DAMPING_REFERENCE_DT),
        )

    recenter_system(molecules, remove_com_velocity=True)
