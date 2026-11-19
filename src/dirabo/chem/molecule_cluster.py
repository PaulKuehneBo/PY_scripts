from __future__ import annotations

from math import cos, isfinite, pi, sin, sqrt
from random import Random, gauss

from .args import CHEM_RENDER_ARGS
from .math import Vec3, add, length, scale, subtract
from .particles import (
    BOND_MINIMUM_LENGTH,
    BROWNIAN_VELOCITY_AMPLITUDE,
    DAMPING_FACTOR,
    DAMPING_REFERENCE_DT,
    LONE_PAIR_BROWNIAN_VELOCITY_FACTOR,
    NUCLEUS_REPULSION_FACTOR,
    Atom,
    Bond,
    CoulombForceModel,
    Electron,
    ElectronPair,
    Molecule,
    ParticleState,
)

PLACEMENT_RADIUS_FACTOR = 3.0
PLACEMENT_SEARCH_LIMIT = 256
PLACEMENT_SEARCH_TOLERANCE = 1e-5


def _default_force_model() -> CoulombForceModel:
    return CoulombForceModel(
        inter_molecular_nucleus_repulsion_factor=CHEM_RENDER_ARGS.nucleus_repulsion_factor,
        inter_molecular_nucleus_repulsion_radius=CHEM_RENDER_ARGS.nucleus_repulsion_radius,
    )


def _rotate_vector(vector: Vec3, pitch: float, yaw: float) -> Vec3:
    pitch_cos = cos(pitch)
    pitch_sin = sin(pitch)
    pitched = (
        pitch_cos * vector[0] + pitch_sin * vector[2],
        vector[1],
        -pitch_sin * vector[0] + pitch_cos * vector[2],
    )
    yaw_cos = cos(yaw)
    yaw_sin = sin(yaw)
    return (
        yaw_cos * pitched[0] - yaw_sin * pitched[1],
        yaw_sin * pitched[0] + yaw_cos * pitched[1],
        pitched[2],
    )


class MoleculeCluster:
    """Molekülcluster mit Geometrie, Platzierung und Simulationsmethoden."""

    def __init__(
        self,
        molecules: Molecule | tuple[Molecule, ...],
        *,
        force_model: CoulombForceModel | None = None,
        rng: Random | None = None,
    ) -> None:
        self.molecules = (molecules,) if isinstance(molecules, Molecule) else tuple(molecules)
        if not self.molecules or not all(isinstance(molecule, Molecule) for molecule in self.molecules):
            raise ValueError("Mindestens ein Molekül muss simuliert werden")
        self.force_model = _default_force_model() if force_model is None else force_model
        self._rng = Random() if rng is None else rng
        self._refresh_particles()

    @staticmethod
    def _validate_unique_ids(molecules: tuple[Molecule, ...]) -> None:
        atom_ids = [atom.atom_id for molecule in molecules for atom in molecule.atoms]
        if len(atom_ids) != len(set(atom_ids)):
            raise ValueError("Atom-IDs müssen über alle Moleküle hinweg eindeutig sein")

        electron_ids = [
            electron.electron_id
            for molecule in molecules
            for electron in (
                *(electron for bond in molecule.bonds for electron in bond.electrons),
                *(electron for pair in molecule.lone_pairs for electron in pair.electrons),
                *molecule.free_electrons,
            )
        ]
        if len(electron_ids) != len(set(electron_ids)):
            raise ValueError("Elektron-IDs müssen über alle Moleküle hinweg eindeutig sein")

        particle_ids = [*(f"{atom_id}:nucleus" for atom_id in atom_ids), *electron_ids]
        if len(particle_ids) != len(set(particle_ids)):
            raise ValueError("Partikel-IDs müssen über alle Moleküle hinweg eindeutig sein")

        pair_ids = [
            pair_id
            for molecule in molecules
            for pair_id in (
                *(bond.bond_id for bond in molecule.bonds),
                *(pair.pair_id for pair in molecule.lone_pairs),
            )
        ]
        if len(pair_ids) != len(set(pair_ids)):
            raise ValueError("Bindungs- und Elektronenpaar-IDs müssen eindeutig sein")

    def _refresh_particles(self) -> None:
        self._validate_unique_ids(self.molecules)
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
                nucleus_repulsion_radius=atom.nucleus_repulsion_radius,
            )
            for item in self.molecules
            for atom in item.atoms
        ]
        bond_minimum_lengths_by_electron_id = {
            electron.electron_id: bond.minimum_length
            for item in self.molecules
            for bond in item.bonds
            for electron in bond.electrons
        }
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
        for state, atom in zip(self._nucleus_particle_states, self.atoms):
            state.position = atom.position
            state.nucleus_repulsion_radius = atom.nucleus_repulsion_radius

        for state, electron, parent_atom in zip(
            self._electron_particle_states,
            self.electrons,
            self._electron_parent_atoms,
        ):
            state.position = electron.position
            state.parent_nucleus_position = parent_atom.position
        return self._particle_states

    @property
    def center(self) -> Vec3:
        return (
            sum(molecule.center[0] for molecule in self.molecules) / len(self.molecules),
            sum(molecule.center[1] for molecule in self.molecules) / len(self.molecules),
            sum(molecule.center[2] for molecule in self.molecules) / len(self.molecules),
        )

    @property
    def size(self) -> float:
        cluster_center = self.center
        return max(length(subtract(molecule.center, cluster_center)) + molecule.radius for molecule in self.molecules)

    @staticmethod
    def placement_radius(atom: Atom) -> float:
        if not isfinite(PLACEMENT_RADIUS_FACTOR) or PLACEMENT_RADIUS_FACTOR < 0:
            raise ValueError("PLACEMENT_RADIUS_FACTOR muss endlich und nicht negativ sein")
        if not isfinite(atom.electron_radius) or atom.electron_radius < 0:
            raise ValueError("Darstellungsradien müssen endlich und nicht negativ sein")
        radius = atom.electron_radius * PLACEMENT_RADIUS_FACTOR
        if not isfinite(radius) or radius < 0:
            raise ValueError("Platzierungsradien müssen endlich und nicht negativ sein")
        return radius

    def _overlaps_at(
        self,
        proposed_center: Vec3,
        new_radii: tuple[float, ...],
        rotated_atom_offsets: tuple[Vec3, ...],
        clearance: float,
    ) -> bool:
        for offset, new_radius in zip(rotated_atom_offsets, new_radii):
            proposed_position = add(proposed_center, offset)
            for existing_atom in self.atoms:
                minimum_distance = new_radius + self.placement_radius(existing_atom) + clearance
                if length(subtract(proposed_position, existing_atom.position)) < minimum_distance:
                    return True
        return False

    def add(self, molecule: Molecule, *, clearance: float = 0.0) -> None:
        """Fügt vor der Visualisierung ein kollisionsfrei platziertes Molekül hinzu."""
        if not isinstance(molecule, Molecule):
            raise TypeError("Es kann nur ein Molekül zum Cluster hinzugefügt werden")
        if not molecule.atoms:
            raise ValueError("Ein Molekül zur Platzierung benötigt mindestens ein Atom")
        molecule_electrons = (
            *(electron for bond in molecule.bonds for electron in bond.electrons),
            *(electron for pair in molecule.lone_pairs for electron in pair.electrons),
            *molecule.free_electrons,
        )
        existing_visuals = any(atom.visual is not None for atom in self.atoms) or any(
            electron.visual is not None for electron in self.electrons
        )
        incoming_visuals = any(atom.visual is not None for atom in molecule.atoms) or any(
            electron.visual is not None for electron in molecule_electrons
        )
        if existing_visuals or incoming_visuals:
            raise RuntimeError("Moleküle können nur vor dem Start der Visualisierung ergänzt werden")
        if not isfinite(clearance) or clearance < 0:
            raise ValueError("Der Platzierungsabstand muss endlich und nicht negativ sein")

        updated_molecules = (*self.molecules, molecule)
        self._validate_unique_ids(updated_molecules)
        molecule_center = molecule.center
        pitch = self._rng.uniform(-pi / 2.0, pi / 2.0)
        yaw = self._rng.uniform(0.0, 2.0 * pi)
        rotated_atom_offsets = tuple(
            _rotate_vector(subtract(atom.position, molecule_center), pitch, yaw) for atom in molecule.atoms
        )
        new_radii = tuple(self.placement_radius(atom) for atom in molecule.atoms)
        existing_radii = tuple(self.placement_radius(atom) for atom in self.atoms)
        maximum_pair_radius = max(new_radii, default=0.0) + max(existing_radii, default=0.0)
        direction_z = self._rng.uniform(-1.0, 1.0)
        direction_angle = self._rng.uniform(0.0, 2.0 * pi)
        direction_xy = sqrt(max(0.0, 1.0 - direction_z * direction_z))
        direction = (
            direction_xy * cos(direction_angle),
            direction_xy * sin(direction_angle),
            direction_z,
        )

        cluster_center = self.center
        distance = self.size + molecule.radius + maximum_pair_radius + clearance + 1e-6
        proposed_center = add(cluster_center, scale(direction, distance))
        for _ in range(32):
            if not self._overlaps_at(
                proposed_center,
                new_radii,
                rotated_atom_offsets,
                clearance,
            ):
                break
            distance = distance * 2.0 + maximum_pair_radius + clearance + 1e-6
            proposed_center = add(cluster_center, scale(direction, distance))
        else:
            raise ValueError("Für das Molekül wurde keine kollisionsfreie Startposition gefunden")

        step = max(maximum_pair_radius / 4.0, distance / PLACEMENT_SEARCH_LIMIT, 1e-3)
        safe_distance = distance
        colliding_distance: float | None = None
        for _ in range(PLACEMENT_SEARCH_LIMIT):
            candidate_distance = max(0.0, safe_distance - step)
            candidate_center = add(cluster_center, scale(direction, candidate_distance))
            if self._overlaps_at(
                candidate_center,
                new_radii,
                rotated_atom_offsets,
                clearance,
            ):
                colliding_distance = candidate_distance
                break
            safe_distance = candidate_distance
            if candidate_distance == 0:
                break

        if colliding_distance is not None:
            for _ in range(48):
                if safe_distance - colliding_distance <= PLACEMENT_SEARCH_TOLERANCE:
                    break
                middle_distance = (safe_distance + colliding_distance) / 2.0
                middle_center = add(cluster_center, scale(direction, middle_distance))
                if self._overlaps_at(
                    middle_center,
                    new_radii,
                    rotated_atom_offsets,
                    clearance,
                ):
                    colliding_distance = middle_distance
                else:
                    safe_distance = middle_distance
            proposed_center = add(cluster_center, scale(direction, safe_distance))

        for atom, offset in zip(molecule.atoms, rotated_atom_offsets):
            atom.position = add(proposed_center, offset)
            atom.velocity = _rotate_vector(atom.velocity, pitch, yaw)
        for electron in molecule_electrons:
            electron.position = add(
                proposed_center,
                _rotate_vector(subtract(electron.position, molecule_center), pitch, yaw),
            )
            electron.velocity = _rotate_vector(electron.velocity, pitch, yaw)
            electron.direction = _rotate_vector(electron.direction, pitch, yaw)
            electron.orbit_normal = _rotate_vector(electron.orbit_normal, pitch, yaw)

        self.molecules = updated_molecules
        self._refresh_particles()

    @staticmethod
    def brownian_velocity_kick(dt: float, velocity_factor: float = 1.0) -> Vec3:
        if BROWNIAN_VELOCITY_AMPLITUDE < 0:
            raise ValueError("Die Brownsche Bewegungsamplitude darf nicht negativ sein")
        if DAMPING_REFERENCE_DT <= 0:
            raise ValueError("Das Dämpfungs-Referenzintervall muss positiv sein")
        if dt < 0:
            raise ValueError("Der Simulationszeitschritt darf nicht negativ sein")
        if velocity_factor < 0:
            raise ValueError("Der Brownsche Geschwindigkeitsfaktor darf nicht negativ sein")

        scale_value = BROWNIAN_VELOCITY_AMPLITUDE * velocity_factor * sqrt(dt / DAMPING_REFERENCE_DT)
        return gauss(0.0, scale_value), gauss(0.0, scale_value), gauss(0.0, scale_value)

    def recenter(self, *, remove_com_velocity: bool) -> None:
        total_mass = sum(atom.core_mass for atom in self.atoms) + sum(electron.mass for electron in self.electrons)
        if total_mass <= 0:
            raise ValueError("Die Gesamtmasse des Molekülsystems muss positiv sein")

        center_of_mass = (0.0, 0.0, 0.0)
        center_velocity = (0.0, 0.0, 0.0)
        for atom in self.atoms:
            mass_fraction = atom.core_mass / total_mass
            center_of_mass = add(center_of_mass, scale(atom.position, mass_fraction))
            center_velocity = add(center_velocity, scale(atom.velocity, mass_fraction))
        for electron in self.electrons:
            mass_fraction = electron.mass / total_mass
            center_of_mass = add(center_of_mass, scale(electron.position, mass_fraction))
            center_velocity = add(center_velocity, scale(electron.velocity, mass_fraction))

        for atom in self.atoms:
            atom.position = subtract(atom.position, center_of_mass)
            if remove_com_velocity:
                atom.velocity = subtract(atom.velocity, center_velocity)
        for electron in self.electrons:
            electron.position = subtract(electron.position, center_of_mass)
            if remove_com_velocity:
                electron.velocity = subtract(electron.velocity, center_velocity)

    def integrate(
        self,
        dt: float,
        *,
        inter_molecular_bond_forces_enabled: bool = True,
    ) -> None:
        lone_pair_electron_ids = {
            electron.electron_id for item in self.molecules for pair in item.lone_pairs for electron in pair.electrons
        }
        atoms_by_id = {atom.atom_id: atom for atom in self.atoms}
        reactions_by_atom = {atom.atom_id: (0.0, 0.0, 0.0) for atom in self.atoms}
        current_totals = self.force_model.calculate(
            self.particle_states,
            inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
        )
        current_forces = current_totals.particle_forces

        for electron in self.electrons:
            atom = atoms_by_id[electron.atom_id]
            if electron.bond_id is not None:
                reactions_by_atom[atom.atom_id] = add(
                    reactions_by_atom[atom.atom_id],
                    current_forces[electron.electron_id],
                )

        old_atom_positions = {atom.atom_id: atom.position for atom in self.atoms}
        old_atom_velocities = {atom.atom_id: atom.velocity for atom in self.atoms}
        for atom in self.atoms:
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
            atom.velocity = add(atom.velocity, self.brownian_velocity_kick(dt))
            atom.position = add(atom.position, scale(atom.velocity, dt))

        for electron in self.electrons:
            atom = atoms_by_id[electron.atom_id]
            atom_displacement = subtract(atom.position, old_atom_positions[atom.atom_id])
            atom_velocity_change = subtract(atom.velocity, old_atom_velocities[atom.atom_id])
            electron.position = add(electron.position, atom_displacement)
            electron.velocity = add(electron.velocity, atom_velocity_change)

        updated_forces = self.force_model.calculate(
            self.particle_states,
            inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
        ).particle_forces
        for electron in self.electrons:
            atom = atoms_by_id[electron.atom_id]
            velocity_factor = (
                LONE_PAIR_BROWNIAN_VELOCITY_FACTOR if electron.electron_id in lone_pair_electron_ids else 1.0
            )
            electron.velocity = add(
                electron.velocity,
                self.brownian_velocity_kick(dt, velocity_factor),
            )
            electron.advance(atom, updated_forces[electron.electron_id], dt)

        next_forces = self.force_model.calculate(
            self.particle_states,
            inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
        ).particle_forces
        for electron in self.electrons:
            atom = atoms_by_id[electron.atom_id]
            electron.complete_step(atom, next_forces[electron.electron_id], dt)
            electron.velocity = scale(
                electron.velocity,
                DAMPING_FACTOR ** (dt / DAMPING_REFERENCE_DT),
            )

        self.recenter(remove_com_velocity=True)
