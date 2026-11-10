from __future__ import annotations

from typing import Any

from collections import deque
from dataclasses import dataclass, field

from random import gauss
from math import sqrt
from .math import Vec3, add, subtract, scale, dot, cross, length, remove_radial_velocity

from vpython import color

ELECTRON_DISPLAY_RADIUS = 0.055
NUCLEUS_REPULSION_FACTOR = 1.0
BOND_MINIMUM_LENGTH = 0.22


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
    orbit_reaction: Vec3 = field(default=(0.0, 0.0, 0.0), init=False)
    applied_force: Vec3 = field(default=(0.0, 0.0, 0.0), init=False)
    velocity: Vec3 = field(default=(0.0, 0.0, 0.0), init=False)
    visual: Any = field(default=None, init=False, repr=False, compare=False)
    force_visual: Any = field(default=None, init=False, repr=False, compare=False)
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
    force_visual: Any = field(default=None, init=False, repr=False, compare=False)
    individual_force_visuals: dict[str, Any] = field(default_factory=dict, init=False, repr=False, compare=False)

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


@dataclass(frozen=True)
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
class ForceContribution:
    source_id: str
    source_atom_id: str
    target_id: str
    target_atom_id: str
    interaction: str
    force: Vec3
    behavior: str


@dataclass(frozen=True)
class AtomForceDecomposition:
    internal_force: Vec3
    external_force: Vec3


@dataclass(frozen=True)
class ForceReport:
    contributions: tuple[ForceContribution, ...]

    def for_atom(self, atom_id: str) -> AtomForceDecomposition:
        internal = (0.0, 0.0, 0.0)
        external = (0.0, 0.0, 0.0)
        for contribution in self.contributions:
            if contribution.target_atom_id != atom_id:
                continue
            if contribution.source_atom_id == atom_id:
                internal = add(internal, contribution.force)
            else:
                external = add(external, contribution.force)
        return AtomForceDecomposition(internal, external)


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
    cutoff_radius = bond_minimum_length + 2 * target_electron_radius
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
    """Coulomb-Wechselwirkungen in skalierten Visualisierungseinheiten."""

    strength: float = 0.02
    softening: float = 0.12
    pair_target_length: float = 0.22
    pair_length_stiffness: float = 1.0

    def __post_init__(self) -> None:
        if self.strength < 0:
            raise ValueError("Die Coulomb-Stärke darf nicht negativ sein")
        if self.softening <= 0:
            raise ValueError("Der Softening-Abstand muss positiv sein")
        if self.pair_target_length <= 0:
            raise ValueError("Die Ziel-Paarlänge muss positiv sein")
        if self.pair_length_stiffness < 0:
            raise ValueError("Die Steifigkeit der Paarlänge darf nicht negativ sein")

    def calculate(
        self,
        particles: tuple[ParticleState, ...],
        *,
        inter_molecular_bond_forces_enabled: bool = True,
    ) -> ForceReport:
        contributions: list[ForceContribution] = []
        softening_squared = self.softening * self.softening

        for index, first in enumerate(particles):
            for second in particles[index + 1 :]:
                interaction = _interaction_type(first, second)
                is_inter_molecular = False
                if interaction == "electron-nucleus":
                    nucleus, electron = (first, second) if first.particle_type == "nucleus" else (second, first)
                    is_inter_molecular = (
                        first.molecule_id is not None
                        and second.molecule_id is not None
                        and first.molecule_id != second.molecule_id
                    )
                    if not is_inter_molecular:
                        if electron.bond_id is not None or nucleus.atom_id == electron.atom_id:
                            continue
                    elif electron.bond_id is not None and not inter_molecular_bond_forces_enabled:
                        continue

                displacement = subtract(first.position, second.position)
                distance_squared = dot(displacement, displacement) + softening_squared
                effective_charge_product = first.charge * second.charge
                behavior = "attraction" if effective_charge_product < 0 else "repulsion"
                if interaction == "electron-nucleus":
                    nucleus, electron = (first, second) if first.particle_type == "nucleus" else (second, first)
                    force_on_first = _nuclear_electron_coulomb_force(
                        nucleus,
                        electron,
                        self.strength,
                        softening_squared,
                        inter_molecular=is_inter_molecular,
                    )
                    if first.particle_type == "electron":
                        force_on_first = scale(force_on_first, -1.0)
                else:
                    scale_value = self.strength * effective_charge_product / (distance_squared * sqrt(distance_squared))
                    force_on_first = scale(displacement, scale_value)
                contributions.extend(
                    (
                        ForceContribution(
                            second.particle_id,
                            second.atom_id,
                            first.particle_id,
                            first.atom_id,
                            interaction,
                            force_on_first,
                            behavior,
                        ),
                        ForceContribution(
                            first.particle_id,
                            first.atom_id,
                            second.particle_id,
                            second.atom_id,
                            interaction,
                            scale(force_on_first, -1.0),
                            behavior,
                        ),
                    )
                )

        paired_electrons: dict[str, list[ParticleState]] = {}
        for particle in particles:
            if particle.particle_type == "electron" and particle.pair_id is not None:
                paired_electrons.setdefault(particle.pair_id, []).append(particle)

        for pair in paired_electrons.values():
            if len(pair) != 2:
                raise ValueError("Eine Bindung oder ein Elektronenpaar muss genau zwei Elektronen enthalten")
            first, second = pair
            displacement = subtract(first.position, second.position)
            distance = length(displacement)
            if distance == 0:
                continue
            deviation = distance - self.pair_target_length
            force_on_first = scale(
                displacement,
                -self.pair_length_stiffness * deviation / distance,
            )
            behavior = "attraction" if deviation > 0 else "repulsion" if deviation < 0 else "balanced"
            contributions.extend(
                (
                    ForceContribution(
                        second.particle_id,
                        second.atom_id,
                        first.particle_id,
                        first.atom_id,
                        "pair-constraint",
                        force_on_first,
                        behavior,
                    ),
                    ForceContribution(
                        first.particle_id,
                        first.atom_id,
                        second.particle_id,
                        second.atom_id,
                        "pair-constraint",
                        scale(force_on_first, -1.0),
                        behavior,
                    ),
                )
            )

        return ForceReport(tuple(contributions))


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


def _force_color(behavior: str | None) -> Any:
    if behavior == "attraction":
        return color.blue
    if behavior == "repulsion":
        return color.red
    return color.gray(0.6)


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


class MoleculeCluster:
    def __init__(self, molecules: Molecule | tuple[Molecule, ...]) -> None:
        self.molecules = (molecules,) if isinstance(molecules, Molecule) else tuple(molecules)
        if not all(isinstance(molecule, Molecule) for molecule in self.molecules):
            raise ValueError("Mindestens ein Molekül muss simuliert werden")

    @property
    def electrons(self) -> tuple[Electron, ...]:
        return tuple(
            electron
            for item in self.molecules
            for electron in (
                *(electron for bond in item.bonds for electron in bond.electrons),
                *(electron for pair in item.lone_pairs for electron in pair.electrons),
                *item.free_electrons,
            )
        )

    @property
    def atoms(self) -> tuple[Atom, ...]:
        return tuple(atom for molecule in self.molecules for atom in molecule.atoms)

    @property
    def bonds(self) -> tuple[Bond, ...]:
        return tuple(bond for molecule in self.molecules for bond in molecule.bonds)

    @property
    def lone_pairs(self) -> tuple[ElectronPair, ...]:
        return tuple(pair for molecule in self.molecules for pair in molecule.lone_pairs)


def _particle_states(molecules: MoleculeCluster) -> tuple[ParticleState, ...]:
    molecule_ids_by_atom = {atom.atom_id: id(item) for item in molecules.molecules for atom in item.atoms}
    pair_ids_by_electron = {
        electron.electron_id: bond.bond_id
        for item in molecules.molecules
        for bond in item.bonds
        for electron in bond.electrons
    }
    pair_ids_by_electron.update(
        {
            electron.electron_id: pair.pair_id
            for item in molecules.molecules
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
        for item in molecules.molecules
        for atom in item.atoms
    ]
    atoms_by_id = {atom.atom_id: atom for molecule in molecules.molecules for atom in molecule.atoms}
    bond_minimum_lengths_by_electron_id: dict[str, float] = {}
    for item in molecules.molecules:
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
            atoms_by_id[electron.atom_id].position,
            atoms_by_id[electron.atom_id].display_radius,
            None,
            bond_minimum_lengths_by_electron_id.get(
                electron.electron_id,
                BOND_MINIMUM_LENGTH,
            ),
        )
        for electron in molecules.electrons
    )
    particle_ids = [particle.particle_id for particle in particles]
    if len(particle_ids) != len(set(particle_ids)):
        raise ValueError("Partikel-IDs müssen über alle Moleküle hinweg eindeutig sein")
    return tuple(particles)


def force_report(
    molecules: MoleculeCluster,
    *,
    inter_molecular_bond_forces_enabled: bool = True,
) -> ForceReport:
    """Berechnet alle paarweisen Kräfte des gesamten Molekülsystems."""
    return force_model.calculate(
        _particle_states(molecules),
        inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
    )


def _forces_by_particle(report: ForceReport) -> dict[str, Vec3]:
    forces: dict[str, Vec3] = {}
    for contribution in report.contributions:
        forces[contribution.target_id] = add(
            forces.get(contribution.target_id, (0.0, 0.0, 0.0)),
            contribution.force,
        )
    return forces


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
    current_report = force_report(
        molecules,
        inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
    )
    current_forces = _forces_by_particle(current_report)

    for electron in electrons:
        atom = atoms_by_id[electron.atom_id]
        if electron.bond_id is not None:
            reactions_by_atom[atom.atom_id] = add(
                reactions_by_atom[atom.atom_id],
                current_forces[electron.electron_id],
            )

    nucleus_repulsions = {atom.atom_id: (0.0, 0.0, 0.0) for atom in atoms}
    nucleus_electron_forces = {atom.atom_id: (0.0, 0.0, 0.0) for atom in atoms}
    for contribution in current_report.contributions:
        if contribution.target_atom_id not in nucleus_repulsions:
            continue
        if contribution.interaction == "nucleus-nucleus":
            nucleus_repulsions[contribution.target_atom_id] = add(
                nucleus_repulsions[contribution.target_atom_id],
                contribution.force,
            )
        elif contribution.interaction == "electron-nucleus":
            nucleus_electron_forces[contribution.target_atom_id] = add(
                nucleus_electron_forces[contribution.target_atom_id],
                contribution.force,
            )

    old_atom_positions = {atom.atom_id: atom.position for atom in atoms}
    old_atom_velocities = {atom.atom_id: atom.velocity for atom in atoms}
    for atom in atoms:
        atom.orbit_reaction = reactions_by_atom[atom.atom_id]
        atom.applied_force = add(
            add(
                atom.orbit_reaction,
                scale(nucleus_repulsions[atom.atom_id], NUCLEUS_REPULSION_FACTOR),
            ),
            nucleus_electron_forces[atom.atom_id],
        )
        atom.velocity = scale(
            add(atom.velocity, scale(atom.applied_force, dt / atom.core_mass)),
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

    updated_forces = _forces_by_particle(
        force_report(
            molecules,
            inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
        )
    )
    for electron in electrons:
        atom = atoms_by_id[electron.atom_id]
        velocity_factor = LONE_PAIR_BROWNIAN_VELOCITY_FACTOR if electron.electron_id in lone_pair_electron_ids else 1.0
        electron.velocity = add(
            electron.velocity,
            _brownian_velocity_kick(dt, velocity_factor),
        )
        electron.advance(atom, updated_forces[electron.electron_id], dt)

    next_forces = _forces_by_particle(
        force_report(
            molecules,
            inter_molecular_bond_forces_enabled=inter_molecular_bond_forces_enabled,
        )
    )
    for electron in electrons:
        atom = atoms_by_id[electron.atom_id]
        electron.complete_step(atom, next_forces[electron.electron_id], dt)
        electron.velocity = scale(
            electron.velocity,
            DAMPING_FACTOR ** (dt / DAMPING_REFERENCE_DT),
        )

    recenter_system(molecules, remove_com_velocity=True)
