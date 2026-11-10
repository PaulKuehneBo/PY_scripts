from math import cos, radians, sin

from dirabo.chem.math import Vec3, add, direction_towards
from dirabo.chem.particles import Bond, Electron, Molecule
from dirabo.chem.psystems import hydrogen, oxygen

# 1 Welt-Einheit entspricht 1 Å. Die Kugelradien sind Darstellungswerte.
OH_DISTANCE = 0.96
HOH_ANGLE = 104.5


def create_water_molecule(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
    *,
    oxygen_lone_pair_count: int | None = None,
    oxygen_free_electron_directions: tuple[Vec3, ...] = (),
) -> Molecule:
    """Erzeugt eine unabhängige H2O-Instanz mit eindeutigen IDs."""
    if not instance_id:
        raise ValueError("Eine Molekülinstanz benötigt eine ID")

    oxygen_id = f"{instance_id}:O1"
    hydrogen1_id = f"{instance_id}:H1"
    hydrogen2_id = f"{instance_id}:H2"
    bond1_id = f"{instance_id}:O1-H1"
    bond2_id = f"{instance_id}:O1-H2"
    oxygen_position = position
    hydrogen1_position = add(position, (OH_DISTANCE, 0.0, 0.0))
    angle = radians(HOH_ANGLE)
    hydrogen2_position = add(
        position,
        (OH_DISTANCE * cos(angle), OH_DISTANCE * sin(angle), 0.0),
    )

    oxygen_particles = oxygen(
        oxygen_id,
        oxygen_position,
        lone_pair_count=oxygen_lone_pair_count,
        free_electron_directions=oxygen_free_electron_directions,
    )
    hydrogen1_particles = hydrogen(hydrogen1_id, hydrogen1_position)
    hydrogen2_particles = hydrogen(hydrogen2_id, hydrogen2_position)

    bond1 = Bond(
        bond1_id,
        (
            Electron(
                f"{instance_id}:eO1",
                oxygen_id,
                direction_towards(oxygen_position, hydrogen1_position, -0.16),
                bond_id=bond1_id,
                orbit_direction=1.0,
            ),
            Electron(
                f"{instance_id}:eH1",
                hydrogen1_id,
                direction_towards(hydrogen1_position, oxygen_position, 0.16),
                bond_id=bond1_id,
                orbit_direction=-1.0,
            ),
        ),
    )
    bond2 = Bond(
        bond2_id,
        (
            Electron(
                f"{instance_id}:eO2",
                oxygen_id,
                direction_towards(oxygen_position, hydrogen2_position, -0.16),
                bond_id=bond2_id,
                orbit_direction=-1.0,
            ),
            Electron(
                f"{instance_id}:eH2",
                hydrogen2_id,
                direction_towards(hydrogen2_position, oxygen_position, 0.16),
                bond_id=bond2_id,
                orbit_direction=1.0,
            ),
        ),
    )
    return Molecule(
        f"H2O {instance_id}",
        atoms=(oxygen_particles.atom, hydrogen1_particles.atom, hydrogen2_particles.atom),
        bonds=(bond1, bond2),
        lone_pairs=oxygen_particles.lone_pairs,
        free_electrons=oxygen_particles.free_electrons,
    )


def create_OH_group(
    instance_id: str,
    position: Vec3 = (0.0, 0.0, 0.0),
) -> Molecule:
    """Erzeugt ein Hydroxid-Ion OH- mit drei freien Elektronenpaaren."""
    if not instance_id:
        raise ValueError("Eine OH-Gruppe benötigt eine ID")

    oxygen_id = f"{instance_id}:O1"
    hydrogen_id = f"{instance_id}:H1"
    bond_id = f"{instance_id}:O1-H1"
    oxygen_position = position
    hydrogen_position = add(position, (OH_DISTANCE, 0.0, 0.0))

    oxygen_particles = oxygen(
        oxygen_id,
        oxygen_position,
        lone_pair_count=3,
        nuclear_charge_override=6.0,
    )
    hydrogen_particles = hydrogen(
        hydrogen_id,
        hydrogen_position,
        nuclear_charge_override=1.0,
    )
    bond = Bond(
        bond_id,
        (
            Electron(
                f"{instance_id}:eO1",
                oxygen_id,
                direction_towards(oxygen_position, hydrogen_position, -0.16),
                bond_id=bond_id,
                orbit_direction=1.0,
            ),
            Electron(
                f"{instance_id}:eH1",
                hydrogen_id,
                direction_towards(hydrogen_position, oxygen_position, 0.16),
                bond_id=bond_id,
                orbit_direction=-1.0,
            ),
        ),
    )
    return Molecule(
        f"OH- {instance_id}",
        atoms=(oxygen_particles.atom, hydrogen_particles.atom),
        bonds=(bond,),
        lone_pairs=oxygen_particles.lone_pairs,
    )
