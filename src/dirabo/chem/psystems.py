"""Elementdaten und generische Atom-Erzeugung für das Periodensystem."""

from types import MappingProxyType
from dataclasses import dataclass

from dirabo.chem.math import Vec3, direction_at_angle, length
from dirabo.chem.particles import Atom, Electron, ElectronPair


@dataclass(frozen=True)
class ElementProperties:
    symbol: str
    nuclear_charge: float
    core_mass: float
    electron_radius: float
    display_radius: float
    color: Vec3 = (1.0, 1.0, 1.0)
    lone_pair_directions: tuple[tuple[Vec3, Vec3], ...] = ()
    default_lone_pair_count: int = 0


@dataclass(frozen=True)
class ElementParticles:
    atom: Atom
    lone_pairs: tuple[ElectronPair, ...] = ()
    free_electrons: tuple[Electron, ...] = ()


PERIODIC_SYSTEM = MappingProxyType(
    {
        "H": ElementProperties(
            symbol="H",
            nuclear_charge=1.0,
            core_mass=1.008,
            electron_radius=0.30,
            display_radius=0.20,
            color=(1.0, 1.0, 1.0),
        ),
        "O": ElementProperties(
            symbol="O",
            nuclear_charge=6.0,
            core_mass=15.999,
            electron_radius=0.48,
            display_radius=0.32,
            color=(1.0, 0.0, 0.0),
            lone_pair_directions=(
                (direction_at_angle(150, -0.38), direction_at_angle(150, 0.38)),
                (direction_at_angle(255, -0.38), direction_at_angle(255, 0.38)),
                (direction_at_angle(30, -0.38), direction_at_angle(30, 0.38)),
            ),
            default_lone_pair_count=2,
        ),
        "Na": ElementProperties(
            symbol="Na",
            nuclear_charge=1.0,
            core_mass=22.990,
            electron_radius=1.66,
            display_radius=0.42,
            color=(1.0, 0.6, 0.0),
        ),
        "Cl": ElementProperties(
            symbol="Cl",
            nuclear_charge=7.0,
            core_mass=35.45,
            electron_radius=1.02,
            display_radius=0.36,
            color=(0.0, 1.0, 0.0),
            lone_pair_directions=tuple(
                (direction_at_angle(angle, -0.38), direction_at_angle(angle, 0.38)) for angle in (45, 135, 225, 315)
            ),
            default_lone_pair_count=4,
        ),
        "N": ElementProperties(
            symbol="N",
            nuclear_charge=5.0,
            core_mass=14.007,
            electron_radius=0.71,
            display_radius=0.30,
            color=(1.0, 1.0, 0.0),
            lone_pair_directions=((direction_at_angle(90, -0.38), direction_at_angle(90, 0.38)),),
            default_lone_pair_count=1,
        ),
        "C": ElementProperties(
            symbol="C",
            nuclear_charge=6.0,
            core_mass=12.011,
            electron_radius=0.76,
            display_radius=0.30,
            color=(0.3, 0.3, 0.3),
        ),
    }
)


def get_color(element: str) -> Vec3:
    """Gibt die RGB-Farbe eines Elements zurück, unbekannte Elemente sind weiß."""
    properties = PERIODIC_SYSTEM.get(element)
    return (1.0, 1.0, 1.0) if properties is None else properties.color


def create_atom(
    element: str,
    atom_id: str,
    position: Vec3,
    *,
    nuclear_charge_override: float | None = None,
) -> Atom:
    """Erzeugt ein Atom anhand der Elementdaten des Periodensystems."""
    try:
        properties = PERIODIC_SYSTEM[element]
    except KeyError:
        raise ValueError(f"Element {element!r} ist nicht im Periodensystem enthalten") from None

    return Atom(
        atom_id=atom_id,
        element=properties.symbol,
        position=position,
        nuclear_charge=(properties.nuclear_charge if nuclear_charge_override is None else nuclear_charge_override),
        core_mass=properties.core_mass,
        electron_radius=properties.electron_radius,
        display_radius=properties.display_radius,
    )


def create_element(
    element: str,
    atom_id: str,
    position: Vec3,
    *,
    lone_pair_count: int | None = None,
    free_electron_directions: tuple[Vec3, ...] = (),
    nuclear_charge_override: float | None = None,
) -> ElementParticles:
    """Erzeugt ein Atom samt optionalen freien Elektronen und Elektronenpaaren.

    Wenn ``lone_pair_count`` nicht angegeben ist, wird die konfigurierte
    Standardzahl des Elements verwendet. Freie Elektronen werden jeweils durch
    eine Richtung relativ zum Atomkern beschrieben.
    """
    try:
        properties = PERIODIC_SYSTEM[element]
    except KeyError:
        raise ValueError(f"Element {element!r} ist nicht im Periodensystem enthalten") from None

    if not atom_id:
        raise ValueError("Ein Atom benötigt eine ID")

    requested_pair_count = properties.default_lone_pair_count if lone_pair_count is None else lone_pair_count
    if requested_pair_count < 0:
        raise ValueError("Die Anzahl der Elektronenpaare darf nicht negativ sein")
    if requested_pair_count > len(properties.lone_pair_directions):
        raise ValueError(
            f"Für Element {element!r} sind nur " f"{len(properties.lone_pair_directions)} Elektronenpaare konfiguriert"
        )
    if any(
        length(direction) == 0 or length((direction[0], direction[1], 0.0)) == 0
        for direction in free_electron_directions
    ):
        raise ValueError(
            "Freie Elektronen benötigen eine Bahnrichtung, die nicht parallel " "zur Standard-Bahnnormale ist"
        )

    atom = create_atom(
        element,
        atom_id,
        position,
        nuclear_charge_override=nuclear_charge_override,
    )
    lone_pairs = tuple(
        ElectronPair(
            f"{atom_id}:lone-pair-{pair_index + 1}",
            (
                Electron(
                    f"{atom_id}:lone-pair-{pair_index + 1}:e1",
                    atom_id,
                    directions[0],
                ),
                Electron(
                    f"{atom_id}:lone-pair-{pair_index + 1}:e2",
                    atom_id,
                    directions[1],
                    orbit_direction=-1.0,
                ),
            ),
        )
        for pair_index, directions in enumerate(properties.lone_pair_directions[:requested_pair_count])
    )
    free_electrons = tuple(
        Electron(
            f"{atom_id}:free-electron-{electron_index + 1}",
            atom_id,
            direction,
        )
        for electron_index, direction in enumerate(free_electron_directions)
    )
    return ElementParticles(atom, lone_pairs, free_electrons)


def hydrogen(
    atom_id: str,
    position: Vec3,
    *,
    lone_pair_count: int | None = None,
    free_electron_directions: tuple[Vec3, ...] = (),
    nuclear_charge_override: float | None = None,
) -> ElementParticles:
    """Erzeugt Wasserstoff samt optionalen freien Elektronen."""
    return create_element(
        "H",
        atom_id,
        position,
        lone_pair_count=lone_pair_count,
        free_electron_directions=free_electron_directions,
        nuclear_charge_override=nuclear_charge_override,
    )


def oxygen(
    atom_id: str,
    position: Vec3,
    *,
    lone_pair_count: int | None = None,
    free_electron_directions: tuple[Vec3, ...] = (),
    nuclear_charge_override: float | None = None,
) -> ElementParticles:
    """Erzeugt Sauerstoff samt standardmäßig zwei freien Elektronenpaaren."""
    return create_element(
        "O",
        atom_id,
        position,
        lone_pair_count=lone_pair_count,
        free_electron_directions=free_electron_directions,
        nuclear_charge_override=nuclear_charge_override,
    )


def sodium(
    atom_id: str,
    position: Vec3,
    *,
    lone_pair_count: int | None = None,
    free_electron_directions: tuple[Vec3, ...] = (),
    nuclear_charge_override: float | None = None,
) -> ElementParticles:
    """Erzeugt Natrium anhand der Elementdaten."""
    return create_element(
        "Na",
        atom_id,
        position,
        lone_pair_count=lone_pair_count,
        free_electron_directions=free_electron_directions,
        nuclear_charge_override=nuclear_charge_override,
    )


def chlorine(
    atom_id: str,
    position: Vec3,
    *,
    lone_pair_count: int | None = None,
    free_electron_directions: tuple[Vec3, ...] = (),
    nuclear_charge_override: float | None = None,
) -> ElementParticles:
    """Erzeugt Chlor anhand der Elementdaten."""
    return create_element(
        "Cl",
        atom_id,
        position,
        lone_pair_count=lone_pair_count,
        free_electron_directions=free_electron_directions,
        nuclear_charge_override=nuclear_charge_override,
    )


def nitrogen(
    atom_id: str,
    position: Vec3,
    *,
    lone_pair_count: int | None = None,
    free_electron_directions: tuple[Vec3, ...] = (),
) -> ElementParticles:
    """Erzeugt Stickstoff anhand der Elementdaten."""
    return create_element(
        "N",
        atom_id,
        position,
        lone_pair_count=lone_pair_count,
        free_electron_directions=free_electron_directions,
    )


def carbon(
    atom_id: str,
    position: Vec3,
    *,
    lone_pair_count: int | None = None,
    free_electron_directions: tuple[Vec3, ...] = (),
) -> ElementParticles:
    """Erzeugt Kohlenstoff anhand der Elementdaten."""
    return create_element(
        "C",
        atom_id,
        position,
        lone_pair_count=lone_pair_count,
        free_electron_directions=free_electron_directions,
    )
