"""Experiment: Wechselwirkung eines Protons mit einer Hydroxidgruppe."""

from dirabo.chem import particles as simulation_parameters
from dirabo.chem.particles import Molecule, MoleculeCluster
from dirabo.chem.psystems import hydrogen
from dirabo.chem.render import show_molecules
from dirabo.chem.water import create_OH_group

EXPERIMENT_BROWNIAN_VELOCITY_AMPLITUDE = 0.0015


def create_proton(
    instance_id: str,
    position: tuple[float, float, float],
) -> Molecule:
    """Erzeugt ein Proton ohne Elektronen."""
    proton = hydrogen(
        f"{instance_id}:H+",
        position,
        nuclear_charge_override=1.0,
    )
    return Molecule(
        name=f"H+ {instance_id}",
        atoms=(proton.atom,),
        bonds=(),
        lone_pairs=(),
    )


def create_oh_proton_system() -> MoleculeCluster:
    """Setzt das Proton an die anfängliche Position des OH-Wasserstoffelektrons."""
    hydroxide = create_OH_group("OH-1", (-1.5, 0.0, 0.0))
    hydroxide_hydrogen = hydroxide.atoms[1]
    hydrogen_bond_electron = next(
        electron for electron in hydroxide.bonds[0].electrons if electron.atom_id == hydroxide_hydrogen.atom_id
    )
    hydrogen_bond_electron.initialize(hydroxide_hydrogen)
    proton = create_proton("H+-1", hydrogen_bond_electron.position)
    return MoleculeCluster((hydroxide, proton))


if __name__ == "__main__":
    # Isoliert die Coulomb-Wechselwirkung, ohne thermische Stöße oder
    # eine künstliche Verstärkung der Kern-Kern-Abstoßung.
    simulation_parameters.BROWNIAN_VELOCITY_AMPLITUDE = EXPERIMENT_BROWNIAN_VELOCITY_AMPLITUDE

    show_molecules(create_oh_proton_system())
