"""Startskript für ein Ethanolmodell (C2H6O) in der 3D-Simulation."""

from dirabo.chem.organic_functional_groups import FunctionalGroupType, create_carbon_chain
from dirabo.chem.particles import Molecule, MoleculeCluster
from dirabo.chem.render import show_molecules


def create_ethanol_molecule(
    instance_id: str = "ethanol-1",
    position: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> Molecule:
    """Erzeugt ein Ethanolmolekül mit acht Bindungen und zwei O-Lone-Pairs."""
    return create_carbon_chain(
        instance_id,
        (3, 2),
        substituents={1: FunctionalGroupType.HYDROXYL},
        position=position,
    ).to_molecule(name=f"Ethanol {instance_id}")


def create_ethanol_cluster(
    instance_id: str = "ethanol-1",
    position: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> MoleculeCluster:
    """Erzeugt den Molekülcluster für die Ethanol-Demonstration."""
    return MoleculeCluster(create_ethanol_molecule(instance_id, position))


if __name__ == "__main__":
    show_molecules(create_ethanol_cluster())
