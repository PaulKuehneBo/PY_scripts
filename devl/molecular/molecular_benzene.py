"""Startskript für ein vereinfachtes Kekulé-Modell von Benzol (C6H6)."""

from dirabo.chem.aromatic_functional_groups import create_aromatic_ring
from dirabo.chem.math import Vec3
from dirabo.chem.particles import Molecule, MoleculeCluster
from dirabo.chem.render import show_molecules


def create_benzol_molecule(
    instance_id: str = "benzol-1",
    position: Vec3 = (0.0, 0.0, 0.0),
) -> Molecule:
    """Erzeugt Benzol über den generischen aromatischen Ring-Erzeuger."""
    ring = create_aromatic_ring(instance_id, ("C",) * 6, position=position)
    return ring.to_molecule(name=f"Benzol {instance_id}")


def create_benzol_cluster(
    instance_id: str = "benzol-1",
    position: Vec3 = (0.0, 0.0, 0.0),
) -> MoleculeCluster:
    """Erzeugt den Molekülcluster für die Benzol-Demonstration."""
    return MoleculeCluster(create_benzol_molecule(instance_id, position))


if __name__ == "__main__":
    show_molecules(create_benzol_cluster())
