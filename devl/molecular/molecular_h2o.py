"""Startskript für die H2O-Demonstration aus dem Paket ``dirabo.chem``."""

from dirabo.chem.render import show_molecules
from dirabo.chem.water import create_water_molecule

from dirabo.chem.particles import MoleculeCluster

if __name__ == "__main__":
    h2o_molecule_1 = create_water_molecule("H2O-1", (-1.6, 0.0, 0.0))
    h2o_molecule_2 = create_water_molecule("H2O-2", (1.6, 0.0, 0.0))

    lake = MoleculeCluster((h2o_molecule_1, h2o_molecule_2))

    show_molecules(lake)
