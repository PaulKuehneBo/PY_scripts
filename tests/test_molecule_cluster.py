import unittest
from random import Random
from types import SimpleNamespace
from unittest.mock import patch

from dirabo.chem.math import length, subtract
from dirabo.chem.keybd_handler import KeyboardControls
from dirabo.chem.molecule_cluster import PLACEMENT_RADIUS_FACTOR, MoleculeCluster
from dirabo.chem.water import create_water_molecule
from devl.molecular.molecular_h2o import _parse_args


class MoleculeClusterTests(unittest.TestCase):
    def test_m_toggles_molecule_bounds(self) -> None:
        bound = SimpleNamespace(visible=False)
        controls = KeyboardControls(
            MoleculeCluster(create_water_molecule("H2O-1")),
            {},
            (bound,),
        )

        controls.handle_keydown(SimpleNamespace(key="m"))
        self.assertTrue(bound.visible)
        self.assertTrue(controls.molecule_bounds_visible)

        controls.handle_keydown(SimpleNamespace(key="M"))
        self.assertFalse(bound.visible)
        self.assertFalse(controls.molecule_bounds_visible)

    def test_stop_at_start_initializes_pause_and_space_resumes(self) -> None:
        self.assertFalse(_parse_args([]).stop_at_start)
        self.assertTrue(_parse_args(["--stop-at-start"]).stop_at_start)

        controls = KeyboardControls(
            MoleculeCluster(create_water_molecule("H2O-1")),
            {},
            start_paused=True,
        )
        self.assertTrue(controls.paused)

        controls.handle_keydown(SimpleNamespace(key=" "))
        self.assertFalse(controls.paused)

    def test_molecule_and_cluster_geometry(self) -> None:
        molecule = create_water_molecule("H2O-1")
        center = tuple(
            sum(atom.position[axis] for atom in molecule.atoms) / len(molecule.atoms)
            for axis in range(3)
        )
        expected_radius = max(length(subtract(atom.position, center)) for atom in molecule.atoms)

        self.assertEqual(molecule.center, center)
        self.assertAlmostEqual(molecule.radius, expected_radius)

        cluster = MoleculeCluster(molecule)
        self.assertEqual(cluster.center, molecule.center)
        self.assertAlmostEqual(cluster.size, molecule.radius)

    def test_add_places_molecule_without_sphere_overlaps(self) -> None:
        first = create_water_molecule("H2O-1")
        second = create_water_molecule("H2O-2")
        cluster = MoleculeCluster(first, rng=Random(7))

        cluster.add(second)

        self.assertEqual(len(cluster.molecules), 2)
        self.assertEqual(len(cluster.atoms), 6)
        self.assertEqual(len(cluster.electrons), 16)
        self.assertEqual(len(cluster.particle_states), 22)
        for first_atom in first.atoms:
            for second_atom in second.atoms:
                minimum_distance = PLACEMENT_RADIUS_FACTOR * (
                    first_atom.display_radius + second_atom.display_radius
                )
                actual_distance = length(subtract(first_atom.position, second_atom.position))
                self.assertGreaterEqual(actual_distance + 1e-8, minimum_distance)

    def test_duplicate_ids_leave_cluster_unchanged(self) -> None:
        first = create_water_molecule("H2O-1")
        cluster = MoleculeCluster(first, rng=Random(7))
        original_molecules = cluster.molecules

        with self.assertRaises(ValueError):
            cluster.add(create_water_molecule("H2O-1"))

        self.assertEqual(cluster.molecules, original_molecules)
        self.assertEqual(len(cluster.atoms), 3)

    def test_add_rejects_changes_after_visualization_started(self) -> None:
        first = create_water_molecule("H2O-1")
        cluster = MoleculeCluster(first, rng=Random(7))
        first.atoms[0].visual = object()

        with self.assertRaises(RuntimeError):
            cluster.add(create_water_molecule("H2O-2"))

        self.assertEqual(len(cluster.molecules), 1)

    def test_add_respects_clearance(self) -> None:
        first = create_water_molecule("H2O-1")
        second = create_water_molecule("H2O-2")
        cluster = MoleculeCluster(first, rng=Random(7))
        clearance = 0.15

        cluster.add(second, clearance=clearance)

        for first_atom in first.atoms:
            for second_atom in second.atoms:
                minimum_distance = PLACEMENT_RADIUS_FACTOR * (
                    first_atom.display_radius + second_atom.display_radius
                ) + clearance
                actual_distance = length(subtract(first_atom.position, second_atom.position))
                self.assertGreaterEqual(actual_distance + 1e-8, minimum_distance)

    def test_placement_radius_factor_scales_display_radii(self) -> None:
        first = create_water_molecule("H2O-1")
        second = create_water_molecule("H2O-2")
        factor = 1.5

        with patch("dirabo.chem.molecule_cluster.PLACEMENT_RADIUS_FACTOR", factor):
            cluster = MoleculeCluster(first, rng=Random(7))
            cluster.add(second)

        for first_atom in first.atoms:
            for second_atom in second.atoms:
                minimum_distance = factor * (
                    first_atom.display_radius + second_atom.display_radius
                )
                actual_distance = length(subtract(first_atom.position, second_atom.position))
                self.assertGreaterEqual(actual_distance + 1e-8, minimum_distance)

    def test_recenter_and_integrate_are_cluster_methods(self) -> None:
        molecule = create_water_molecule("H2O-1")
        cluster = MoleculeCluster(molecule, rng=Random(7))
        atoms_by_id = {atom.atom_id: atom for atom in cluster.atoms}
        for electron in cluster.electrons:
            electron.initialize(atoms_by_id[electron.atom_id])

        cluster.integrate(0.0)
        cluster.recenter(remove_com_velocity=False)


if __name__ == "__main__":
    unittest.main()
