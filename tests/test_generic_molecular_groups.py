import json
import unittest
from pathlib import Path

from dirabo.chem.generic_molecular_groups import GenericMolecularGroups
from devl.molecular.molecular_soup import _parse_args, create_cluster

CONFIG_PATH = Path(__file__).parents[1] / "devl" / "molecular" / "molecular_generic.json"


class GenericMolecularGroupsTests(unittest.TestCase):
    def test_ethanol_is_assembled_from_json_group_components(self) -> None:
        library = GenericMolecularGroups.load(CONFIG_PATH)
        molecule = library.create_molecule("ethanol", instance_id="ethanol-test")

        self.assertEqual(
            library.molecule_names,
            ("ethanol", "dimethyl_ether", "acetone"),
        )
        self.assertEqual(len(molecule.atoms), 9)
        self.assertEqual(len(molecule.bonds), 8)
        self.assertEqual(len(molecule.lone_pairs), 2)
        electron_count = (
            sum(len(bond.electrons) for bond in molecule.bonds)
            + sum(len(pair.electrons) for pair in molecule.lone_pairs)
            + len(molecule.free_electrons)
        )
        self.assertEqual(electron_count, 20)
        self.assertEqual(len({atom.atom_id for atom in molecule.atoms}), 9)
        self.assertEqual(len({bond.bond_id for bond in molecule.bonds}), 8)

    def test_dimethyl_ether_and_acetone_are_composed_from_groups(self) -> None:
        library = GenericMolecularGroups.load(CONFIG_PATH)
        ether = library.create_molecule("dimethyl_ether", "ether-test")
        acetone = library.create_molecule("acetone", "acetone-test")

        self.assertEqual(len(ether.atoms), 9)
        self.assertEqual(len(ether.bonds), 8)
        self.assertEqual(len(ether.lone_pairs), 2)

        self.assertEqual(len(acetone.atoms), 10)
        self.assertEqual(len(acetone.bonds), 10)
        self.assertEqual(len(acetone.lone_pairs), 2)
        carbonyl_bonds = [bond for bond in acetone.bonds if ":carbonyl:C=O:" in bond.bond_id]
        self.assertEqual(len(carbonyl_bonds), 2)
        self.assertEqual(
            {electron.atom_id for bond in carbonyl_bonds for electron in bond.electrons},
            {"acetone-test:carbonyl:C", "acetone-test:carbonyl:O"},
        )

    def test_cli_molecule_count_options_come_from_json_names(self) -> None:
        library = GenericMolecularGroups.load(CONFIG_PATH)
        options = _parse_args(
            ["--ethanol", "2", "--acetone", "1"],
            library.molecule_names,
        )

        self.assertEqual(
            options.molecule_counts,
            {"ethanol": 2, "dimethyl_ether": 0, "acetone": 1},
        )

    def test_selected_molecule_counts_build_a_mixed_cluster(self) -> None:
        library = GenericMolecularGroups.load(CONFIG_PATH)
        cluster = create_cluster(
            library,
            {"ethanol": 2, "dimethyl_ether": 1, "acetone": 1},
        )

        self.assertEqual(
            [molecule.name.split()[0] for molecule in cluster.molecules],
            ["Ethanol", "Ethanol", "Dimethylether", "Aceton"],
        )
        self.assertEqual(len(cluster.molecules), 4)
        self.assertEqual(len(cluster.atoms), 37)
        self.assertEqual(len({atom.atom_id for atom in cluster.atoms}), 37)
        self.assertEqual(len({electron.electron_id for electron in cluster.electrons}), 84)

    def test_cli_rejects_all_zero_molecule_counts(self) -> None:
        library = GenericMolecularGroups.load(CONFIG_PATH)

        with self.assertRaises(SystemExit):
            _parse_args(["--ethanol", "0"], library.molecule_names)

    def test_instances_receive_independent_ids(self) -> None:
        library = GenericMolecularGroups.load(CONFIG_PATH)
        first = library.create_molecule("ethanol", instance_id="first")
        second = library.create_molecule("ethanol", instance_id="second")

        self.assertTrue({atom.atom_id for atom in first.atoms}.isdisjoint(atom.atom_id for atom in second.atoms))
        self.assertTrue({bond.bond_id for bond in first.bonds}.isdisjoint(bond.bond_id for bond in second.bonds))

    def test_unknown_connection_port_is_rejected(self) -> None:
        definition = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        definition["molecules"]["ethanol"]["connections"][0]["to"]["port"] = "missing"

        with self.assertRaisesRegex(ValueError, "keinen Port 'missing'"):
            GenericMolecularGroups(definition)

    def test_unknown_configuration_field_is_rejected(self) -> None:
        definition = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        definition["groups"]["ethyl"]["atoms"][0]["lone_pair"] = 0

        with self.assertRaisesRegex(ValueError, "unbekannte Felder: lone_pair"):
            GenericMolecularGroups(definition)

    def test_component_rotation_preserves_internal_bond_lengths(self) -> None:
        definition = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        definition["molecules"]["ethanol"]["components"][0]["rotation_degrees"] = [90, 0, 0]
        molecule = GenericMolecularGroups(definition).create_molecule("ethanol", "rotated")

        carbon_bond = next(bond for bond in molecule.bonds if bond.bond_id.endswith("ethyl:C1-C2"))
        first_atom_id = carbon_bond.electrons[0].atom_id
        second_atom_id = carbon_bond.electrons[1].atom_id
        atoms_by_id = {atom.atom_id: atom for atom in molecule.atoms}
        distance = (
            sum(
                (atoms_by_id[first_atom_id].position[axis] - atoms_by_id[second_atom_id].position[axis]) ** 2
                for axis in range(3)
            )
            ** 0.5
        )

        self.assertAlmostEqual(distance, 1.54)


if __name__ == "__main__":
    unittest.main()
