"""Startskript für ein NaCl-Ionenmodell mit sichtbaren Chlorid-Elektronen."""

from dirabo.chem.particles import Molecule, MoleculeCluster
from dirabo.chem.psystems import chlorine, sodium
from dirabo.chem.render import show_molecules

# NaCl-Gitterabstand (Natrium-Chlor-Nachbarn), in den Modell-Angström.
ION_DISTANCE = 2.82
HALF_ION_DISTANCE = ION_DISTANCE / 2


def create_sodium_chloride_cluster() -> MoleculeCluster:
    """Erzeugt zwei Na+ und zwei Cl- mit vier sichtbaren Elektronenpaaren je Cl-."""
    ion_specs = (
        ("Na+", sodium, (-HALF_ION_DISTANCE, -HALF_ION_DISTANCE, 1.0), 1.0, 0),
        ("Cl-", chlorine, (HALF_ION_DISTANCE, -HALF_ION_DISTANCE, -1.0), 7.0, 4),
        ("Cl-", chlorine, (-HALF_ION_DISTANCE, 3, 0.0), 7.0, 4),
        ("Na+", sodium, (-HALF_ION_DISTANCE, HALF_ION_DISTANCE, 0.75), 1.0, 0),
        # ("Na+", sodium, (-HALF_ION_DISTANCE, -HALF_ION_DISTANCE, 0.0), 1.0, 0),
        # ("Cl-", chlorine, (HALF_ION_DISTANCE, -HALF_ION_DISTANCE, 0.0), 7.0, 4),
        # ("Cl-", chlorine, (-HALF_ION_DISTANCE, HALF_ION_DISTANCE, 0.0), 7.0, 4),
        # ("Na+", sodium, (HALF_ION_DISTANCE, HALF_ION_DISTANCE, 0.0), 1.0, 0),
    )
    ions = []
    for index, (ion_name, factory, position, nucleus_charge, lone_pair_count) in enumerate(ion_specs, start=1):
        atom_id = f"{ion_name.replace('+', 'plus').replace('-', 'minus')}-{index}"
        ion = factory(
            atom_id,
            position,
            lone_pair_count=lone_pair_count,
            nuclear_charge_override=nucleus_charge,
        )
        ions.append(
            Molecule(
                name=f"{ion_name} {index}",
                atoms=(ion.atom,),
                bonds=(),
                lone_pairs=ion.lone_pairs,
            )
        )

    return MoleculeCluster(tuple(ions))


if __name__ == "__main__":
    show_molecules(create_sodium_chloride_cluster())
