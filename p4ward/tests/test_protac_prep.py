from rdkit import Chem
from p4ward.run.protac_prep import make_indices_link, check_linker_size
from p4ward.tools.classes import Protac


def test_make_indices_link():
    """Verify that linker heavy atoms are correctly identified between ligands"""
    # Propane: C0 - C1 - C2. If atom 0 and atom 2 are ligands, atom 1 is the linker.
    mol = Chem.MolFromSmiles("CCC")
    linker_indices = make_indices_link(mol, indices_ligs=[0, 2])

    assert linker_indices == [1]


def test_check_linker_size_extension():
    """Verify that small linkers extend flexibility to adjacent neighbor atoms"""
    # Pentane: C0 - C1 - C2 - C3 - C4. Linker is atom 2 (length 1).
    mol = Chem.MolFromSmiles("CCCCC")
    matches = {
        "receptor_lig_indices": [0, 1],
        "pose_lig_indices": [3, 4],
        "receptor_lig_coords": [0, 1],
        "pose_lig_coords": [3, 4],
    }

    # With min_linker_length=2, flexibility should expand into neighbors
    new_matches = check_linker_size(
        mol,
        indices_ligs=[0, 1, 3, 4],
        min_linker_length=2,
        neighbour_number=1,
        matches=matches,
    )

    # One atom from each ligand side was converted to flexible linker
    assert new_matches["receptor_lig_indices"] == [0]
    assert new_matches["pose_lig_indices"] == [4]


def test_protac_conformer_sampling():
    """Verify that rdkit can generate 3D conformers and compute mmff energies"""
    protac = Protac(smiles="CCO", name="ethanol", number=1)
    protac.sample_unbound_confs(num_unbound_confs=5)

    assert protac.num_confs == 5
    assert protac.unbound_energy is not None
    assert protac.unbound_energy > 0
