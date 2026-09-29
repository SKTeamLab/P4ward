from pathlib import Path
from p4ward.run.structure_tools import load_biopython_structures
from p4ward.tools.classes import Protein

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_load_biopython_structure():
    """Verify that biopython can load and parse a pdb file into atom coordinates"""
    pdb_path = REPO_ROOT / "tutorial" / "receptor.pdb"
    struct = load_biopython_structures(str(pdb_path))

    atoms = list(struct.get_atoms())
    assert len(atoms) > 0

    coord = atoms[0].get_coord()
    assert len(coord) == 3


def test_protein_class_initialization():
    """Verify that Protein class initializes and keeps proper file references"""
    protein = Protein(
        ptn_type="receptor",
        file=REPO_ROOT / "tutorial" / "receptor.pdb",
        lig_file=REPO_ROOT / "tutorial" / "receptor_ligand.mol2",
    )

    assert protein.type == "receptor"
    assert protein.file.name == "receptor.pdb"
    assert protein.lig_file.name == "receptor_ligand.mol2"
    assert protein.active_file == protein.file
