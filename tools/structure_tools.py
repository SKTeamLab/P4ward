from ..tools.logger import logger
import subprocess
import os


def load_biopython_structures(structure_file, mol2=False):
    """
    load a structure as biopython object. If file in mol2 format,
    convert first using rdkit. Intended for use with the ligands.
    """

    from Bio.PDB.PDBParser import PDBParser as pdbp
    parser = pdbp(PERMISSIVE=1)

    if mol2:
        from rdkit import Chem
        from io import StringIO

        rdkit_obj = Chem.MolFromMol2File(structure_file)
        pdb_block = Chem.MolToPDBBlock(rdkit_obj)
        pdbfile = StringIO(pdb_block)

        structure_obj = parser.get_structure('structure', pdbfile)

    else:
        structure_obj = parser.get_structure('structure', structure_file)

    return(structure_obj)


def structure_proximity(struct1, struct2, dist_cutoff=None):
    """
    Calculate the distance between the center of mass of two biopython struct objects
    and check if it is <= than a distance cutoff, returns Bool.
    """

    import numpy as np

    distance = np.linalg.norm(
        struct1.center_of_mass() - struct2.center_of_mass()
    )
    print(struct1.center_of_mass())
    print(struct2.center_of_mass())
    return(distance, distance <= dist_cutoff)


def get_rmsd(obj1, obj2, ca=False):
    """
    Use biopython to get the atomic coordinates of both objects
    and the pip tool rmsd.rmsd to calulate the rmsd. Returns rmsd value
    """

    #from Bio.PDB.PDBParser import PDBParser as pdbp
    import numpy as np
    import rmsd

    objects = {
        'obj1': obj1,
        'obj2': obj2,
        'obj1_coords':[],
        'obj2_coords':[]
    }
    
    for obj in ('obj1', 'obj2'):
        for res in objects[obj].get_residues():
            for atom in res.get_atoms():
                if ca:
                    if atom.get_name() == 'CA':
                        objects[f'{obj}_coords'].append(atom.get_coord())
                    else: continue
                else:
                    objects[f'{obj}_coords'].append(atom.get_coord())
    
    obj1_coords = np.asarray(objects['obj1_coords'])
    obj2_coords = np.asarray(objects['obj2_coords'])

    rmsd = rmsd.rmsd(obj1_coords, obj2_coords)
    return(rmsd)


def reduce(protein_obj_list, protein_only=False):
    """
    use chimerax to add hydrogens to proteins
    """

    for protein_obj in protein_obj_list:
        protein_file = protein_obj.file
        reduced_protein_file = os.path.splitext(protein_file)[0] + '_h.pdb'

        if protein_only:
            del_nonstd = 'del ~protein & H;'
        else: 
            del_nonstd = ''

        command = (
             f"open {protein_file};"
            + "addh;"
            +  del_nonstd
            +f"save {os.path.join(reduced_protein_file)}"
        )
        subprocess.run(['chimerax', '--nogui'], input=command, encoding='ascii')

        protein_obj.reduced_file = reduced_protein_file
        logger.info(f"Added hydrogens to {protein_file}")


def smiles2smarts(smiles_code):
    """
    turn a smiles code into a smarts code that is a wildcard,
    all atoms that are not carbon are changed into a *
    """
    smiles_code = smiles_code.split('\t')[0]

    atomlist = ['s','h','o','n','p','f','i','Cl','Br']

    for atom in atomlist:
        smiles_code = smiles_code.replace(atom, '*')
        smiles_code = smiles_code.replace(atom.upper(), '*')
        smiles_code = smiles_code.replace('.', '')
    
    return(smiles_code)