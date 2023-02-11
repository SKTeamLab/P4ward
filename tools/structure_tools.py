from ..tools.logger import logger


def load_biopython_structures(protein=None, protein_ligand_chain=None, protein_ligand_resnum=None):
    """
    get a protein structure and the protein_ligand as biopython objects
    returns (protein_struct, protein_ligand_struct)
    """

    from Bio.PDB.PDBParser import PDBParser as pdbp
    parser = pdbp(PERMISSIVE=1)

    # protein and its ligand
    protein_struct = parser.get_structure('protein_struct', protein)

    # find ligand, if required:
    if protein_ligand_chain is not None and protein_ligand_resnum is not None:
        for i in protein_struct[0][protein_ligand_chain].get_residues():
            if i.get_id()[1] == protein_ligand_resnum:
                protein_ligand_struct = i
                logger.debug('found residue match')
            else: continue
        # make sure that the ligand was found:
        try: 
            protein_ligand_struct.get_id()[1] == protein_ligand_resnum
        except:
            err = f'protein ligand with resnum {protein_ligand_resnum} not found'
            logger.error(err)
            raise Exception(err)
    # if no ligand resnum and chain specified, we get just the protein:
    else:
        protein_ligand_struct = None

    # only return the structures that were asked for
    if protein_ligand_struct is not None:
        return(protein_struct, protein_ligand_struct)
    else:
        return(protein_struct)


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

    logger.debug(f'obj1 has {len(obj1_coords)} points')
    logger.debug(f'obj2 has {len(obj2_coords)} points')

    rmsd = rmsd.rmsd(obj1_coords, obj2_coords)
    return(rmsd)
