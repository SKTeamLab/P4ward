from ..tools.logger import logger
from ..tools import decorators
from pathlib import Path
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
        from io import StringIO

        run = subprocess.run(
            ['obabel','-imol2',structure_file,'-opdb'],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True
        )
        pdb_block = run.stdout
        pdbfile = StringIO(pdb_block)

        structure_obj = parser.get_structure('structure', pdbfile)

    else:
        structure_obj = parser.get_structure('structure', structure_file)

    return(structure_obj)


@decorators.user_choice
@decorators.track_run
def get_protac_dist_cuttoff(
        protac_obj,
        reclig_file,
        liglig_file,
        dist_cutoff
):

    if dist_cutoff == 'auto':

        logger.info("Ligands distance cutoff set to automatic.")
        logger.info('Sampling unbound protac conformations to determine distace cutoff.')

        import numpy as np
        from rdkit import Chem
        from rdkit.Chem import rdFMCS, AllChem
        
        def center_of_mass(atom_indices, conf):

            total_mass = 0.0
            xs = []; ys = []; zs = []
            for j in atom_indices:
                atom = conf.GetAtomWithIdx(j)
                mass = atom.GetMass()
                total_mass += mass
                x,y,z = conf.GetConformer().GetAtomPosition(j)
                xs.append(mass*x)
                ys.append(mass*y)
                zs.append(mass*z)
            center = np.array([np.sum(xs)/total_mass, np.sum(ys)/total_mass, np.sum(zs)/total_mass])
            return(center)

        protac = Chem.MolFromSmiles(protac_obj.smiles)
        reclig = Chem.MolFromMol2File(str(reclig_file), sanitize=False, cleanupSubstructures=False)
        liglig = Chem.MolFromMol2File(str(liglig_file), sanitize=False, cleanupSubstructures=False)

        receptor_lig_smarts_ = rdFMCS.FindMCS([protac, reclig]).smartsString
        pose_lig_smarts_ = rdFMCS.FindMCS([protac, liglig]).smartsString
        receptor_lig_smarts = Chem.MolFromSmarts(receptor_lig_smarts_)
        pose_lig_smarts = Chem.MolFromSmarts(pose_lig_smarts_)
        receptor_lig_indices = protac.GetSubstructMatches(receptor_lig_smarts)[0]
        pose_lig_indices = protac.GetSubstructMatches(pose_lig_smarts)[0]

        protac = Chem.AddHs(protac)

        params = AllChem.ETKDGv3()
        AllChem.EmbedMultipleConfs(protac, numConfs=100, params=params)
       
        distances = []

        for i in range(100):

            conf = Chem.Mol(protac, confId=i)
            center_reclig = center_of_mass(receptor_lig_indices, conf)
            center_liglig = center_of_mass(pose_lig_indices, conf)

            distance = np.linalg.norm(center_liglig - center_reclig)
            distances.append(distance)

        cutoff = np.mean(distances)
        protac_obj.dist_cutoff = cutoff

        logger.info(f"Setting distance cutoff to {cutoff}")
    
    else:
        protac_obj.dist_cutoff = float(dist_cutoff)
        logger.info(f"Ligands distance cutoff set to {dist_cutoff}.")



def structure_proximity(struct1, struct2, dist_cutoff=None):
    """
    Calculate the distance between the center of mass of two biopython struct objects
    and check if it is <= than a distance cutoff, returns Bool.
    """

    import numpy as np

    distance = np.linalg.norm(
        struct1.center_of_mass() - struct2.center_of_mass()
    )
    return(distance, distance <= dist_cutoff)


def get_rmsd(obj1, obj2, ca=False):
    """
    Use biopython to get the atomic coordinates of both objects
    and the pip tool rmsd.rmsd to calulate the rmsd. Returns rmsd value
    """

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


def reduce(protein_obj_list, file_attribute_name, protein_only=False):
    """
    use chimerax to add hydrogens to proteins.
    takes the objects as arguments so that they can be updated on the fly
    `file_attribute_name` specifies which file attribute to capture from the obj
    """

    for protein_obj in protein_obj_list:
        protein_file = getattr(protein_obj, file_attribute_name)
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

        setattr(protein_obj, file_attribute_name+'_reduced', reduced_protein_file)
        logger.info(f"Added hydrogens to {protein_file}")

# TODO remove
def smiles2smarts(smiles_code):
    """
    turn a smiles code into a smarts code that is a wildcard,
    all atoms that are not carbon are changed into a *
    """
    smiles_code = smiles_code.split('\t')[0]

    atomlist = ['s','h','o','n','p','f','i','Cl','Br']

    for atom in atomlist:

        # replace atoms with number
        for i in range(9):
            smiles_code = smiles_code.replace(f'[{atom}-{i}]', '[*]')
            smiles_code = smiles_code.replace(f'[{atom.upper()}-{i}]', '[*]')

        smiles_code = smiles_code.replace(f'[{atom}-]', '[*]')
        smiles_code = smiles_code.replace(f'[{atom}-]', '[*]')
        smiles_code = smiles_code.replace(atom, '*')
        smiles_code = smiles_code.replace(f'[{atom.upper()}-]', '[*]')
        smiles_code = smiles_code.replace(f'[{atom.upper()}-]', '[*]')
        smiles_code = smiles_code.replace(atom.upper(), '*')
        smiles_code = smiles_code.replace('.', '')
    
    return(smiles_code)


def obabel_convert(file_path, input_format, output_format, split=False, split_folder=''):

    if split:
        split_str = '--split'
        new_file = os.path.join(split_folder, f'template.{output_format}')
    else:
        split_str = ''
        new_file = file_path.replace(input_format, output_format)

    command = [
        'obabel', f'-i{input_format}', file_path, split_str, '-d'
        f'-o{output_format}', '-O', new_file
    ]
    subprocess.run(command)
    return(new_file)


def pymol_combine(*args, out_filename='combined.pdb'):

    import pymol

    basenames = [Path(filename).stem for filename in args]
    for filename in args:
        pymol.cmd.load(filename)
    
    pymol.cmd.create('combined', ' '.join(basenames))
    pymol.cmd.save(out_filename, 'combined')
