from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry.rdGeometry import Point3D


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

reclig = Chem.MolFromMol2File('/home/paula/Documents/projects/protacs/benchmark-2303/7KHH/receptor_ligand.mol2', sanitize=False, cleanupSubstructures=False)
liglig = Chem.MolFromMol2File('/home/paula/Documents/projects/protacs/benchmark-2303/7KHH/ligase_ligand.mol2', sanitize=False, cleanupSubstructures=False)
reference_ligs = Chem.CombineMols(liglig, reclig)

protac = Chem.MolFromSmiles('Cc1c(scn1)c2ccc(cc2)CNC(=O)[C@@H]3C[C@H](CN3C(=O)[C@H](C(C)(C)C)NC(=O)CCCCCCCCCCNC(=O)c4cc5c(cc4CS(=O)(=O)C)C6=CN(C(=O)c7c6c(c[nH]7)CN5c8c(cc(cn8)F)F)C)O', sanitize=False)
protac = Chem.AddHs(protac)

pose_lig_smiles = Chem.MolToSmiles(liglig, isomericSmiles=False, canonical=False)
receptor_lig_smiles = Chem.MolToSmiles(reclig, isomericSmiles=False, canonical=False)

# for each of their smiles, transform into wildcard smarts
receptor_lig_smarts_ = smiles2smarts(receptor_lig_smiles)
pose_lig_smarts_ = smiles2smarts(pose_lig_smiles)
receptor_lig_smarts = Chem.MolFromSmarts(receptor_lig_smarts_)
pose_lig_smarts = Chem.MolFromSmarts(pose_lig_smarts_)

# find atoms in the protac that match the wildcards:
receptor_lig_indices = protac.GetSubstructMatches(receptor_lig_smarts)[0]
pose_lig_indices = protac.GetSubstructMatches(pose_lig_smarts)[0]
# find atoms in the reference structure that match the wildcards,
# they are going to give the coordinates we have to match:
receptor_lig_coords = reference_ligs.GetSubstructMatches(receptor_lig_smarts)[0]
pose_lig_coords = reference_ligs.GetSubstructMatches(pose_lig_smarts)[0]
