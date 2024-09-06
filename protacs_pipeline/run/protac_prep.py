import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS
from ..tools.logger import logger


def get_matches(protac, reclig, liglig, reference_ligs):

    # args:
    #   protac: as read by rdkit, local
    #   reclig: as read by rdkit, local
    #   liglig: as read by rdkit, local
    #   reference_ligs: combined mols from rdkit, local

    receptor_lig_smarts_ = rdFMCS.FindMCS([protac, reclig]).smartsString
    pose_lig_smarts_ = rdFMCS.FindMCS([protac, liglig]).smartsString
    receptor_lig_smarts = Chem.MolFromSmarts(receptor_lig_smarts_)
    pose_lig_smarts = Chem.MolFromSmarts(pose_lig_smarts_)
    
    # find atoms in the protac that match the wildcards:
    receptor_lig_indices = protac.GetSubstructMatches(receptor_lig_smarts)[0]
    pose_lig_indices = protac.GetSubstructMatches(pose_lig_smarts)[0]
    # find atoms in the reference structure that match the wildcards,
    # they are going to give the coordinates we have to match:
    receptor_lig_coords = reference_ligs.GetSubstructMatches(receptor_lig_smarts)[0]
    pose_lig_coords = reference_ligs.GetSubstructMatches(pose_lig_smarts)[0]

    matches = {
        'receptor_lig_indices' : list(receptor_lig_indices),
        'pose_lig_indices'     : list(pose_lig_indices),
        'receptor_lig_coords'  : list(receptor_lig_coords),
        'pose_lig_coords'      : list(pose_lig_coords)
    }
    return(matches)


def make_indices_link(protac, indices_ligs):
    indices_link = []
    for atom in Chem.RemoveHs(protac).GetAtoms():
        ix = atom.GetIdx()
        if ix not in indices_ligs:
            indices_link.append(ix)
    return(indices_link)


def check_linker_size(protac, indices_ligs, min_linker_length, neighbour_number, matches):
   
    indices_link = make_indices_link(protac, indices_ligs)
    
    if len(indices_link) <= min_linker_length:
        logger.warning(f"Linker length is {len(indices_link)}, extending flexibility to neighbouring atoms.")

        for _ in range(neighbour_number):
        
            neighbours = []
            protac_noh = Chem.RemoveHs(protac)
            for linker_atom in indices_link:                    # for each linker atom:
                atom = protac_noh.GetAtomWithIdx(linker_atom)   # get the atom obj
                nbs = [i.GetIdx() for i in atom.GetNeighbors()] # get idx of the atom's neighbours (nbs)
                nbs = [i for i in nbs if i not in indices_link] # only keep nbs that are not in the linker already
                neighbours.extend(nbs)                          # add them to neighbour list
            neighbours = list(np.unique(neighbours))            # remove duplicates

            receptor_lig_indices = dict(zip(range(len(matches['receptor_lig_indices'])), matches['receptor_lig_indices']))
            receptor_lig_coords = dict(zip(range(len(matches['receptor_lig_coords'])), matches['receptor_lig_coords']))            
            for i in range(len(receptor_lig_indices)):    # remove nbs from receptor_lig_indices and coords
                if receptor_lig_indices[i] in neighbours:
                    receptor_lig_indices.pop(i)
                    receptor_lig_coords.pop(i)
            matches['receptor_lig_indices'] = list(receptor_lig_indices.values())
            matches['receptor_lig_coords'] = list(receptor_lig_coords.values())

            pose_lig_indices = dict(zip(range(len(matches['pose_lig_indices'])), matches['pose_lig_indices']))
            pose_lig_coords = dict(zip(range(len(matches['pose_lig_coords'])), matches['pose_lig_coords']))
            for i in range(len(pose_lig_indices)):        # remove nbs from pose_lig_indices and coords
                if pose_lig_indices[i] in neighbours:
                    pose_lig_indices.pop(i)
                    pose_lig_coords.pop(i)
            matches['pose_lig_indices'] = list(pose_lig_indices.values())
            matches['pose_lig_coords'] = list(pose_lig_coords.values())

            indices_ligs = list(matches['receptor_lig_indices'])
            indices_ligs.extend(list(matches['pose_lig_indices']))

            indices_link = make_indices_link()
        
        logger.info(f"There are now {len(indices_link)} atoms treated as flexible.")
    
    else:
        logger.info(f"Linker length is {len(indices_link)}, not extending flexibility.")
    
    return(matches)


def check_ligand_matches(protac_objs, receptor_obj, ligase_obj, rdkit_ligands_cleanup):

    from rdkit.Chem.Draw import rdMolDraw2D
    from ..definitions import CWD

    reclig = Chem.MolFromMol2File(str(receptor_obj.lig_file), sanitize=rdkit_ligands_cleanup, cleanupSubstructures=rdkit_ligands_cleanup)
    liglig = Chem.MolFromMol2File(str(ligase_obj.lig_file), sanitize=rdkit_ligands_cleanup, cleanupSubstructures=rdkit_ligands_cleanup)
    reference_ligs = Chem.CombineMols(liglig, reclig)

    logger.info(
        f"Testing ligand match between the protac smiles codes " +
        f"and the ligand structures at {receptor_obj.lig_file.name} and {ligase_obj.lig_file.name}"
    )

    for protac_obj in protac_objs:
        
        protac = Chem.MolFromSmiles(protac_obj.smiles)
        matches_ = get_matches(protac, reclig, liglig, reference_ligs)
        matches = [*matches_['receptor_lig_indices'], *matches_['pose_lig_indices']]

        indices_ligs = list(matches_['receptor_lig_indices'])
        indices_ligs.extend(list(matches_['pose_lig_indices']))

        # draw the protac with the highlighted matches
        drawer = rdMolDraw2D.MolDraw2DCairo(1366, 1366)
        drawer.drawOptions().addAtomIndices = True
        highlight_atoms = [atom for atom in matches]
        drawer.DrawMolecule(protac, highlightAtoms=highlight_atoms)
        drawer.FinishDrawing()

        with open(CWD / f'ligand_matches-{protac_obj.name}.png', 'wb') as img:
            img.write(drawer.GetDrawingText())
        
        indices_link = make_indices_link(protac, indices_ligs)

        logger.info(
            f"Wrote image for {protac_obj.name} matches at ligand_matches-{protac_obj.name}.png" +
            f"Number of linker atoms found: {len(indices_link)}, at indices: {indices_link}"
        )


#~~~~~~~~~~~~~~~~~~~

def protac_prep(
                    receptor_obj,
                    ligase_obj,
                    protac_obj,
                    extend_flexible_small_linker,
                    min_linker_length,
                    neighbour_number,
                    rdkit_ligands_cleanup
):

    # open receptor ligand - will not change, get its rotation info
    reclig = Chem.MolFromMol2File(str(receptor_obj.lig_file), sanitize=rdkit_ligands_cleanup, cleanupSubstructures=rdkit_ligands_cleanup)
    ref_rotation  = ligase_obj.rotate
    # open ligase ligand - raw initial position
    liglig = Chem.MolFromMol2File(str(ligase_obj.lig_file), sanitize=rdkit_ligands_cleanup, cleanupSubstructures=rdkit_ligands_cleanup)
    # combine reclig and liglig into a single molecule
    reference_ligs = Chem.CombineMols(liglig, reclig)
    # open protac and addHs
    protac = Chem.MolFromSmiles(protac_obj.smiles)
    protac = Chem.AddHs(protac)

    # get substructure matches between protac 2d structure and ligands' structure files
    matches = get_matches(protac, reclig, liglig, reference_ligs)

    # get the indices for the ligands, in protac 2d
    indices_ligs = list(matches['receptor_lig_indices'])
    indices_ligs.extend(list(matches['pose_lig_indices']))

    if extend_flexible_small_linker:
        # check if linker size is too small and we need to extend flexibility to ligs
        matches = check_linker_size(protac, indices_ligs, min_linker_length, neighbour_number, matches)
        indices_ligs = list(matches['receptor_lig_indices'])
        indices_ligs.extend(list(matches['pose_lig_indices']))
    
    # make instructions to later align sampled protac 3d with reference ligs
    alignment_pairs = [(matches['receptor_lig_indices'][i], matches['receptor_lig_coords'][i]) for i in range(len(matches['receptor_lig_indices']))]
    alignment_pairs.extend([(matches['pose_lig_indices'][i], matches['pose_lig_coords'][i]) for i in range(len(matches['pose_lig_indices']))])

    # update parameters and attributes
    ## save protac_obj ligand and linker atom indices if not done before
    if protac_obj.indices_ligs == None:
        protac_obj.indices_ligs = indices_ligs
    ## save necessary variables to the parameters dict
    parameters = {
        'ref_rotation'    : ref_rotation,
        'reclig'          : reclig,
        'liglig'          : liglig,
        'protac'          : protac,
        'matches'         : matches,
        'alignment_pairs' : alignment_pairs
    }
    return(parameters)
