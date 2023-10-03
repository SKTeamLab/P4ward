from pathlib import Path
from ..tools import decorators
from ..tools.logger import logger
from ..tools.structure_tools import load_biopython_structures
from ..definitions import ROOT_DIR

@decorators.user_choice
@decorators.track_run
def ligand_distances(receptor_obj, ligase_obj, protac_objs):
    """
    Use Biopython to filter the megadock poses which satisfy a
    distance cutoff for both binding sites. Takes the a Protein object
    and handles the rest by accessing its attributes.
    """

    import numpy as np
    dist_cutoff = np.max([i.dist_cutoff for i in protac_objs])
    logger.info(f'Filtering megadock poses with cuttoff {dist_cutoff}')

    from ..tools.structure_tools import structure_proximity

    pose_objs = ligase_obj.active_confs()

    for pose_obj in pose_objs:
        receptor_lig_obj = receptor_obj.get_ligand_struct()
        ligand_lig_obj = pose_obj.get_rotated_struct(struct_type='ligand')

        distance, proximity = structure_proximity(
            receptor_lig_obj,
            ligand_lig_obj,
            dist_cutoff=dist_cutoff
        )

        if proximity:
            pose_obj.active = True
            pose_obj.filtered = True
            logger.info(f'Activating filtered pose {pose_obj.pose_number}, with distance of {distance} between ligands.')
        else:
            pose_obj.filtered = False
            pose_obj.active = False
    
    # make a protac_obj for each filtered 
    # for protac_obj in protac_objs:
    #     for pose_obj in pose_objs:
    #         if pose_obj.filtered:
    #             classes.ProtacPose(parent=protac_obj, protein_parent=pose_obj)

    # TODO if there are no filtered poses, quit the program


@decorators.user_choice
@decorators.track_run
def crl_filters(
                        receptor_obj,
                        ligase_obj,
                        crl_model_clash,
                        clash_threshold,
                        clash_count_tol,
                        accessible_lysines,
                        lys_sasa_cutoff,
                        overlap_dist_cutoff,
                        vhl_ubq_dist_cutoff,
                        crbn_ubq_dist_cutoff,
                        e3,
):
    
    from copy import deepcopy
    import numpy as np
    from Bio.PDB import Superimposer


    model_info = {
        'vhl':{
            'model_numbers' : [1],
            'ubq_point'     : np.array([-1.400, 47.989, 16.144]),
            'dist_cutoff'   : vhl_ubq_dist_cutoff
        },
        'crbn':{
            'model_numbers' : [1,2,3,4,5,9,10,11,12],
            'ubq_point'     : np.array([-33.329, -5.636, 15.041]),
            'dist_cutoff'   : crbn_ubq_dist_cutoff
        }
    }

    def get_e3_modelfile(e3, model_number=1, subrec_only=False):
        """
        return the a model file corresponding to the e3 (vhl or crbn) and the model number
        if subrec_only == True, we get the pdb that contains only the vhl or crbn proteins
        if False, we get the pdb which contains the full crl model
        """
        if subrec_only:
            filename = f'model{model_number}_{e3}.pdb'
        else:
            filename = f'model{model_number}.pdb' 
        
        modelfile =  ROOT_DIR / 'structures' / 'crl_models' / e3 / filename
        return(modelfile)
        
    def prep_alignment(modelfile, outfilename, ligase_obj):
        """
        Uses pymol to align the active ligase file, which was the one actually used for 
        docking to the ligase structure in the full CRL model. This way the operations that are
        performed downstream will use identical structures.
        """
            
        import pymol

        pymol.cmd.load(modelfile)
        pymol.cmd.load(ligase_obj.active_file)
        pymol.cmd.align(ligase_obj.active_file.stem, modelfile.stem)
        pymol.cmd.save(outfilename, ligase_obj.active_file.stem)

        return(outfilename)
    
    def check_model_clash(
            full_model_struct, rec_struct_align,
            clash_threshold, clash_count_tol
    ):
        
        from Bio.PDB import Selection, NeighborSearch

        ns = NeighborSearch(list(full_model_struct.get_atoms()))
        receptor_atoms = Selection.unfold_entities(rec_struct_align, 'A')

        clash_count = 0
        for atom in receptor_atoms:
            close_atoms = ns.search(atom.coord, clash_threshold)
            if len(close_atoms) > 0:
                # then this receptor atom is clashing with at least one atom for the model
                clash_count+=1
            if clash_count >= clash_count_tol:
                return(True)
                
        return(False)

    def check_access_lys(
            rec_struct_align,
            ubq,
            lys_sasa_cutoff,
            dist_cutoff,
            overlap_dist_cutoff
    ):
        
        # calculate the solvent exposed surface area for the residues in the rec
        from Bio.PDB.SASA import ShrakeRupley
        sr = ShrakeRupley()
        sr.compute(rec_struct_align, level='R')

        # now I must grab all the lysines in the receptor that have sasa > cutoff
        lysines = []
        for res in rec_struct_align.get_residues():
            if res.resname == 'LYS' and res.sasa > lys_sasa_cutoff:
                lysines.append(res)

        lysines_accepted = []
        for lysine in lysines:

            ## the lysine point will be its nitrogen:
            for atom in lysine.get_atoms():
                if atom.id == 'NZ':
                    x,y,z = atom.get_vector()
                    lys = np.array([x,y,z])
                    break
            
            # make atoms list for the protein that does not include the lysine itself
            atoms = []
            for atom in rec_struct_align.get_atoms():
                if not atom.get_parent()._id[1] == lysine._id[1]:
                    x,y,z = atom.get_vector()
                    atoms.append([x,y,z])
                else:
                    continue
            atoms = np.asarray(atoms)

            ubq2lys = lys - ubq
            ubq2atoms = atoms - ubq

            if np.linalg.norm(ubq2lys) > dist_cutoff:
                final_verdict = True

            else:
                ## project the receptor atoms to the line between ubq and lys
                proj = (np.sum(ubq2atoms * ubq2lys, axis=1) / np.dot(ubq2lys, ubq2lys))[:, np.newaxis] * ubq2lys
                proj_vector = ubq + proj

                proj_dist = np.linalg.norm(ubq2atoms - proj, axis=1)

                is_positive = np.dot((atoms - ubq), ubq2lys) > 0
                is_less_sqr = np.dot((atoms - ubq), ubq2lys) < np.linalg.norm(ubq2lys)**2
                in_segment = np.logical_and(is_positive, is_less_sqr)

                final_check = np.logical_and((proj_dist < overlap_dist_cutoff), in_segment)
                final_verdict = np.any(final_check)

            #the final verdict is True if there is an atom in the way, False if not
            lysines_accepted.append(not final_verdict)

        accessible_lys = np.any(lysines_accepted)
        lys_count = np.array(lysines_accepted).sum()
        return(accessible_lys, lys_count)


    receptor_struct = receptor_obj.get_protein_struct()

    # make aligned ligase file and load the structures for each crl model
    models_folder = Path('./crl_models')
    models_folder.mkdir(exist_ok=True)

    aligned_ligase_structs = {}

    for model_number in model_info[e3]['model_numbers']:
        model_file = get_e3_modelfile(e3, model_number, subrec_only=True)
        aligned_ligase_file = models_folder / f'model{model_number}_aligned_ligase.pdb'
        prep_alignment(model_file, str(aligned_ligase_file), ligase_obj)

        aligned_ligase_struct = load_biopython_structures(aligned_ligase_file)
        aligned_ligase_structs[model_number] = aligned_ligase_struct

    ##
    ## we check each model number for each pose_obj
    ##
    pose_objs = ligase_obj.active_confs()
    for pose_obj in pose_objs:

        pose_struct = pose_obj.get_rotated_struct('protein')

        ##              ##
        ## begin checks ##
        ##              ##        
        accepted_models = []

        for model_number in model_info[e3]['model_numbers']:

            model_file = get_e3_modelfile(e3, model_number, subrec_only=True)
            full_model_file = get_e3_modelfile(e3, model_number)
            aligned_ligase_struct = aligned_ligase_structs[model_number]

            # now that the ligase in CRL model conformation is identical to the ligase that was docked
            # (thanks to prep_alignment()), they can be aligned by biopython
            superimposer = Superimposer()
            superimposer.set_atoms(
                list(aligned_ligase_struct.get_atoms()),
                list(pose_struct.get_atoms())
            )    
            ## apply the rotation to the receptor
            rec_struct_align = deepcopy(receptor_struct)
            superimposer.apply(rec_struct_align)

            if crl_model_clash:

                full_model_struct = load_biopython_structures(full_model_file)
                clash = check_model_clash(
                    full_model_struct,
                    rec_struct_align,
                    clash_threshold,
                    clash_count_tol
                )

                if clash:
                    accepted_models.append(False)
                    pose_obj.crl = -1
                    continue

                else:
                    if accessible_lysines:
                        accessible_lys, lys_count = check_access_lys(
                            rec_struct_align,
                            ubq=model_info[e3]['ubq_point'],
                            lys_sasa_cutoff=lys_sasa_cutoff,
                            dist_cutoff=model_info[e3]['dist_cutoff'],
                            overlap_dist_cutoff=overlap_dist_cutoff
                        )
                        pose_obj.crl = lys_count

                        if accessible_lys:
                            accepted_models.append(True)
                        else:
                            accepted_models.append(False)
                    
                    else:
                        accepted_models.append(True)
                        
        accepted_pose = np.any(accepted_models)
        if accepted_pose:
            pose_obj.filtered = True
            pose_obj.active = True
        else:
            pose_obj.filtered = False
            pose_obj.active = False
        
        logger.debug(f'Pose number {pose_obj.pose_number} filtered: {pose_obj.filtered}')