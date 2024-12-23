import logging
from pathlib import Path
from multiprocessing import Process, JoinableQueue
from ..tools import decorators
from .structure_tools import load_biopython_structures, pymol_align
from ..definitions import ROOT_DIR

logger = logging.getLogger('p4ward')

@decorators.user_choice
@decorators.track_run
def ligand_distances(receptor_obj, ligase_obj, protac_objs):
    """
    Use Biopython to filter the megadock poses which satisfy a
    distance cutoff for both binding sites. Takes the a Protein object
    and handles the rest by accessing its attributes.
    """

    import numpy as np
    from .megadock import rotate_atoms

    dist_cutoff = np.max([i.dist_cutoff for i in protac_objs])
    logger.info(f'Filtering megadock poses with cuttoff {dist_cutoff}')

    ref_rotate = ligase_obj.rotate
    # ligase_com = ligase_obj.get_protein_struct().center_of_mass()
    ligase_lig_com = ligase_obj.get_ligand_struct().center_of_mass()
    receptor_lig_com = receptor_obj.get_ligand_struct().center_of_mass()

    pose_objs = ligase_obj.active_confs()

    for pose_obj in pose_objs:

        pose_rotate = pose_obj.rotate
        pose_lig_com = rotate_atoms(
            tuple(ligase_lig_com),
            ref_rotation=ref_rotate,
            pose_rotation=pose_rotate
        )

        distance = np.linalg.norm(receptor_lig_com - pose_lig_com)
        proximity = distance <= dist_cutoff
        logger.debug(f'Pose {pose_obj.pose_number}: distance {distance}, filtered: {proximity}')

        pose_obj.active = proximity
        pose_obj.filtered = proximity
        #        ^ need this attr for management even though the info is going to filter_info dict
        pose_obj.filter_info['dist_filter'] = proximity
        pose_obj.filter_info['distance'] = distance

    results = [i for i in pose_objs if i.filtered]
    if any(results):
        logger.info(f'Finished filtering {len(results)} protein poses according to ligand distance criteria.')
    else:
        logger.info('There are no poses which satisfy the ligand distance filtering criteria. Exiting now.')
        exit(0)
    


@decorators.user_choice
@decorators.track_run
def crl_filters(
                        receptor_obj,
                        ligase_obj,
                        crl_model_clash,
                        clash_threshold,
                        clash_count_tol,
                        accessible_lysines,
                        lysine_count,
                        lys_sasa_cutoff,
                        overlap_dist_cutoff,
                        vhl_ubq_dist_cutoff,
                        crbn_ubq_dist_cutoff,
                        e3,
                        num_procs,
                        pose_objs=None
):
    
    from copy import deepcopy
    import numpy as np
    from Bio.PDB import Superimposer
    from .structure_tools import get_e3_modelfile
    from ..structures.model_info import model_info

    model_info['vhl']['dist_cutoff'] = vhl_ubq_dist_cutoff
    model_info['crbn']['dist_cutoff'] = crbn_ubq_dist_cutoff

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
            
            # make atoms list for the protein that does not include the lysine itself or hydrogens
            atoms = []
            for atom in rec_struct_align.get_atoms():
                if not atom.get_parent()._id[1] == lysine._id[1] and not atom.element == 'H':
                    x,y,z = atom.get_vector()
                    atoms.append([x,y,z])
                else:
                    continue
            atoms = np.asarray(atoms)

            lys_dist = lys - ubq
            lys_dist_scalar = np.linalg.norm(lys_dist)

            if lys_dist_scalar > dist_cutoff:
                lys_is_accessible = False

            else:

                atoms_dist = lys - atoms
                lys_dist_unit = lys_dist / lys_dist_scalar
                ## project the receptor atoms to the line between ubq and lys
                atoms_proj_scalars = np.dot(atoms_dist, lys_dist_unit).reshape(-1,1)
                atoms_proj = lys_dist_unit * atoms_proj_scalars
                ## check distance of atoms to the line
                lys_occlusions = atoms_proj - atoms_dist
                occlusion_scalars = np.linalg.norm(lys_occlusions, axis=1).reshape(-1,1)
                ## find which atoms (if any) are in segment and too close to segment therefore occluding the lys
                atoms_in_segment = (atoms_proj_scalars < lys_dist_scalar) & (atoms_proj_scalars > 0)
                atoms_occluding = atoms_in_segment & (occlusion_scalars < overlap_dist_cutoff)

                if atoms_occluding.any():
                    lys_is_accessible = False
                else:
                    lys_is_accessible = True

            if lys_is_accessible:
                lys = {'resname':lysine._id[1], 'distance':lys_dist_scalar, 'sasa':lysine.sasa, 'accessible':True}
                lysines_accepted.append(lys)

        lys_count = len(lysines_accepted)
        return(lys_count, lysines_accepted)

    receptor_struct = receptor_obj.get_protein_struct()

    # make aligned ligase file and load the structures for each crl model
    models_folder = Path('./crl_models')
    models_folder.mkdir(exist_ok=True)

    aligned_ligase_structs = {}

    for model_number in model_info[e3]['model_numbers']:
        model_file = get_e3_modelfile(e3, model_number, subrec_only=True)
        aligned_ligase_file = models_folder / f'model{model_number}_aligned_ligase.pdb'

        pymol_align(
            target_file=model_file,
            moving_file=ligase_obj.active_file,
            outfilename=str(aligned_ligase_file)
        )
        # ^ Use pymol to align the active ligase file, which was the one actually used for 
        # docking to the ligase structure in the full CRL model. This way the operations that are
        # performed downstream will use identical structures.

        aligned_ligase_struct = load_biopython_structures(aligned_ligase_file)
        aligned_ligase_structs[model_number] = aligned_ligase_struct

    ##
    ## we check each model number for each pose_obj
    ##

    def worker(pose_obj):
        """
        The worker function in this case takes a pose_obj and does its filtering for each
        crl model number
        """

        pose_struct = pose_obj.get_rotated_struct('protein')

        ##              ##
        ## begin checks ##
        ##              ##  
        accepted_models = []
        crls = []
        lys_info = []

        for model_number in model_info[e3]['model_numbers']:

            model_file = get_e3_modelfile(e3, model_number, subrec_only=True)
            full_model_file = get_e3_modelfile(e3, model_number)
            aligned_ligase_struct = aligned_ligase_structs[model_number]

            # now that the ligase in CRL model conformation is identical to the ligase that was docked
            # they can be aligned by biopython
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
                    crl = -1
                    continue

                else:
                    crl = None
                    if accessible_lysines:
                        lys_count, lysines_accepted = check_access_lys(
                            rec_struct_align,
                            ubq=model_info[e3]['ubq_point'],
                            lys_sasa_cutoff=lys_sasa_cutoff,
                            dist_cutoff=model_info[e3]['dist_cutoff'],
                            overlap_dist_cutoff=overlap_dist_cutoff
                        )
                        crl = lys_count

                        if lys_count >= lysine_count:
                            accepted_models.append(True)
                        else:
                            accepted_models.append(False)
                        lys_info.append(lysines_accepted)
                    
                    else:
                        accepted_models.append(True)
                
                crls.append(crl)

        accepted_pose = np.any(accepted_models)
        pose_obj.filter_info['crl_filter'] = accepted_pose
        pose_obj.filter_info['crls'] = crls
        pose_obj.filter_info['lys_info'] = lys_info

        if accepted_pose:
            pose_obj.filtered = True
            pose_obj.active = True
        else:
            pose_obj.filtered = False
            pose_obj.active = False



        # while True:

        #     pose_obj = inQ.get()
        #     pose_struct = pose_obj.get_rotated_struct('protein')

        #     ##              ##
        #     ## begin checks ##
        #     ##              ##  
        #     accepted_models = []
        #     crls = []
        #     lys_info = []

        #     for model_number in model_info[e3]['model_numbers']:

        #         model_file = get_e3_modelfile(e3, model_number, subrec_only=True)
        #         full_model_file = get_e3_modelfile(e3, model_number)
        #         aligned_ligase_struct = aligned_ligase_structs[model_number]

        #         # now that the ligase in CRL model conformation is identical to the ligase that was docked
        #         # they can be aligned by biopython
        #         superimposer = Superimposer()
        #         superimposer.set_atoms(
        #             list(aligned_ligase_struct.get_atoms()),
        #             list(pose_struct.get_atoms())
        #         )    
        #         ## apply the rotation to the receptor
        #         rec_struct_align = deepcopy(receptor_struct)
        #         superimposer.apply(rec_struct_align)


        #         if crl_model_clash:

        #             full_model_struct = load_biopython_structures(full_model_file)
        #             clash = check_model_clash(
        #                 full_model_struct,
        #                 rec_struct_align,
        #                 clash_threshold,
        #                 clash_count_tol
        #             )

        #             if clash:
        #                 accepted_models.append(False)
        #                 crl = -1
        #                 continue

        #             else:
        #                 crl = None
        #                 if accessible_lysines:
        #                     lys_count, lysines_accepted = check_access_lys(
        #                         rec_struct_align,
        #                         ubq=model_info[e3]['ubq_point'],
        #                         lys_sasa_cutoff=lys_sasa_cutoff,
        #                         dist_cutoff=model_info[e3]['dist_cutoff'],
        #                         overlap_dist_cutoff=overlap_dist_cutoff
        #                     )
        #                     crl = lys_count

        #                     if lys_count >= lysine_count:
        #                         accepted_models.append(True)
        #                     else:
        #                         accepted_models.append(False)
                        
        #                 else:
        #                     accepted_models.append(True)
                    
        #             crls.append(crl)
        #             lys_info.append(lysines_accepted)

        #     accepted_pose = np.any(accepted_models)
        #     outQ.put((pose_obj.pose_number, accepted_pose, crls, lys_info))


    ######################

    if pose_objs == None:
        pose_objs = ligase_obj.active_confs()

    done_count = 0
    for pose_obj in pose_objs:

        worker(pose_obj)
        done_count += 1
        
        logger.debug(f'Pose number {pose_obj.pose_number} filtered: {pose_obj.filtered}')

        if done_count == len(pose_objs):
            break
    
    results = [i.filtered for i in pose_objs]
    if any(results):
        logger.info(f'Finished filtering {len(results)} protein poses according to CRL models.')
    else:
        logger.info('There are no poses which satisfy the CRL model filtering criteria. Exiting now.')
        exit(0)


    # inQ = JoinableQueue()
    # outQ = JoinableQueue()

    # for i in pose_objs:
    #     inQ.put(i)

    # procs = []
    # for i in range(num_procs):
    #     p = Process(name=i, target=worker, args=(inQ, outQ))
    #     p.daemon = True
    #     p.start()
    #     procs.append(p)

    # done_count = 0
    # while True:

    #     pose_number, accepted_pose, crls, lys_info = outQ.get()
    #     done_count += 1

    #     pose_obj = [i for i in pose_objs if i.pose_number == pose_number][0]
    #     pose_obj.filter_info['crl_filter'] = accepted_pose
    #     pose_obj.filter_info['crls'] = crls
    #     pose_obj.filter_info['lys_info'] = lys_info
    #     # pose_obj.crl = crl
    #     if accepted_pose:
    #         pose_obj.filtered = True
    #         pose_obj.active = True
    #     else:
    #         pose_obj.filtered = False
    #         pose_obj.active = False
        
    #     outQ.task_done()
    #     logger.debug(f'Pose number {pose_obj.pose_number} filtered: {pose_obj.filtered}')

    #     if done_count == len(pose_objs):
    #         break
    
    # results = [i.filtered for i in pose_objs]
    # if any(results):
    #     logger.info(f'Finished filtering {len(results)} protein poses according to CRL models.')
    # else:
    #     logger.info('There are no poses which satisfy the CRL model filtering criteria. Exiting now.')
    #     exit(0)