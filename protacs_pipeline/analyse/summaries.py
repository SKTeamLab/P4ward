import pandas as pd
from pathlib import Path
from ..tools import decorators
from ..tools.logger import logger
from ..definitions import CWD

def summary_csv(protac_objs, ligase_obj, benchmark, cluster_trend):

    results_folder = Path("./results_summaries")
    results_folder.mkdir(exist_ok=True)

    for protac_obj in protac_objs:

        # pose_objs = [i for i in ligase_obj.active_confs() if i in protac_obj.protein_poses]
        # # ^ this ensures that the pose_objs list is still ranked as done by rank.protein_poses()
        pose_objs = protac_obj.protein_poses
        # ^ this ensures that the pose_objs list is ranked by rank.rescore() if it was used

        ## get general protein pose attributes
        data_dict = {
            'pose_number':[],
            'megadock_score':[],
            'crl':[]
        }

        if benchmark:
            data_dict = {**data_dict, **{
                'l_rms':[],
                'i_rms':[],
                'fnat':[],
                'capri_rank':[]
            }}

        for pose_obj in pose_objs:
            
            for attr in data_dict.keys():
                try:
                    if attr == 'crl':
                        attr_value = pose_obj.filter_info['crls']
                    # elif attr == 'cl_size':
                    #     cln = pose_obj.cluster_redund.get_cl_from_pose(pose_obj)
                    #     attr_value = pose_obj.cluster_redund.get_cl_size(cln)
                    # elif attr == 'cluster_centr':
                    #     if pose_obj in protac_obj.cluster.repr_centr:
                    #         attr_value = True
                    #     else:
                    #         attr_value = False
                    else:
                        attr_value = getattr(pose_obj, attr)
                except:
                    attr_value = None
                data_dict[attr].append(attr_value)


        ## get cluster trend attributes if clustering ran successfully
        if cluster_trend and hasattr(protac_obj.cluster, 'repr_centr'):

            data_dict['cluster_number'] = []
            data_dict['cluster_centr'] = []
            data_dict['cluster_best'] = []
            data_dict['cluster_size'] = []

            for pose_obj in pose_objs:

                repr_centr = pose_obj in protac_obj.cluster.repr_centr
                data_dict['cluster_centr'].append(repr_centr)

                repr_best = pose_obj in [i.protein_parent for i in protac_obj.cluster.repr_best]
                data_dict['cluster_best'].append(repr_best)

                cln = protac_obj.cluster.get_cl_from_pose(pose_obj)
                data_dict['cluster_number'].append(cln)

                if repr_centr or repr_best:
                    cl_size = protac_obj.cluster.get_cl_size(cln)
                    data_dict['cluster_size'].append(cl_size)
                else:
                    data_dict['cluster_size'].append(None)


        # get protac pose and linkers attributes
        data_dict['protac_pose'] = []
        data_dict['active_linkers'] = []
        data_dict['top_protac_score'] = []
        data_dict['final_score'] = []

        for pose_obj in pose_objs:

            protac_pose = protac_obj.get_pose(pose_obj)
            data_dict['protac_pose'].append(protac_pose.active)

            if hasattr(protac_pose, 'rescore'):
                data_dict['final_score'].append(protac_pose.rescore)
            else:
                data_dict['final_score'].append(None)
                
            if protac_pose.active:
                active_linkers = [i for i in protac_pose.linker_confs if i.active]
                if len(active_linkers) == 0:
                    active_linkers = None
                    top_protac_score = None
                else:
                    top_protac_score = active_linkers[0].rx_score
                    active_linkers = ','.join([str(i.conf_number) for i in active_linkers])
            else:
                active_linkers = None
                top_protac_score = None
            data_dict['active_linkers'].append(active_linkers)
            data_dict['top_protac_score'].append(top_protac_score)
        
        data = pd.DataFrame.from_dict(data_dict)
        data.to_csv(results_folder / f'summary-{protac_obj.name}.csv')


@decorators.user_choice
def chimerax_view(receptor_obj, protac_objs, generated_poses_folder, protac_poses_folder, benchmark=False, ref_ligase=None):
    """
    Make a chimerax visualization script to see the final successful poses
    """

    import seaborn as sns
    
    results_folder = Path("./results_summaries")
    results_folder.mkdir(exist_ok=True)

    for protac_obj in protac_objs:

        if hasattr(protac_obj.cluster, 'repr_centr'):
            pose_objs = protac_obj.cluster.repr_centr
        else:
            pose_objs = protac_obj.protein_poses

        script = f'open ../{receptor_obj.file} name receptor;\ncolor ##name="receptor" #afafaf target asr;\n'
        palette = sns.color_palette('viridis', len(pose_objs)).as_hex()

        for i in range(len(pose_objs)):
            pose_obj = pose_objs[i]
            script += (
                f'open ../{generated_poses_folder}/pose{pose_obj.pose_number}.pdb name pose{pose_obj.pose_number};\n'
                +f'open ../{protac_poses_folder}/protac_{protac_obj.name}/protein_pose_{pose_obj.pose_number}/protac_embedded_confs.sdf name pose{pose_obj.pose_number}_protac;\n'
                +f'color ##name="pose{pose_obj.pose_number}" | ##name="pose{pose_obj.pose_number}_protac" {palette[i]} target asr;\n'
                +f'color ##name="pose{pose_obj.pose_number}_protac" byhet;\n'
            )
        
        if benchmark:
            script += (
                 f'open ../{ref_ligase} name ref_ligase\n'
                + 'color ##name="ref_ligase" #b44b49 target asr\n'
                +f'transparency ~(##name="receptor" | ##name="ref_ligase") 50 target asr\n'
            )

            for pose_obj in protac_obj.protein_poses:
                if pose_obj.capri_rank in ['acceptable', 'medium', 'high']:
                    script += (
                        f'transparency ##name="pose{pose_obj.pose_number}" 0 target asr\n'
                    )

        script += (
            f"hide H;\n"
            +f"lighting full;\n"
            +f"set bgcolor white;\n"
            +f"view;\n"
        )

        with open(results_folder / f'summary-{protac_obj.name}.cxc', 'w+') as cmx:
            cmx.write(script)


def protac_summaries(protac_objs, cluster_trend):
    """
    Write protac cluster information
    """

    import numpy as np

    for protac_obj in protac_objs:

        # ratio = np.mean(cl_sizes)/cl_count
        n_ptn_poses = len(protac_obj.protein_poses)
        
        logger.info(f'-> Protac {protac_obj.name}:')
        logger.info(f'Number of protein poses = {n_ptn_poses}')
        
        
        if cluster_trend:

            cl_count = protac_obj.cluster.clusterer.n_clusters_
            scores = [i.megadock_score for i in protac_obj.cluster.get_all_confs()]
            avg_cl_size = n_ptn_poses/cl_count
        
            logger.info(f'Avg cluster size = {avg_cl_size}')
            logger.info(f'Number of clusters = {cl_count}')
            logger.info(f'Avg cluster size / NClusters = {avg_cl_size/cl_count}')
        
        else:

            scores = [i.megadock_score for i in protac_obj.protein_poses]
        
        mean_score = np.mean(scores)
        logger.info(f'Mean ptn score = {mean_score}')


@decorators.user_choice
def write_crl_complex(receptor_obj, protac_objs, e3, protac_poses_folder, linker_scoring_folder, cluster_rep_only=True):

    import pymol2
    from ..structures.model_info import model_info
    from ..run.structure_tools import get_e3_modelfile

    crl_folder = CWD / 'crl_models'

    for protac_obj in protac_objs:

        if cluster_rep_only and hasattr(protac_obj.cluster, 'repr_centr'):
            pose_objs = protac_obj.cluster.repr_centr
        else:
            pose_objs = protac_obj.protein_poses
        
        for pose_obj in pose_objs:

            output_folder = crl_folder / f'protac_{protac_obj.name}' / f'protein_pose{pose_obj.pose_number}'
            output_folder.mkdir(parents=True, exist_ok=True)

            for model_num in model_info[e3]['model_numbers']:

                crl_file = get_e3_modelfile(e3, model_number=model_num, subrec_only=False)
                receptor_file = receptor_obj.active_file
                pose_file = pose_obj.file

                protac_file = linker_scoring_folder / f'protac_{protac_obj.name}' / f'protein_pose_{pose_obj.pose_number}' / 'protac_scored_confs.sdf'
                if not protac_file.is_file():
                    protac_file = protac_poses_folder / f'protac_{protac_obj.name}' / f'protein_pose_{pose_obj.pose_number}' / 'protac_embedded_confs.sdf'
                
                output_file = output_folder / f'model{model_num}.pdb'

                with pymol2.PyMOL() as pm:

                    pm.cmd.load(crl_file, 'target')
                    pm.cmd.load(pose_file, 'moving')
                    pm.cmd.load(receptor_file, 'bring_receptor')
                    pm.cmd.load(protac_file, 'bring_protac')

                    pm.cmd.align('moving', 'target')
                    pm.cmd.matrix_copy('moving','bring_receptor')
                    pm.cmd.matrix_copy('moving', 'bring_protac', target_state=0)

                    pm.cmd.save(output_file, state=0)
