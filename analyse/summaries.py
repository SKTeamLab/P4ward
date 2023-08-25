import pandas as pd
from pathlib import Path

def summary_csv(protac_objs):

    results_folder = Path("./results_summaries")
    results_folder.mkdir(exist_ok=True)

    for protac_obj in protac_objs:

        pose_objs = protac_obj.protein_poses
        # get protein pose attributes
        data_dict = {
            'pose_number':[],
            'megadock_score':[],
            'z_score':[],
            'cluster':[],
            'clrep':[],
        }
        for pose_obj in pose_objs:
            for attr in data_dict.keys():
                try:
                    attr_value = getattr(pose_obj, attr)
                except:
                    attr_value = None
                data_dict[attr].append(attr_value)
        
        # get protac pose and linkers attributes
        data_dict['protac_pose'] = []
        data_dict['active_linkers'] = []
        data_dict['top_protac_score'] = []

        for pose_obj in pose_objs:
            protac_pose = protac_obj.get_pose(pose_obj)
            data_dict['protac_pose'].append(protac_pose.active)
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


def chimerax_view(receptor_obj, protac_objs, pose_objs, generated_poses_folder, protac_poses_folder):
    """
    Make a chimerax visualization script to see the final successful poses
    """

    import seaborn as sns
    
    results_folder = Path("./results_summaries")
    results_folder.mkdir(exist_ok=True)

    for protac_obj in protac_objs:

        script = f'open ../{receptor_obj.file} name receptor;\ncolor ##name="receptor" #afafaf target asr;\n'
        palette = sns.color_palette('viridis', len(pose_objs)).as_hex()

        for i in range(len(protac_obj.protein_poses)):
            pose_obj = protac_obj.protein_poses[i]
            script += (
                f'open ../{generated_poses_folder}/pose{pose_obj.pose_number}.pdb name pose{pose_obj.pose_number};\n'
                +f'open ../{protac_poses_folder}/protac_{protac_obj.name}/protein_pose_{pose_obj.pose_number}/protac_embedded_confs.sdf name pose{pose_obj.pose_number}_protac;\n'
                +f'color ##name="pose{pose_obj.pose_number}" | ##name="pose{pose_obj.pose_number}_protac" {palette[i]} target asr;\n'
                +f'color ##name="pose{pose_obj.pose_number}_protac" byhet;\n'
            )
        
        script += (
            f"hide H;\n"
            +f"lighting full;\n"
            +f"set bgcolor white;\n"
            +f"view;\n"
        )

        with open(results_folder / f'summary-{protac_obj.name}.cxc', 'w+') as cmx:
            cmx.write(script)
