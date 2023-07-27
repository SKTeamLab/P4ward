import pandas as pd

def summary_csv(pose_objects):

    # get protein pose attributes
    data_dict = {
        'pose_number':[],
        'megadock_score':[],
        'z_score':[],
        'cluster':[],
        'clrep':[],
    }
    for pose_obj in pose_objects:
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

    for pose_obj in pose_objects:
        data_dict['protac_pose'].append(pose_obj.protac_pose.active)
        if pose_obj.protac_pose.active:
            active_linkers = [i for i in pose_obj.protac_pose.linker_confs if i.active]
            if len(active_linkers) == 0:
                active_linkers = '0'
                top_protac_score = None
            else:
                top_protac_score = active_linkers[0].rx_score
                active_linkers = ','.join([str(i.conf_number) for i in active_linkers])
        else:
            active_linkers = '0'
            top_protac_score = None
        data_dict['active_linkers'].append(active_linkers)
        data_dict['top_protac_score'].append(top_protac_score)
    
    data = pd.DataFrame.from_dict(data_dict)
    data.to_csv('summary.csv')


def chimerax_view(receptor_obj, pose_objs):
    """
    Make a chimerax visualization script to see the final successful poses
    """
    
    import seaborn as sns

    script = f"open {receptor_obj.file};\n color #1 #afafaf target asr;\n"
    model_num = 1

    models = [i for i in pose_objs if i.protac_pose.active and len(i.protac_pose.linker_confs) > 0]
    models = sorted(models, key=lambda x: getattr(x, 'megadock_score'), reverse=True)
    palette = sns.color_palette('viridis', len(models)).as_hex()

    for i in range(len(models)):
        pose_obj = models[i]
        model_num += 2
        script += (
             f"open protein_docking/pose{pose_obj.pose_number}.pdb;\n"
            +f"open ligand_sampling/protein_pose_{pose_obj.pose_number}/protac_embedded_confs.sdf;\n"
            +f"color #{model_num-1},{model_num} {palette[i]} target asr;\n"
            +f"color #{model_num} byhet;\n"
        )
        for i in pose_obj.protac_pose.linker_confs:
            if not i.active:
                script += (f'close #{model_num} &  ##name="conf_{i.conf_number}";\n')
    
    script += (
         f"hide H;\n"
        +f"lighting full;\n"
        +f"set bgcolor white;\n"
        +f"view;\n"
    )

    with open('summary.cxc', 'w+') as cmx:
        cmx.write(script)
