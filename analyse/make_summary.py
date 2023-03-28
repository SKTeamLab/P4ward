import pandas as pd

def summary_csv(pose_objects):

    data_dict = {
        'pose_number':[],
        'megadock_score':[],
        'z_score':[],
        'cluster':[],
        'centroid':[],
        'linker_gen':[],
        'active_linkers':[],
    }
    
    for pose_obj in pose_objects:
        for attr in data_dict.keys():
            try:
                attr_value = getattr(pose_obj, attr)
            except:
                attr_value = None
            data_dict[attr].append(attr_value)
    
    data = pd.DataFrame.from_dict(data_dict)
    data.to_csv('summary.csv')