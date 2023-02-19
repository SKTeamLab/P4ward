import pandas as pd

def summary_csv(pose_objects):

    data_dict = {
        'pose_number':[],
        'megadock_score':[],
        'z_score':[],
        'cluster':[],
        'centroid':[]
    }
    
    for pose_obj in pose_objects:
        for attr in data_dict.keys():
            data_dict[attr].append(getattr(pose_obj, attr))
    
    data = pd.DataFrame.from_dict(data_dict)
    data.to_csv('summary.csv')