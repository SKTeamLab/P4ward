import logging
from pathlib import Path
import numpy as np
import plotly.express as px
import plotly.io as pio
from ..tools import decorators
from ..tools.classes import *

logger = logging.getLogger('p4ward')

plot_colors = ['#34345e', '#548687', '#FCAA67', '#eee220']


def make_step_check(conf:object, protac_obj:Protac) -> dict:
    """
    Checks which modelling steps were run so that their results can be plotted (or not).
    """

    steps = {
        'all'         : conf.getboolean('megadock', 'run_docking'),
        'cl_red'      : conf.getboolean('protein_ranking', 'cluster_poses_redundancy'),
        'dist_filter' : conf.getboolean('protein_filter', 'ligand_distances'),
        'crl_noclash' : conf.getboolean('protein_filter', 'crl_model_clash'),
        'crl_lys'     : conf.getboolean('protein_filter', 'accessible_lysines'),
        'protac'      : conf.getboolean('linker_sampling', 'rdkit_sampling'),
        'cl_trend'    : conf.getboolean('protein_ranking', 'cluster_poses_trend') and hasattr(protac_obj, 'cluster'),
        'filter'      : (
            conf.getboolean('protein_filter', 'ligand_distances') or 
            conf.getboolean('protein_filter', 'crl_model_clash') or
            conf.getboolean('protein_filter', 'rdkit_sampling')
        )
    }
    return(steps)


def plot_funnel(ligase_obj:Protein, protac_obj:Protac, steps:dict) -> str:
    """
    Make the funnel interactive plot.
    This plot shows how many protein poses passed each filtering stage.
    """

    plot_labels = ['all','dist_filter','cl_red','crl_noclash','crl_lys','protac','cl_trend']
    labels = [i for i in plot_labels if steps[i]]
    crl_filtered = [i for i in ligase_obj.conformations if i.filter_info.get('crls') is not None]
    all_counts = {
        'all'         : len(ligase_obj.conformations),
        'dist_filter' : len([i for i in ligase_obj.conformations if i.filter_info['dist_filter']]),
        'cl_red'      : (ligase_obj.cluster.clusterer.n_clusters_ if hasattr(ligase_obj,'cluster') else None),
        'crl_noclash' : len([i for i in crl_filtered if any(j >= 0 for j in i.filter_info['crls'])]),
        'crl_lys'     : len([i for i in crl_filtered if any(j >= 1 for j in i.filter_info['crls'])]),
        'protac'      : len(protac_obj.protein_poses),
        'cl_trend'    : (protac_obj.cluster.clusterer.n_clusters_ if hasattr(protac_obj, 'cluster') else None)
    }
    counts = {i:all_counts[i] for i in plot_labels if steps[i]}

    data = {'labels':counts.keys(), 'counts':counts.values()}

    fig = px.funnel(
        data,x='counts',y='labels',
        title='Number of Filtered Poses',
        color_discrete_sequence=[plot_colors[0]]
    )
    fig.update_layout(
        # width=500, height=500,
        template="ggplot2"
    )
    # fig.write_html('funnel.html')
    fig_html = pio.to_html(fig)

    return(fig_html)


def plot_ppi(ligase_obj:Protein, protac_obj:Protac, steps:dict) -> str:
    """
    Make the ppi distribution interactive plot.
    This shows, for each main filtering step, what is the distribution of the ppi scores.
    """

    import plotly.figure_factory as ff

    plot_labels = ['all', 'filter', 'protac', 'cl_trend']
    labels = [i for i in plot_labels if steps[i]]

    all_scores = {
        'all'         : [i.megadock_score for i in ligase_obj.conformations],
        'filter'      : [i.megadock_score for i in ligase_obj.active_confs()],
        'protac'      : [i.megadock_score for i in protac_obj.protein_poses]
    }
    if hasattr(protac_obj, 'cluster'):
        all_scores['cl_trend'] = [i.megadock_score for i in protac_obj.cluster.repr_centr]

    hist_data = []
    labels_ = []
    for i in all_scores.keys():
        if i in labels and len(all_scores[i])>1:
            hist_data.append(all_scores[i])
            labels_.append(i)

    fig = ff.create_distplot(
        hist_data, labels_, show_rug=False, bin_size=100,
        colors=plot_colors[::-1]
    )
    fig.update_layout(
        # width=500, height=500,
        template="ggplot2",
        title="Distribution of PPI scores by filter"
    )
    # fig.write_html('ppi.html')
    fig_html = pio.to_html(fig)

    return(fig_html)


def plot_pca(ligase_obj:Protein, receptor_obj:Protein, protac_obj:Protac, steps:dict) -> str:
    """
    Make PCA plot.
    This runs PCA on the triad coordinates of the protein poses, then plots PC1xPC2.
    """

    import pandas as pd
    from sklearn.decomposition import PCA
    from ..run.structure_tools import get_coords_array

    pose_objs = ligase_obj.active_confs()
    coords = get_coords_array(pose_objs, ligase_obj)
    is_protac_pose = [i in protac_obj.protein_poses for i in pose_objs]
    if steps['cl_trend']:
        is_trend_repr = [i in protac_obj.cluster.repr_centr for i in pose_objs]
   
    poses = []
    for i in range(len(is_protac_pose)):
        if hasattr(protac_obj,'cluster') and is_trend_repr[i]:
            poses.append('centroid')
        elif is_protac_pose[i]:
            poses.append('protac_pose')
        elif steps['filter']:
            poses.append('filtered')
        else:
            poses.append('all')

    # get the receptor's coordinates
    a,c,d = receptor_obj.get_triad_points()
    rec_coords = [[*a, *c, *d]]

    # run PCA
    pca = PCA(n_components=2)
    pca.fit(coords)
    coords_pca = pca.transform(coords)
    rec_coords_pca = pca.transform(rec_coords)
    # organize data
    pca_data = pd.DataFrame(coords_pca)
    # pca_data['cluster'] = protacs[0].cluster.clusterer.labels_
    pca_data['set'] = poses
    row = [rec_coords_pca[0][0], rec_coords_pca[0][1], 'receptor']
    pca_data.loc[len(pca_data)] = row

    # make fig

    fig = px.scatter(
        pca_data, x=0, y=1,
        color='set',
        color_discrete_sequence=['#c2c2c2','#eda62b', plot_colors[0], plot_colors[1]]
                # plot_colors = ['#34345e', '#548687', '#FCAA67', '#eee220']
    )
    fig.update_layout(
        # width=500, height=500,
        template="ggplot2",
        title='Protac interaction energies vs PPI'
    )
    # fig.write_html('pca.html')
    fig_html = pio.to_html(fig)
    return(fig_html)


def plot_scatter(protac_obj:Protac) -> str:
    """
    Make scatter plot that shows ppi score vs protein-protac score.
    The datapoints are coloured according to rescore, if rescore was not calculated,
    then they are colored according to ppi score.
    """

    import pandas as pd

    rescore = len(protac_obj.active_poses()) > 0 and hasattr(protac_obj.active_poses()[0], 'rescore') 

    data_dict = {
        'pose_number':[],
        'ppi':[],
        'interaction':[],
    }
    if rescore:
        data_dict['rescore'] = []
    for protac_pose in protac_obj.active_poses():
        interactions = [i.rx_score for i in protac_pose.active_confs()]

        data_dict['pose_number'].append(protac_pose.protein_parent.pose_number)
        data_dict['ppi'].append(protac_pose.protein_parent.megadock_score)
        data_dict['interaction'].append(np.mean(interactions))
        if rescore:
            data_dict['rescore'].append(protac_pose.rescore)
    
    data_scores = pd.DataFrame.from_dict(data_dict)

    if rescore:
        color = 'rescore'
    else:
        color = None

    fig = px.scatter(
        data_frame=data_scores,
        x='ppi',
        y='interaction',
        color=color,
        hover_data=data_scores.columns,
        color_continuous_scale='Cividis'
    )
    fig.update_layout(
        # width=500, height=500,
        template="ggplot2",
        title='Protac calculated scores'
    )
    # fig.write_html('protacs.html')
    fig_html = pio.to_html(fig)
    return(fig_html)

@decorators.user_choice
def interactive_plots(
        protacs:list[Protac],
        ligase_obj:Protein,
        receptor_obj:Protein,
        conf:object
) -> None:
    """
    Runs the functions that generate the html for the interactive plots:
    plot_funnel, plot_ppi, plot_pca, plot_scatter.
    Then generates an overarching html file that combines these plots and writes
    this html to disk.
    """
    
    results_folder = Path("./results_summaries")

    for protac_obj in protacs:

        # get run steps
        steps = make_step_check(conf, protac_obj)
        if not steps['all']:
            logger.info(f'Cannot generate plots for protac {protac_obj.name} because no docking was run.')
            return(None)

        # make plots
        fig_funnel_html  = plot_funnel(ligase_obj, protac_obj, steps)
        fig_ppi_html     = plot_ppi(ligase_obj, protac_obj, steps)
        fig_pca_html     = plot_pca(ligase_obj, receptor_obj, protac_obj, steps)
        fig_scatter_html = plot_scatter(protac_obj)

        # combine plots
        html_template = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>P4ward Result Plots</title>
    <!-- <script src="https://cdn.plot.ly/plotly-latest.min.js"></script> -->
</head>
<body>
    <div style="display: flex; flex-wrap: wrap; justify-content: space-around;">
        <div style="width: 45%; margin: 10px;">{plot1}</div>
        <div style="width: 45%; margin: 10px;">{plot2}</div>
        <div style="width: 45%; margin: 10px;">{plot3}</div>
        <div style="width: 45%; margin: 10px;">{plot4}</div>
    </div>
</body>
</html>
        """

        combined_html = html_template.format(
            plot1=fig_funnel_html,
            plot2=fig_ppi_html,
            plot3=fig_pca_html,
            plot4=fig_scatter_html
        )

        protac_file = results_folder / f"plots-{protac_obj.name}.html"
        with open(protac_file, 'w') as file:
            file.write(combined_html)