
def plot_ptn_poses(ligase_obj):

    import plotly.graph_objects as go
    
    active_confs = ligase_obj.active_confs()

    counts = {
        'total': len(ligase_obj.conformations),
        'dist_filter': len([i for i in ligase_obj.conformations if i.filter_info['dist_filter']]),
        'crl_noclash': len([i for i in active_confs if any(j >= 0 for j in i.filter_info['crls'])]),
        'lys': len([i for i in active_confs if any(j >= 1 for j in i.filter_info['crls'])]),
        '1lys': len([i for i in active_confs if any(j == 1 for j in i.filter_info['crls'])]),
        '2lys': len([i for i in active_confs if any(j == 2 for j in i.filter_info['crls'])]),
        '3lys': len([i for i in active_confs if any(j == 3 for j in i.filter_info['crls'])]),
        'other_lys': len([i for i in active_confs if any(j >= 4 for j in i.filter_info['crls'])]),
    }

    labels = list(counts.keys())
    parents = ["", "total", "dist_filter", "crl_noclash", "lys", "lys", "lys", "lys"]
    values = list(counts.values())

    fig = go.Figure()

    fig.add_trace(go.Treemap(
        labels = labels,
        parents = parents,
        values = values,
        textinfo = "label+value+percent parent+percent entry+percent root",
        root_color="lightgrey"
    ))
    fig.write_html('results_summaries/protein_filter.html')
