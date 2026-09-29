import pandas as pd
import plotly.express as px
import plotly.io as pio


def test_summary_dataframe_creation():
    """Verify that Pandas can construct and filter summary tables for poses"""
    data = {
        "pose_number": [1, 2, 3],
        "megadock_score": [1200.5, 950.2, 810.0],
        "crl": [True, True, False],
        "cluster_centr": [True, False, False],
    }
    df = pd.DataFrame(data)

    assert len(df) == 3

    # filter by cluster representative
    top_poses = df[df["cluster_centr"] == True]
    assert len(top_poses) == 1
    assert top_poses.iloc[0]["pose_number"] == 1


def test_plotly_funnel_html_generation():
    """Verify that px builds a funnel figure and exports to html"""
    data = {
        "stages": ["Initial Poses", "Distance Filter", "Clustered"],
        "counts": [1000, 250, 15],
    }
    fig = px.funnel(data, x="counts", y="stages", title="Pose Funnel Test")
    html_output = pio.to_html(fig)

    assert "Pose Funnel Test" in html_output
    assert "<html" in html_output.lower() or "<div" in html_output.lower()
