import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import NearestCentroid


def test_hierarchical_clustering_ward():
    """Verify that sklearn hierarchical clustering groups nearby 3D points"""
    # Two distinct clusters of 3D points: one near (0,0,0) and one near (10,10,10)
    coords = np.array([
        [0.0, 0.0, 0.0],
        [0.1, 0.1, 0.1],
        [10.0, 10.0, 10.0],
        [10.2, 10.1, 9.9],
    ])

    clusterer = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=1.0,
        metric="euclidean",
        linkage="ward",
    ).fit(coords)

    # Expect exactly 2 clusters
    assert clusterer.n_clusters_ == 2

    # First two points share a label, second two points share another label
    assert clusterer.labels_[0] == clusterer.labels_[1]
    assert clusterer.labels_[2] == clusterer.labels_[3]
    assert clusterer.labels_[0] != clusterer.labels_[2]

    # Compute centroids using NearestCentroid
    nc = NearestCentroid().fit(coords, clusterer.labels_)
    assert len(nc.centroids_) == 2
