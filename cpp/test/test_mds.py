import numpy as np
from scipy.spatial.distance import pdist
import s_gd2

def test_mds():
    X = np.random.normal(0, 1, (100, 10))
    D = pdist(X)
    X_mds = s_gd2.mds_direct(X.shape[0], D)
    assert X_mds.shape == (X.shape[0], 2)
    X_mds2 = s_gd2.mds_direct(X.shape[0], D, init=X_mds)
    np.testing.assert_allclose(X_mds, X_mds2)

def test_mds_seed():
    X = np.random.normal(0, 1, (100, 10))
    D = pdist(X)
    X_mds = s_gd2.mds_direct(X.shape[0], D, random_seed=42)
    X_mds2 = s_gd2.mds_direct(X.shape[0], D, random_seed=42)
    np.testing.assert_array_equal(X_mds, X_mds2)
    X_mds2 = s_gd2.mds_direct(X.shape[0], D, random_seed=41)
    assert not np.all(X_mds == X_mds2)
   