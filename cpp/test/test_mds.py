import numpy as np
import s_gd2


def pdist(X):
    n = X.shape[0]
    i_idx = np.concatenate([np.repeat(i, n - i - 1) for i in range(n)])
    j_idx = np.concatenate([np.arange(i + 1, n) for i in range(n)])
    return np.sqrt(np.sum((X[i_idx] - X[j_idx]) ** 2, axis=1))


def test_mds():
    X = np.random.normal(0, 1, (100, 10))
    D = pdist(X)
    X_mds = s_gd2.mds_direct(X.shape[0], D)
    assert X_mds.shape == (X.shape[0], 2)
    # X_mds2 = s_gd2.mds_direct(X.shape[0], D, init=X_mds
    # np.testing.assert_allclose(X_mds, X_mds2)


def test_mds_3D():
    X = np.random.normal(0, 1, (100, 10))
    D = pdist(X)
    X_mds = s_gd2.mds_direct(X.shape[0], D, num_dimensions=3)
    assert X_mds.shape == (X.shape[0], 3)
    # X_mds2 = s_gd2.mds_direct(X.shape[0], D, init=X_mds, num_dimensions=3)
    # np.testing.assert_allclose(X_mds, X_mds2)


def test_mds_seed():
    X = np.random.normal(0, 1, (100, 10))
    D = pdist(X)
    X_mds = s_gd2.mds_direct(X.shape[0], D, random_seed=42)
    X_mds2 = s_gd2.mds_direct(X.shape[0], D, random_seed=42)
    np.testing.assert_array_equal(X_mds, X_mds2)
    X_mds2 = s_gd2.mds_direct(X.shape[0], D, random_seed=41)
    assert not np.all(X_mds == X_mds2)
