import layout
import numpy as np

def draw(n, I, J, V=None, t_max=15, eps=.1):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress"""

    # initialize positions
    X = np.random.rand(n, 2)
    
    if V is None:
        layout.layout_unweighted(X, I, J, t_max, eps)
    else:
        layout.layout_weighted(X, I, J, V, t_max, eps)
    return X

def draw_convergent(n, I, J, V=None, t_max=30, eps=.1, delta=.03, t_maxmax=200):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress
    at a guaranteed stationary point."""

    # initialize positions
    X = np.random.rand(n, 2)

    if V is None:
        layout.layout_unweighted_convergent(X, I, J, t_max, eps, delta, t_maxmax)
    else:
        layout.layout_weighted_convergent(X, I, J, V, t_max, eps, delta, t_maxmax)

    return X

def draw_focus(n, I, J, focus, V=None, t_max=30, eps=.1, delta=.03, t_maxmax=200):
    """takes a list of indices I and J with single index f
    and returns a n-by-2 matrix of positions X with a focus on node f
    at a guaranteed stationary point."""

    # initialize positions
    X = np.random.rand(n, 2)

    if V is None:
        layout.layout_unweighted_focus(X, I, J, focus, t_max, eps, delta, t_maxmax)
    else:
        layout.layout_weighted_focus(X, I, J, focus, V, t_max, eps, delta, t_maxmax)

    return X

def draw_horizontal(n, I, J, Y, V=None, t_max=30, eps=.1, delta=.03, t_maxmax=200):
    """takes a list of indices I and J with vertical positions y
    and returns a n-by-2 matrix of positions X with a focus on node f
    at a guaranteed stationary point."""

    # initialize positions
    X = np.random.rand(n, 2)
    X[:,1] = Y

    if V is None:
        layout.layout_unweighted_horizontal(X, I, J, t_max, eps, delta, t_maxmax)
    else:
        layout.layout_weighted_horizontal(X, I, J, V, t_max, eps, delta, t_maxmax)
    return X

def mds_direct(n, d, w, eta):
    """takes nC2 vectors d and w with a vector of step sizes eta
    and returns a n-by-2 matrix of positions X"""

    # initialize positions
    X = np.random.rand(n, 2)

    layout.mds_direct(X, d, w, eta)
    return X


def concentric_circles(n_circ, offset=[0,0], radius=1, n_seg=100):
    """Returns coordinates to plot concentric circles using pyplot.plot in a numpy matrix
    where the first dimension is which circle, the second is x/y axis, third is the coordinates
    an optional offset may also be used. Useful for layout_focus()."""

    circles = np.empty(shape=(n_circ, 2, n_seg+1))
    for circ in range(n_circ):
        r = circ+1 * radius
        for seg in range(n_seg):
            angle = seg/n_seg * 2*np.pi
            circles[circ,0,seg] = r * np.cos(angle) + offset[0]
            circles[circ,1,seg] = r * np.sin(angle) + offset[1]

        # join the end back to the start
        circles[circ,0,n_seg] = circles[circ,0,0]
        circles[circ,1,n_seg] = circles[circ,1,0]

    return circles

