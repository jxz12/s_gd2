import layout
import numpy as np

def test(X):
    layout.test(X)

def layout_unweighted(n, I, J, t_max=15, eps=.1):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress"""

    # initialize positions
    X = np.random.rand(n, 2)
    
    layout.layout_unweighted(X, I, J, t_max, eps)
    return X

def layout_weighted(n, I, J, V, t_max=15, eps=.1):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress"""

    # initialize positions
    X = np.random.rand(n, 2)

    layout.layout_weighted(X, I, J, V, t_max, eps)
    return X

def mds_direct(n, d, w, eta):
    """takes nC2 vectors d and w, and a vector of step sizes eta
    and returns a n-by-2 matrix of positions X with minimized stress"""

    # initialize positions
    X = np.random.rand(n, 2)

    layout.mds_direct(X, d, w, eta)
    return X


# def layout_scipy_focus(S, f, weighted=False, t_max=20, eps=.01):
#     """takes a scipy sparse matrix S with n vertices
#     and returns a n-by-2 matrix of positions X with minimized stress
#     and a focus on vertex f
#     usually requires larger t_max and lower eps to guarantee a tidy layout"""
# 
#     n,_ = S.shape
#     if f < 0 or f >= n:
#         raise Exception("focus has index not within bounds")
# 
#     d, w = init_dw(S, weighted)
#     eta = init_eta(w, t_max, eps)
# 
#     # set relevant weights in distance vector to effectively infinity for focus
#     w_f = 1/min(eta)/min(w)
#     
#     ij = f-1
#     for i in range(f):
#         w[ij] = w_f
#         ij += n-i-2
#     ij += 1
#     for i in range(f, n-1):
#         w[ij] = w_f
#         ij += 1
# 
#     # initialize positions
#     X = np.random.rand(n, 2)
# 
#     # perform SGD
#     layout.sgd_direct(X, d, w, eta)
#     return X

# def layout_scipy_y_constrained(S, Y, weighted=False, t_max=15, eps=.1):
#     """takes a scipy sparse matrix S with n vertices
#     and returns a n-by-2 matrix of positions X with minimized stress
#     with y-axis positions constrained to the values in Y"""
# 
#     n,_ = S.shape
#     d, w = init_dw(S, weighted)
#     eta = init_eta(w, t_max, eps)
# 
#     # initialize positions
#     X = np.random.rand(n, 2)
# 
#     # constrain y axis positions
#     for i in range(n):
#         X[i,1] = Y[i];
# 
#     # perform SGD only on x-axis
#     layout.sgd_direct_horizontal(X, d, w, eta)
#     return X

def concentric_circles(n_circ, radius=1, n_seg=100, offset=[0,0]):
    """returns coordinates to plot concentric circles using pyplot.plot in a numpy matrix
    where the first dimension is which circle, the second is x/y axis, third is the coordinates
    an optional offset may also be used"""

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

