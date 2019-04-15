import layout as cpp
import numpy as np

def layout(n, I, J, V=None, t_max=15, eps=.1):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress"""

    # initialize positions
    X = np.random.rand(n, 2)
    
    if V is None:
        cpp.layout_unweighted(X, I, J, t_max, eps)
    else:
        cpp.layout_weighted(X, I, J, V, t_max, eps)
    return X

def layout_convergent(n, I, J, V=None, t_max=30, eps=.1, delta=.03, t_maxmax=200):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress
    at a guaranteed stationary point."""

    # initialize positions
    X = np.random.rand(n, 2)

    if V is None:
        cpp.layout_unweighted_convergent(X, I, J, t_max, eps, delta, t_maxmax)
    else:
        cpp.layout_weighted_convergent(X, I, J, V, t_max, eps, delta, t_maxmax)

    return X

def layout_focus(n, I, J, focus, V=None, t_max=30, eps=.1, delta=.03, t_maxmax=200):
    """takes a list of indices I and J with single index f
    and returns a n-by-2 matrix of positions X with a focus on node f
    at a guaranteed stationary point."""

    # initialize positions
    X = np.random.rand(n, 2)

    if V is None:
        cpp.layout_unweighted_focus(X, I, J, focus, t_max, eps, delta, t_maxmax)
    else:
        cpp.layout_weighted_focus(X, I, J, focus, V, t_max, eps, delta, t_maxmax)

    return X

def layout_horizontal(n, I, J, Y, V=None, t_max=30, eps=.1, delta=.03, t_maxmax=200):
    """takes a list of indices I and J with vertical positions y
    and returns a n-by-2 matrix of positions X with a focus on node f
    at a guaranteed stationary point."""

    # initialize positions
    X = np.random.rand(n, 2)
    X[:,1] = Y

    if V is None:
        cpp.layout_unweighted_horizontal(X, I, J, t_max, eps, delta, t_maxmax)
    else:
        cpp.layout_weighted_horizontal(X, I, J, V, t_max, eps, delta, t_maxmax)
    return X

def mds_direct(n, d, w, eta):
    """takes nC2 vectors d and w with a vector of step sizes eta
    and returns a n-by-2 matrix of positions X"""

    # initialize positions
    X = np.random.rand(n, 2)

    cpp.mds_direct(X, d, w, eta)
    return X



### helper functions below, no c++ bindings ###

def draw_svg(X, I, J, filepath=None, noderadius=0, linkwidth=3, scale=30, border=30, linkopacity=1):
    """Takes a n-by-2 matrix of positions X and index pairs I and J
    and returns a string in svg format.
    The style at the top of the output svg may also be edited as necessary.
    Note that the parameters are in svg pixel units.
    The svg is printed on the command line if filename is empty"""

    n = len(X)
    m = len(I)

    X_min = [min(X[i,0] for i in range(n)), min(X[i,1] for i in range(n))]
    X_max = [max(X[i,0] for i in range(n)), max(X[i,1] for i in range(n))]

    X_svg = np.empty((n,2))
    for i in range(n):
        X_svg[i] = (X[i] - X_min) * scale
        X_svg[i] += [border, border]

    width = (X_max[0] - X_min[0])*scale + 2*border
    height = (X_max[1] - X_min[1])*scale + 2*border

    svg_list = []
    svg_list.append('<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg">'.format(width, height))
    svg_list.append('<style type="text/css">')
    svg_list.append('line{{stroke:black;stroke-width:{};stroke-opacity:{};stroke-linecap:round;}}'.format(linkwidth,linkopacity))
    svg_list.append('circle{{r:{}}}'.format(noderadius))
    svg_list.append('</style>');

    for ij in range(m):
        i = I[ij]
        j = J[ij]
        X_i = X_svg[i]
        X_j = X_svg[j]
        svg_list.append('<line x1="{}" x2="{}" y1="{}" y2="{}"/>'.format(X_i[0], X_j[0], X_i[1], X_j[1]));

    if noderadius > 0:
        # draw nodes
        for i in range(n):
            svg_list.append('<circle cx="{}" cy="{}"/>'.format(X_svg[i][0], X_svg[i][1]))

    svg_list.append("</svg>")

    if filepath == None or filepath == '':
        print('\n'.join(svg_list))
    else:
        f = open(filepath, 'w')
        f.write('\n'.join(svg_list))
        f.close()


def concentric_circles(num_circles, offset=[0,0], radius=1, n_seg=100):
    """Returns coordinates to plot concentric circles using pyplot.plot in a numpy matrix
    where the first dimension is which circle, the second is x/y axis, third is the coordinates
    an optional offset may also be used. Useful for layout_focus()."""

    circles = np.empty(shape=(num_circles, 2, n_seg+1))
    for circ in range(num_circles):
        r = circ+1 * radius
        for seg in range(n_seg):
            angle = seg/n_seg * 2*np.pi
            circles[circ,0,seg] = r * np.cos(angle) + offset[0]
            circles[circ,1,seg] = r * np.sin(angle) + offset[1]

        # join the end back to the start
        circles[circ,0,n_seg] = circles[circ,0,0]
        circles[circ,1,n_seg] = circles[circ,1,0]

    return circles

