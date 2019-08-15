import layout as cpp
import numpy as np

__all__ = ['layout','layout_convergent','layout_sparse','mds_direct','draw_svg','draw_png']

def layout(I, J, V=None, t_max=15, eps=.1):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress."""

    # initialize positions
    X = random_init(I, J)
    
    if V is None:
        cpp.layout_unweighted(X, I, J, t_max, eps)
    else:
        cpp.layout_weighted(X, I, J, V, t_max, eps)
    return X

def layout_convergent(I, J, V=None, t_max=30, eps=.1, delta=.03, t_maxmax=200):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress
    at a guaranteed stationary point."""

    # initialize positions
    X = random_init(I, J)

    if V is None:
        cpp.layout_unweighted_convergent(X, I, J, t_max, eps, delta, t_maxmax)
    else:
        cpp.layout_weighted_convergent(X, I, J, V, t_max, eps, delta, t_maxmax)

    return X

def layout_sparse(I, J, npivots, V=None, t_max=30, eps=.01):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress
    using the sparse approximation of Ortmann et al. (2017)"""

    # initialize positions
    X = random_init(I, J)

    if (npivots > X.shape[0]):
        raise "number of pivots exceeds number of vertices"

    if V is None:
        cpp.layout_sparse_unweighted(X, I, J, npivots, t_max, eps);
    else:
        cpp.layout_sparse_weighted(X, I, J, V, npivots, t_max, eps);

    return X


def mds_direct(n, d, w, etas):
    """takes nC2 vectors d and w with a vector of step sizes eta
    and returns a n-by-2 matrix of positions X"""

    # initialize positions
    X = np.random.rand(n, 2)
    nC2 = (n*(n-1))/2
    if len(d) != nC2 or len(w) != nC2:
        raise "d and w are not correct length condensed distance matrices"

    cpp.mds_direct(X, d, w, eta)
    return X



### no c++ bindings for functions below ###

def random_init(I, J):
    if len(I) != len(J):
        raise "length of edge indices I and J not equal"

    n = max(max(I), max(J)) + 1
    X = np.random.rand(n,2)
    return X

#def draw_png(X, I, J, filepath, noderadius=.2, linkwidth=.05, width=100, border=50, linkopacity=1):
#    """Takes a n-by-2 matrix of positions X and index pairs I and J
#    and draws it to a .png file.
#    It is a direct conversion of draw_svg(), but using the cairosvg
#    package to convert to render the image. It may therefore be slow
#    for large graphs."""
#
#    svg = draw_svg(X,I,J,None,noderadius,linkwidth,width,border,linkopacity)
#    import cairosvg as cairo
#    cairo.svg2png(bytestring=svg.encode('ASCII'), write_to=filepath)
    

def draw_svg(X, I, J, filepath=None, noderadius=.2, linkwidth=.05, width=1000, border=50, linkopacity=1):
    """Takes a n-by-2 matrix of positions X and index pairs I and J
    and writes the equivalent picture in svg format.
    The drawing will be expanded into a width*width square
    Note that the parameters are in svg pixel units.
    The style at the top of the output svg may also be edited as necessary.
    The svg is returned as a string if filename is empty."""

    n = len(X)
    m = len(I)

    X_min = [min(X[i,0] for i in range(n)), min(X[i,1] for i in range(n))]
    X_max = [max(X[i,0] for i in range(n)), max(X[i,1] for i in range(n))]

    range_max = max(X_max[0]-X_min[0], X_max[1]-X_min[1]) # taller or wider
    range_max += 2*noderadius # guarantee no nodes are cut off at the edges
    scale = (width-2*border) / range_max

    X_svg = np.empty((n,2))
    for i in range(n):
        X_svg[i] = (X[i] - X_min) * scale
        X_svg[i] += [border + scale*noderadius, border + scale*noderadius]

    svg_list = []
    svg_list.append('<svg width="{:.0f}" height="{:.0f}" xmlns="http://www.w3.org/2000/svg">'.format(width, width))
    svg_list.append('<style type="text/css">')
    svg_list.append('line{{stroke:black;stroke-width:{:.3f};stroke-opacity:{:.3f};stroke-linecap:round;}}'.format(scale*linkwidth,linkopacity))
    svg_list.append('circle{{r:{};fill:black}}'.format(scale*noderadius))
    svg_list.append('</style>');

    for ij in range(m):
        i = I[ij]
        j = J[ij]
        X_i = X_svg[i]
        X_j = X_svg[j]
        svg_list.append('<line x1="{:.1f}" x2="{:.1f}" y1="{:.1f}" y2="{:.1f}"/>'.format(X_i[0], X_j[0], X_i[1], X_j[1]));

    if noderadius > 0:
        # draw nodes
        for i in range(n):
            svg_list.append('<circle cx="{:.1f}" cy="{:.1f}"/>'.format(X_svg[i][0], X_svg[i][1]))

    svg_list.append("</svg>")

    if filepath == None or filepath == '':
        return '\n'.join(svg_list)
    else:
        f = open(filepath, 'w')
        f.write('\n'.join(svg_list))
        f.close()
