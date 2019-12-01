from .swig import layout as cpp
import numpy as np

__all__ = ['layout','layout_convergent','layout_sparse','mds_direct','default_schedule','draw_svg','draw_png']

def layout(I, J, V=None, t_max=30, eps=.01, random_seed=None, init=None):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress."""

    if len(I) == 0 or len(I) != len(J):
        raise ValueError("length of edge indices I and J not equal or zero")

    # seed random state
    random_seed = _check_random_seed(random_seed)

    X = random_init(I, J, random_seed, init)
    
    if V is None:
        cpp.layout_unweighted(X, I, J, t_max, eps, random_seed)
    else:
        cpp.layout_weighted(X, I, J, V, t_max, eps, random_seed)
    return X

def layout_convergent(I, J, V=None, t_max=30, eps=.01, delta=.03, t_maxmax=200, random_seed=None, init=None):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress
    at a guaranteed stationary point."""

    if len(I) == 0 or len(I) != len(J):
        raise ValueError("length of edge indices I and J not equal or zero")

    # seed random state
    random_seed = _check_random_seed(random_seed)

    X = random_init(I, J, random_seed, init)

    if V is None:
        cpp.layout_unweighted_convergent(X, I, J, t_max, eps, delta, t_maxmax, random_seed)
    else:
        cpp.layout_weighted_convergent(X, I, J, V, t_max, eps, delta, t_maxmax, random_seed)

    return X

def layout_sparse(I, J, npivots, V=None, t_max=30, eps=.01, random_seed=None, init=None):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress
    using the sparse approximation of Ortmann et al. (2017)"""

    # seed random state
    random_seed = _check_random_seed(random_seed)

    X = random_init(I, J, random_seed, init)

    if (npivots > X.shape[0]):
        raise ValueError("number of pivots exceeds number of vertices")

    if V is None:
        cpp.layout_sparse_unweighted(X, I, J, npivots, t_max, eps, random_seed)
    else:
        cpp.layout_sparse_weighted(X, I, J, V, npivots, t_max, eps, random_seed)

    return X


def mds_direct(n, d, w=None, etas=None, random_seed=None, init=None):
    """takes nC2 vectors d (distance) and w (weight) with a vector of step sizes eta
    and returns a n-by-2 matrix of positions X"""

    nC2 = int((n*(n-1))/2)
    if w is None:
        w = np.ones(nC2) # standard metric MDS

    if len(d) != nC2 or len(w) != nC2:
        raise ValueError("d and/or w are not correct length condensed distance matrices")

    if etas is None:
        etas = default_schedule(w)

    # seed random state
    random_seed = _check_random_seed(random_seed)

    X = _random_init(n, random_seed, init)

    # do mds
    cpp.mds_direct(X, d, w, etas, random_seed)
    return X

def default_schedule(w, t_max=30, eps=.01):
    eta_max = 1/min(w)
    eta_min = eps/max(w)
    lambd = np.log(eta_max / eta_min) / (t_max-1)
    etas = eta_max * np.exp(-lambd * np.arange(t_max))
    return etas


### no c++ bindings for functions below ###

def _check_random_seed(random_seed):
    if random_seed is not None:
        random_seed = np.random.randint(65536)
    return random_seed

def _random_init(n, random_seed, init=None):
    # initialize positions
    if init is not None:
        if len(init) != n:
            raise ValueError("initial layout has incorrect shape")
        X = np.ascontiguousarray(init)
    else:
        np.random.seed(random_seed)
        X = np.random.rand(n, 2)
    return X


def random_init(I, J, random_seed, init=None):
    if len(I) == 0 or len(I) != len(J):
        raise ValueError("length of edge indices I and J not equal or zero")

    n = max(max(I), max(J)) + 1

    return _random_init(n, random_seed, init)


def draw_png(X, I, J, filepath, noderadius=.2, linkwidth=.05, width=1000, border=50, nodeopacity=1, linkopacity=1):
    """Takes a n-by-2 matrix of positions X and index pairs I and J
    and draws it to a .png file, using the PyCairo library."""

    n = len(X)
    m = len(I)

    X_min = [min(X[i,0] for i in range(n)), min(X[i,1] for i in range(n))]
    X_max = [max(X[i,0] for i in range(n)), max(X[i,1] for i in range(n))]

    range_max = max(X_max[0]-X_min[0], X_max[1]-X_min[1]) # taller or wider
    range_max += 2*noderadius # guarantee no nodes are cut off at the edges
    scale = (width-2*border) / range_max

    X_png = np.empty((n,2))
    for i in range(n):
        X_png[i] = (X[i] - X_min) * scale
        X_png[i] += [border + scale*noderadius, border + scale*noderadius]

    # use cairo to draw
    import cairo
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, width)
    ctx = cairo.Context(surface)

    ctx.set_line_width(linkwidth * scale)
    ctx.set_line_cap(cairo.LineCap.ROUND)
    ctx.set_source_rgba(0,0,0,linkopacity)

    # draw links
    for ij in range(m):
        i = I[ij]
        j = J[ij]
        X_i = X_png[i]
        X_j = X_png[j]
        ctx.move_to(X_i[0], X_i[1])
        ctx.line_to(X_j[0], X_j[1])
        ctx.stroke()

    # draw nodes
    if noderadius > 0 and nodeopacity > 0:
        ctx.set_source_rgba(0,0,0,nodeopacity)
        radius = noderadius * scale
        for i in range(n):
            ctx.arc(X_png[i][0], X_png[i][1], radius, 0, 7)
            ctx.fill()

    surface.write_to_png(filepath)


def draw_svg(X, I, J, filepath=None, noderadius=.2, linkwidth=.05, width=1000, border=50, nodeopacity=1, linkopacity=1):
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
    svg_list.append('circle{{r:{};fill:black;fill-opacity:{:.3f}}}'.format(scale*noderadius,nodeopacity))
    svg_list.append('</style>')

    # draw links
    for ij in range(m):
        i = I[ij]
        j = J[ij]
        X_i = X_svg[i]
        X_j = X_svg[j]
        svg_list.append('<line x1="{:.1f}" x2="{:.1f}" y1="{:.1f}" y2="{:.1f}"/>'.format(X_i[0], X_j[0], X_i[1], X_j[1]))

    # draw nodes
    if noderadius > 0:
        for i in range(n):
            svg_list.append('<circle cx="{:.1f}" cy="{:.1f}"/>'.format(X_svg[i][0], X_svg[i][1]))

    svg_list.append("</svg>")

    if filepath == None or filepath == '':
        return '\n'.join(svg_list)
    else:
        f = open(filepath, 'w')
        f.write('\n'.join(svg_list))
        f.close()
