from .swig import layout as cpp
import numpy as np

__all__ = [
    "layout",
    "layout_convergent",
    "layout_sparse",
    "mds_direct",
    "default_schedule",
    "draw_svg",
    "draw_png",
]


def layout(I, J, V=None, t_max=30, eps=0.01, random_seed=None, init=None):
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


def layout_convergent(
    I,
    J,
    V=None,
    t_max=30,
    eps=0.01,
    delta=0.03,
    t_maxmax=200,
    random_seed=None,
    init=None,
):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress
    at a guaranteed stationary point."""

    if len(I) == 0 or len(I) != len(J):
        raise ValueError("length of edge indices I and J not equal or zero")

    # seed random state
    random_seed = _check_random_seed(random_seed)

    X = random_init(I, J, random_seed, init)

    if V is None:
        cpp.layout_unweighted_convergent(
            X, I, J, t_max, eps, delta, t_maxmax, random_seed
        )
    else:
        cpp.layout_weighted_convergent(
            X, I, J, V, t_max, eps, delta, t_maxmax, random_seed
        )

    return X


def layout_sparse(
    I, J, npivots, V=None, t_max=30, eps=0.01, random_seed=None, init=None
):
    """takes a list of indices I and J
    and returns a n-by-2 matrix of positions X with minimized stress
    using the sparse approximation of Ortmann et al. (2017)"""

    # seed random state
    random_seed = _check_random_seed(random_seed)

    X = random_init(I, J, random_seed, init)

    if npivots > X.shape[0]:
        raise ValueError("number of pivots exceeds number of vertices")

    if V is None:
        cpp.layout_sparse_unweighted(X, I, J, npivots, t_max, eps, random_seed)
    else:
        cpp.layout_sparse_weighted(X, I, J, V, npivots, t_max, eps, random_seed)

    return X


def mds_direct(n, d, w=None, etas=None, num_dimensions=2, random_seed=None, init=None):
    """takes nC2 vectors d (distance) and w (weight) with a vector of step sizes eta
    and returns a n-by-2 matrix of positions X"""

    nC2 = int((n * (n - 1)) / 2)
    if w is None:
        w = np.ones(nC2)  # standard metric MDS

    if len(d) != nC2 or len(w) != nC2:
        raise ValueError(
            "d and/or w are not correct length condensed distance matrices"
        )

    if etas is None:
        etas = default_schedule(w)

    # seed random state
    random_seed = _check_random_seed(random_seed)

    X = _random_init(n, random_seed, init, num_dimensions)

    # do mds
    cpp.mds_direct(X, d, w, etas, random_seed)
    return X


def default_schedule(w, t_max=30, eps=0.01):
    eta_max = 1 / min(w)
    eta_min = eps / max(w)
    lambd = np.log(eta_max / eta_min) / (t_max - 1)
    etas = eta_max * np.exp(-lambd * np.arange(t_max))
    return etas


### no c++ bindings for functions below ###


def _check_random_seed(random_seed=None):
    if random_seed is None:
        random_seed = np.random.randint(65536)
    return random_seed


def _random_init(n, random_seed, init=None, num_dimensions=2):
    # initialize positions
    if init is not None:
        if len(init) != n:
            raise ValueError("initial layout length does not match number of points")
        X = np.ascontiguousarray(init)
        if X.shape[1] != num_dimensions:
            raise ValueError(
                "initial layout dimensionality does not match num_dimensions"
            )
    else:
        np.random.seed(random_seed)
        if num_dimensions == 2:
            X = np.random.rand(n, 2)
        elif num_dimensions == 3:
            X = np.random.rand(n, 3)
        else:
            raise ValueError("only 2D and 3D layouts are supported")
    return X


def random_init(I, J, random_seed, init=None):
    if len(I) == 0 or len(I) != len(J):
        raise ValueError("length of edge indices I and J not equal or zero")

    n = max(max(I), max(J)) + 1

    return _random_init(n, random_seed, init)


def draw_svg(
    X,
    I,
    J,
    filepath=None,
    noderadius=0.2,
    linkwidth=0.05,
    width=1000,
    border=50,
    nodeopacity=1,
    linkopacity=1,
):
    """Takes a n-by-2 matrix of positions X and index pairs I and J
    and writes the equivalent picture in svg format.
    The drawing will be expanded into a width*width square
    Note that the parameters are in svg pixel units.
    The style at the top of the output svg may also be edited as necessary.
    The svg is returned as a string if filename is empty."""

    n = len(X)
    m = len(I)

    X_min = [min(X[:, 0]), min(X[:, 1])]
    X_max = [max(X[:, 0]), max(X[:, 1])]
    range_max = max(X_max[0] - X_min[0], X_max[1] - X_min[1])  # taller or wider
    range_max += 2 * noderadius  # guarantee no nodes are cut off at the edges
    scale = (width - 2 * border) / range_max

    X_svg = np.empty((n, 2))
    for i in range(n):
        X_svg[i] = (X[i] - X_min) * scale
        X_svg[i] += [border + scale * noderadius, border + scale * noderadius]

    svg_list = []
    svg_list.append(
        '<svg width="{:.0f}" height="{:.0f}" xmlns="http://www.w3.org/2000/svg">'.format(
            width, width
        )
    )
    svg_list.append('<style type="text/css">')
    svg_list.append(
        "line{{stroke:black;stroke-width:{:.3f};stroke-opacity:{:.3f};stroke-linecap:round;}}".format(
            scale * linkwidth, linkopacity
        )
    )
    svg_list.append(
        "circle{{r:{};fill:black;fill-opacity:{:.3f}}}".format(
            scale * noderadius, nodeopacity
        )
    )
    svg_list.append("</style>")

    # draw links
    for ij in range(m):
        i = I[ij]
        j = J[ij]
        X_i = X_svg[i]
        X_j = X_svg[j]
        svg_list.append(
            '<line x1="{:.1f}" x2="{:.1f}" y1="{:.1f}" y2="{:.1f}"/>'.format(
                X_i[0], X_j[0], X_i[1], X_j[1]
            )
        )

    # draw nodes
    if noderadius > 0:
        for i in range(n):
            svg_list.append(
                '<circle cx="{:.1f}" cy="{:.1f}"/>'.format(X_svg[i][0], X_svg[i][1])
            )

    svg_list.append("</svg>")

    if filepath is None:
        return "\n".join(svg_list)
    else:
        f = open(filepath, "w")
        f.write("\n".join(svg_list))
        f.close()


import warnings

try:
    ModuleNotFoundError
except NameError:
    # doesn't exist in py2.7
    ModuleNotFoundError = ImportError


def draw_png(
    X,
    I,
    J,
    filepath,
    noderadius=0.2,
    linkwidth=0.05,
    width=1000,
    border=50,
    nodeopacity=1,
    linkopacity=1,
    backend=None,
):
    """Takes a n-by-2 matrix of positions X and index pairs I and J
    and draws it to a .png file with the same format as draw_svg(),
    but using the PyCairo library instead of text svg.
    If PyCairo is not available, defaults to a matplotlib backend instead."""

    if backend is None:
        try:
            import cairo
        except ModuleNotFoundError:
            warnings.warn(
                "pycairo is not installed. Using `backend='matplotlib'` instead.",
                UserWarning,
            )
            try:
                import matplotlib
            except:
                warnings.warn(
                    "matplotlib is not installed either. Either install a valid backend from `['pycairo', 'matplotlib']` or try `draw_svg()` instead."
                )
            else:
                _draw_png_matplotlib(
                    X,
                    I,
                    J,
                    filepath,
                    noderadius,
                    linkwidth,
                    width,
                    border,
                    nodeopacity,
                    linkopacity,
                )
        else:
            _draw_png_cairo(
                X,
                I,
                J,
                filepath,
                noderadius,
                linkwidth,
                width,
                border,
                nodeopacity,
                linkopacity,
            )
    elif backend == "pycairo":
        try:
            import cairo
        except ModuleNotFoundError:
            warnings.warn(
                "Backend .png renderer `pycairo` is not installed. Either install `pycairo` or try `draw_svg()` instead."
            )
        else:
            _draw_png_cairo(
                X,
                I,
                J,
                filepath,
                noderadius,
                linkwidth,
                width,
                border,
                nodeopacity,
                linkopacity,
            )
    elif backend == "matplotlib":
        try:
            import matplotlib
        except ModuleNotFoundError:
            warnings.warn(
                "Backend .png renderer `matplotlib` is not installed. Either install `matplotlib` or try `draw_svg()` instead."
            )
        else:
            _draw_png_matplotlib(
                X,
                I,
                J,
                filepath,
                noderadius,
                linkwidth,
                width,
                border,
                nodeopacity,
                linkopacity,
            )
    else:
        raise ValueError(
            "Expected ['pycairo', 'matplotlib'] in `backend`. Got {}".format(backend)
        )


def _draw_png_cairo(
    X,
    I,
    J,
    filepath,
    noderadius=0.2,
    linkwidth=0.05,
    width=1000,
    border=50,
    nodeopacity=1,
    linkopacity=1,
):
    n = len(X)
    m = len(I)

    X_min = [min(X[:, 0]), min(X[:, 1])]
    X_max = [max(X[:, 0]), max(X[:, 1])]
    range_max = max(X_max[0] - X_min[0], X_max[1] - X_min[1])  # taller or wider
    range_max += 2 * noderadius  # guarantee no nodes are cut off at the edges
    scale = (width - 2 * border) / range_max

    X_png = np.empty((n, 2))
    for i in range(n):
        X_png[i] = (X[i] - X_min) * scale
        X_png[i] += [border + scale * noderadius, border + scale * noderadius]

    # use cairo to draw
    import cairo

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, width)
    ctx = cairo.Context(surface)

    ctx.set_line_width(linkwidth * scale)
    ctx.set_line_cap(cairo.LineCap.ROUND)
    ctx.set_source_rgba(0, 0, 0, linkopacity)

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
        ctx.set_source_rgba(0, 0, 0, nodeopacity)
        radius = noderadius * scale
        for i in range(n):
            ctx.arc(X_png[i][0], X_png[i][1], radius, 0, 7)
            ctx.fill()

    surface.write_to_png(filepath)


# def _draw_png_matplotlib(X, I, J, filepath, noderadius=7, linkwidth=.5, border=1, dpi=None, nodeopacity=1, linkopacity=1):
def _draw_png_matplotlib(
    X,
    I,
    J,
    filepath,
    noderadius=0.2,
    linkwidth=0.05,
    width=1000,
    border=50,
    nodeopacity=1,
    linkopacity=1,
):
    import matplotlib.pyplot as plt
    import matplotlib.collections as mc

    fig = plt.figure(figsize=(width, width))
    ax = plt.axes()
    ax.set_xlim(min(X[:, 0]), max(X[:, 0]))
    ax.set_ylim(min(X[:, 1]), max(X[:, 1]))
    ax.axis("off")
    ax.set_aspect("equal", "box")

    # convert input data widths to display coordinates
    x_display_min, y_display_min = ax.transData.inverted().transform((0, 0))
    noderadius_disp, linkwidth_disp = ax.transData.transform(
        (x_display_min + noderadius, y_display_min + linkwidth)
    )

    links = zip((X[i] for i in I), (X[j] for j in J))
    lc = mc.LineCollection(
        links, linewidths=linkwidth_disp, colors=(0, 0, 0, linkopacity)
    )
    ax.add_collection(lc)

    cc = mc.CircleCollection(
        np.full(len(X), noderadius_disp * noderadius_disp),
        offsets=X,
        transOffset=ax.transData,
        linewidths=0,
        facecolors=(0, 0, 0, nodeopacity),
    )
    ax.add_collection(cc)

    border_data, _ = ax.transLimits.transform((border / width, 0))
    ax.set_xlim(
        min(X[:, 0]) - noderadius - border_data, max(X[:, 0]) + noderadius + border_data
    )
    ax.set_ylim(
        min(X[:, 1]) - noderadius - border_data, max(X[:, 1]) + noderadius + border_data
    )

    fig.savefig(
        filepath,
        dpi=1,
        bbox_inches="tight",
        format="png",
        transparent=True,
        origin="upper",
    )
