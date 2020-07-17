[![Build Status](https://travis-ci.com/jxz12/s_gd2.svg?branch=master)](https://travis-ci.com/jxz12/s_gd2)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/5h4mp93ku0ft62ha/branch/master?svg=true)](https://ci.appveyor.com/project/jxz12/s-gd2/branch/master)


# (sgd)²
Stress-based Graph Drawing by Stochastic Gradient Descent [[arXiv:1710.04626]](https://arxiv.org/abs/1710.04626). A video of the corresponding talk (given at IEEE VIS 2019) can be viewed at [vimeo.com/373015168](https://vimeo.com/373015168).

![image](https://media.giphy.com/media/JoaboGdTq1sXNnnIND/giphy.gif)

We recommend using the available python package, implemented in C++ using SWIG to generate bindings. It may be installed through conda using
```
conda install -c conda-forge s_gd2
```
or through pip with
```
pip install s_gd2
```
and an example use case looks like
```python
import s_gd2
I = [0,1,2,3,4]
J = [1,2,3,4,0]
X = s_gd2.layout(I, J)
s_gd2.draw_svg(X, I, J, 'C5.svg')
```

Useful functions include the following:
```python
layout(I, J, V=None, t_max=30, eps=.01, random_seed=None, init=None)
```
takes two lists `I` and `J` as edge indices for a graph, and lays it out using stochastic gradient descent. If `V` is provided, the graph is treated as weighted. `t_max` and `eps` are parameters used to determine the running time of the algorithm, as in Section 2.1.1 of the paper. `random_seed` is an optional integer used to seed random number generation to produce the same layouts in multiple runs, and `init` may be used if an initial layout is provided.
```python
layout_convergent(I, J, V=None, t_max=30, eps=.01, delta=.03, t_maxmax=200, random_seed=None, init=None)
```
implements the convergent schedule from the paper, Section 2.1.2. It produces layouts that are slightly tidier than `layout()`, but requires more iterations. `t_maxmax` is a maximum iteration regardless of the stopping threshold `delta`.
```python
layout_sparse(I, J, npivots, V=None, t_max=30, eps=.01, random_seed=None, init=None)
```
is an implementation of the sparse stress approximation of Ortmann et al. (2017), described in Section 4.3 in the paper. It allows the algorithm to draw large graphs of up to millions of vertices, by strategically cutting terms from the loss function. A standard value for `npivots` is 200.
```python
mds_direct(n, d, w=None, etas=None, num_dimensions=2, random_seed=None, init=None)
```
directly optimises the stress function (Equation (1) in the paper) given n vertices, and condensed distance matrices `d` and `w` (see `scipy.spatial.distance.squareform`). If `w` is not provided, it is initialised to an array full of 1s. If `etas` is not provided, it defaults to the same schedule as used in `layout()`. `num_dimensions` may also be set to `3` if a 3-dimensional layout is desired.
```python
draw_svg(X, I, J, filepath, noderadius=.2, linkwidth=.05, framewidth=1000, border=50, nodeopacity=1, linkopacity=1)
```
renders a given layout into .svg format.
```python
draw_png(X, I, J, filepath, noderadius=.2, linkwidth=.05, framewidth=1000, border=50, nodeopacity=1, linkopacity=1)
```
draws the same image as `draw_svg()` and has the same parameters, but uses the `pycairo` library to draw it directly onto a .png file. This function is especially useful when drawing large graphs, as the equivalent .svg files can become too large to render in common web browsers.


## Building from source
To build the package from source, the easiest way is through Python `setuptools`. To install from source with precompiled C++ wrapper code:

```shell
git clone https://github.com/jxz12/s_gd2
cd s_gd2/cpp
python setup.py install
```
To recompile C++ wrappers (required if there are any API changes):
```shell
git clone https://github.com/jxz12/s_gd2
cd s_gd2/cpp
swig -python -c++ s_gd2/swig/layout.i
python setup.py install
```

## Code used for the paper
The (old) code used for timing experiments in the paper is in C#, run as a command line application that takes paths as command line arguments: input .txt file, output stress trajectory, output .svg layout.

The graph data used comes from <https://sparse.tamu.edu>, downloaded in .mat format. A matlab script is provided in the same folder to convert it into the unweighted .txt format that the C# code uses. If a graph has disconnected subgraphs, it will select the largest connected component.

A Jupyter notebook is also provided, which is rather slow, but serves as a good introduction to the algorithm. It can be rendered directly in GitHub, showing the output for the graph `qh882`, as well as `commanche_dual` in an implementation of the Sparse Stress Model of Ortmann et al. (2017).
