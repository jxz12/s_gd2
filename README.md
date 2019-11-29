# s(gd)²
Stochastic Gradient Descent for Graph Drawing [<https://arxiv.org/abs/1710.04626>]

![image](comparison.gif)

We recommend using the available python package, implemented in C++ using SWIG to generate bindings. It may be installed using
```
$ pip install s_gd2
```
and an example use case looks like
```python
import s_gd2
I = [0,1,2,3,4]
J = [1,2,3,4,0]
X = s_gd2.layout(I, J)
s_gd2.draw_svg(X, I, J, 'C5.svg')
```

Other useful functions include the following:
```python
layout_convergent(I, J, V=None, t_max=30, eps=.1, delta=.03, t_maxmax=200)
```
This function implements the convergent schedule from the paper, Section 2.1.2. If `V` is provided, the graph is treated as weighted (same as the non-convergent case). It produces layouts that are slightly tidier than `layout()`, but requires more iterations. `t_maxmax` is a maximum iteration regardless of delta.
```python
layout_sparse(I, J, npivots, V=None, t_max=30, eps=.01)
```
This is an implementation of the sparse stress approximation of Ortmann et al. (2017), described in Section 4.3 in the paper. It allows the algorithm to draw large graphs of up to millions of vertices, by strategically cutting terms from the loss function. A standard value for `npivots` is 200.
```python
draw_png(X, I, J, filepath, noderadius=.2, linkwidth=.05, framewidth=1000, border=50, nodeopacity=1, linkopacity=1)
```
This function draws the same image as `draw_svg()` and has the same parameters, but uses the `pycairo` library to draw it directly onto a .png file. This function is especially useful when drawing large graphs, as the equivalent .svg files can become too large to render in common web browsers.
```python
mds_direct(n, d, w, etas=None)
```
This function directly optimises the stress function (Equation (1) in the paper) given n vertices, and condensed distance matrices `d` and `w` (see `scipy.spatial.distance.squareform`). The step sizes are given as input `etas`, which defaults to the same schedule as `layout()` if not provided.

## Code used for the paper
The (old) code used for timing experiments in the paper is in C#, run as a command line application that takes paths as command line arguments: input .txt file, output stress trajectory, output .svg layout.

The graph data used comes from <https://sparse.tamu.edu>, downloaded in .mat format. A matlab script is provided in the same folder to convert it into the unweighted .txt format that the C# code uses. If a graph has disconnected subgraphs, it will select the largest connected component.

A Jupyter notebook is also provided, which is rather slow, but serves as a good introduction to the algorithm. It can be rendered directly in GitHub, showing the output for the graph `qh882`, as well as `commanche_dual` in an implementation of the Sparse Stress Model of Ortmann et al. (2017).
