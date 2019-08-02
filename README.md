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
s_gd2.draw(X, I, J, 'C5.svg')
```

Other useful functions include the following:
```python
layout_convergent(I, J, V=None, t_max=30, eps=.1, delta=.03, t_maxmax=200)
```
This function implements the convergent schedule from the paper, Section 2.1.2. If `V` is provided, the graph is treated as weighted (same as the non-convergent case). It produces layouts that are slightly tidier than `layout()`, but requires more iterations. `t_maxmax` is a maximum iteration regardless of delta.
```python
mds_direct(n, d, w, etas)
```
This function directly optimises the stress function (Equation (1) in the paper) given n vertices, and condensed distance matrices `d` and `w` (see `scipy.spatial.distance.squareform`). The step sizes are given as input `etas`.

TODO: a fast implementation of the sparse model, from Section 4.3 of the paper, is in the works.

## Code used for the paper
The code used for timing experiments in the paper is in C#, run as a command line application that takes paths as command line arguments: input .txt file, output stress trajectory, output .svg layout.

The graph data used comes from <https://sparse.tamu.edu>, downloaded in .mat format. A matlab script is provided in the same folder to convert it into the unweighted .txt format that the C# code uses. If a graph has disconnected subgraphs, it will select the largest connected component.

A Jupyter notebook is also provided, which is rather slow, but serves as a good introduction to the algorithm. It can even be rendered directly in GitHub, showing the output for the graph `qh882`, as well as `commanche_dual` in an implementation of the Sparse Stress Model of Ortmann et al. [<https://arxiv.org/abs/1608.08909>].
