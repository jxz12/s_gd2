# s(gd)²
Stochastic Gradient Descent for Graph Drawing [<https://arxiv.org/abs/1710.04626>]

![image](comparison.gif)

We recommend using the available python package, implemented in C++ using SWIG to generate bindings. It may be installed using
```
$ pip install s_gd2
```
and an example use case looks like
```python
import networkx as nx
G = nx.balanced_tree(2,9)
S = nx.to_scipy_sparse_matrix(G)

import s_gd2
X = s_gd2.layout_scipy(S)
nx.draw(G, pos=X, node_size=0)

import matplotlib.pyplot as plt
plt.axis('equal')
plt.show()
```

Other useful functions include the following:
```python
layout_scipy_y_constrained(S, Y, weighted=False, t_max=20, eps=.1)
```
is used to constrain the y-axis positions to the values in `Y`. This is useful for generating Sugiyama-style layouts, or for embedding metadata within the layout itself e.g. the trophic level of species in a food web.
```python
layout_scipy_focus(S, f, weighted=False, t_max=20, eps=.01)
```
can be used to focus on a vertex with index `f`, as in Section 4.1 of the paper. Usually requires larger t_max and smaller eps in order to guarantee a tidy layout.
```python
concentric_circles(n_circ, radius=1, n_seg=100, offset=[0,0])
```
returns sets of coordinates to draw concentric circles that can be used with matplotlib.pyplot.plt. It is in the form of a 3D numpy matrix, where the first dimension is which circle, the second is x/y axis, and the third is the coordinates of the circle.

## Code used for the paper
The code used for timing experiments in the paper is in C#, run as a command line application that takes paths as command line arguments: input .txt file, output stress trajectory, output .svg layout.

The graph data used comes from <https://sparse.tamu.edu>, downloaded in .mat format. A matlab script is provided in the same folder to convert it into the unweighted .txt format that the C# code uses. If a graph has disconnected subgraphs, it will select the largest connected component.

A Jupyter notebook is also provided, which is rather slow, but serves as a good introduction to the algorithm. It can even be rendered directly in GitHub, showing the output for the graph `qh882`, as well as `commanche_dual` in an implementation of the Sparse Stress Model of Ortmann et al. [<https://arxiv.org/abs/1608.08909>].
