# s(gd)²
Stochastic Gradient Descent for Graph Drawing [<https://arxiv.org/abs/1710.04626>]

![image](comparison.gif)

The main code is in C#, run as a command line application that takes paths as command line arguments: input .txt file, output stress trajectory, output .svg layout.

The graph data used comes from <https://sparse.tamu.edu>, downloaded in .mat format. A matlab script is provided in the same folder to convert it into the unweighted .txt format that the C# code uses. If a graph has disconnected subgraphs, it will select the largest weakly connected component.

A Jupyter notebook is also provided, which is rather slow, but serves as a good introduction to the algorithm. It can even be rendered directly in GitHub, showing the output for the graph `qh882`, as well as `commanche\_dual` in an implementation of the Sparse Stress Model of Ortmann et al. [<https://arxiv.org/abs/1608.08909>].
