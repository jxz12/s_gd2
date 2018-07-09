# get graph in sparse matrix format
import networkx as nx

G = nx.karate_club_graph()
S = nx.to_scipy_sparse_matrix(G)

# calculate distances and weights
import scipy.sparse.csgraph as csg
import scipy.spatial.distance as dist

d = csg.shortest_path(S, directed=False, unweighted=True)
d = dist.squareform(d)
w = 1/d**2

# prepare annealing schedule
import numpy as np
t_max = 15
eta_max = 1/min(w)
eta_min = 0.1/max(w)
lambd = np.log(eta_min/eta_max) / (t_max-1)

eta = np.arange(t_max)
eta = eta_max * np.exp(lambd*eta)

# initialize positions
n = G.number_of_nodes()
X = np.random.rand(n, 2)
