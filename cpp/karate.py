print('starting wrangle')
# get graph in sparse matrix format
import networkx as nx

G = nx.karate_club_graph()
# G = nx.cycle_graph(10)
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

print('starting sgd')

# optimize
import s_gd2
s_gd2.sgd(X, d, w, eta)

print('sgd done')

# draw
import matplotlib.pyplot as plt
nx.draw(G, pos=X)
plt.axis('equal')
plt.savefig('karate.png')

