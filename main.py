import igraph as ig 
import numpy as np
import pandas as pd

nb_nodes = 1000
rho0 = 0

J = np.full((nb_nodes, nb_nodes), -1)
array = np.random.choice([-1, 1], size=nb_nodes, p=[rho0, 1-rho0])