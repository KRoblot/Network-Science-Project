import igraph as ig 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import random
from calculate import initialize_J, initialize_S, initialize_Hs_Hp, calculate_Hs_Hp, calculate_H

nb_nodes = 15
rho0 = 0.375
alpha = 0.5
n_iterations = 30000
proportion_ones = 0.5

J = initialize_J(nb_nodes, proportion_ones)

S = initialize_S(nb_nodes, rho0)

g = ig.Graph.Erdos_Renyi(n=nb_nodes, p=1)

triangles = g.cliques(min=3, max=3)

Hs, Hp = initialize_Hs_Hp(triangles, J, S)

H = calculate_H(nb_nodes, triangles, Hp, Hs)

# Perform the simulated annealing
for n in range(n_iterations):
    # Save the current state
    J_old = J.copy()
    S_old = S.copy()
    H_old = H
    
    probability = np.random.uniform(0,1.0)

    # Randomly choose a pair connection
    i, j = random.sample(range(nb_nodes), 2) 
    if S[i] == S[j]:
            # Rule 1: Both nodes are susceptible or infected, change the sign of the edge
        J[i, j] = -J[i, j]
        J[j, i] = J[i, j]
    elif S[i] != S[j] and J[i, j] == 1:
            # Rule 2 and 3: One node is susceptible and the other is infected, and the edge is friendly
        if probability < alpha:
                # Disease spreads
            S[i] = S[j] = -1
            # Change the sign of the edge
        J[i, j] = -J[i, j]
        J[j, i] = J[i, j]
        # Rule 4: If one node is susceptible and the other is infected and the edge is unfriendly, nothing happens

    Hs, Hp = calculate_Hs_Hp(i, j, Hs, Hp, S, J, triangles)
    H = calculate_H(nb_nodes, triangles, Hp, Hs)
    
    # Calculate the difference in H
    delta_H = H - H_old
    
    probability = np.random.uniform(0,1.0)

    # Apply the global constraint
    if H < H_old:
        # Accept the new configuration
        print(f"Iteration {n+1}, H: {H} (Accepted), ΔH: {delta_H}")
    elif H == H_old and probability < 0.5:
        # Accept the new configuration with probability 1/2
        print(f"Iteration {n+1}, H: {H} (Accepted with probability 1/2), ΔH: {delta_H}")
    else:
        # Revert to the old configuration
        J = J_old
        S = S_old
        H = H_old
        print(f"Iteration {n+1}, H: {H} (Reverted), ΔH: {delta_H}")

    # Check for global or local minimum
    if H == -1:
        print("Global minimum reached (H = -1)")
        break
    elif n == n_iterations - 1:
        print("Local minimum reached (H > -1)")
        
node_colors = ["green" if S[i] == 1 else "red" for i in range(nb_nodes)]
g.vs["color"] = node_colors

edge_colors = ["green" if J[edge.source, edge.target] == 1 else "red" for edge in g.es]
g.es["color"] = edge_colors

layout = g.layout("kk")

ig.plot(g, layout=layout, bbox=(1000, 1000), target = "graph.png")