import igraph as ig 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import random
from calculate import calculate_Hs_Hp, calculate_H

nb_nodes = 5
rho0 = 0.375
alpha = 0.2
n_iterations = 3000

J = np.full((nb_nodes, nb_nodes), -1)
    
    # Set a proportion of 1s in J while keeping the diagonal -1
proportion_ones = 0.1  # Adjust this proportion as needed
num_ones = int(proportion_ones * nb_nodes * (nb_nodes - 1) / 2)
indices = [(i, j) for i in range(nb_nodes) for j in range(i + 1, nb_nodes)]
selected_indices = random.sample(indices, num_ones)
for i, j in selected_indices:
    J[i, j] = 1
    J[j, i] = 1

num_neg_ones = int(nb_nodes * rho0)
num_ones = nb_nodes - num_neg_ones
S = np.array([-1] * num_neg_ones + [1] * num_ones)
np.random.shuffle(S)

g = ig.Graph.Erdos_Renyi(n=nb_nodes, p=1)

triangles = g.cliques(min=3, max=3)

Hs = {}
Hp = {}

# Initialize Hs, Hp, and H
for triangle in triangles:
    i, j, k = triangle
    Hs[tuple(triangle)] = - J[i, j] * J[j, k] * J[k, i]
    Hp[tuple(triangle)] = -(
        (S[i] * S[j])**((3 - J[i, j])/2) * (S[j] * S[k])**((3 - J[j, k])/2) +
        (S[j] * S[k])**((3 - J[j, k])/2) * (S[k] * S[i])**((3 - J[k, i])/2) +
        (S[k] * S[i])**((3 - J[k, i])/2) * (S[i] * S[j])**((3 - J[j, i])/2) - 1
    ) / 2

H = calculate_H(nb_nodes, triangles, Hp, Hs)

# Perform the simulated annealing
for n in range(n_iterations):
    # Save the current state
    J_old = J.copy()
    S_old = S.copy()
    H_old = H
    
    probability = np.random.uniform(0,1)

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
    
    probability = np.random.uniform(0,1)

    # Apply the global constraint
    if H < H_old:
        # Accept the new configuration
        print(f"Iteration {n+1}, H: {H} (Accepted), ΔH: {delta_H}")
    elif H == H_old and random.random() < 0.5:
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

# Plot the results
plt.figure(figsize=(10, 6))

for a in alpha:
    plt.plot(np.arange(n_iterations) / 10, rho_t[str(a)] - rho0, label=f'alpha = {a}')

plt.xlabel('Iteration (x10)', fontsize=14)
plt.ylabel(r'$\rho_t - \rho_0$', fontsize=14)
plt.legend()
plt.title('Evolution of Infected Nodes Density Over Time', fontsize=16)
plt.show()
        
node_colors = ["green" if S[i] == 1 else "red" for i in range(nb_nodes)]
g.vs["color"] = node_colors

edge_colors = ["red" if J[edge.source, edge.target] == 1 else "black" for edge in g.es]
g.es["color"] = edge_colors

layout = g.layout("kk")

ig.plot(g, layout=layout, bbox=(1000, 1000), target = "graph.png")