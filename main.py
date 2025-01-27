import igraph as ig 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import random

nb_nodes = 20
rho0 = 0.3
alpha = 0.15
n_iterations = 100
p = 0.1 

J = np.full((nb_nodes, nb_nodes), -1)

random_matrix = np.random.rand(nb_nodes, nb_nodes)
J[random_matrix < p] = 1

np.fill_diagonal(J, -1)

num_neg_ones = int(nb_nodes * rho0)
num_ones = nb_nodes - num_neg_ones
S = np.array([-1] * num_neg_ones + [1] * num_ones)
np.random.shuffle(S)

g = ig.Graph.Erdos_Renyi(n=nb_nodes, p=1)

triangles = g.cliques(min=3, max=3)

Hs = {}
Hp = {}

def calculate_Hs_Hp():
    for triangle in triangles:
        i, j, k = triangle
        Hs[tuple(triangle)] = J[i, j] * J[j, k] * J[k, i]
        Hp[tuple(triangle)] = -(
            (S[i] * S[j])**((3 - J[i, j])/2) * (S[j] * S[k])**((3 - J[j, k])/2) +
            (S[j] * S[k])**((3 - J[j, k])/2) * (S[k] * S[i])**((3 - J[k, i])/2) +
            (S[k] * S[i])**((3 - J[k, i])/2) * (S[i] * S[j])**((3 - J[j, i])/2) - 1
        ) / 2

def calculate_H():
    return -(1 / m.comb(nb_nodes, 3)) * sum(
        Hp[triangle] * (Hs[triangle] - 1) - (Hs[triangle] + 1) for triangle in triangles
    ) / 2

calculate_Hs_Hp()

# Debug prints for initial values
print("Initial Hs:", Hs)
print("Initial Hp:", Hp)

H = calculate_H()
print(f"Initial H: {H}")

for n in range(n_iterations):
    # Save the current state
    J_old = J.copy()
    S_old = S.copy()
    H_old = H

    # Randomly choose a pair connection
    i, j = random.sample(range(nb_nodes), 2)
    if g.are_adjacent(i, j):
        if S[i] == S[j]:
            # Rule 1: Both nodes are susceptible or infected, change the sign of the edge
            J[i, j] = -J[i, j]
            J[j, i] = J[i, j]
        elif S[i] != S[j] and J[i, j] == 1:
            # Rule 2 and 3: One node is susceptible and the other is infected, and the edge is friendly
            if random.random() < alpha:
                # Disease spreads
                S[i] = S[j] = -1
            # Change the sign of the edge
            J[i, j] = -J[i, j]
            J[j, i] = J[i, j]
        # Rule 4: If one node is susceptible and the other is infected and the edge is unfriendly, nothing happens

    calculate_Hs_Hp()
    H = calculate_H()
    
    # Calculate the difference in H
    delta_H = H - H_old

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