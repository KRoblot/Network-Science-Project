import igraph as ig 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import random
from calculate import calculate_Hs_Hp, calculate_H

nb_nodes = 5
n_iterations = 100
n_calculations = 100
proportion_ones = 0.1 


rho0 = np.linspace(0, 0.5, 25)
alpha = np.linspace(0, 1, 50)
RHO0, ALPHA = np.meshgrid(rho0, alpha)


proba = np.zeros((len(rho0), len(alpha)))

def probability(rho0, alpha):
    
    n_jammed = 0

    for i in range(n_calculations):
        J = np.full((nb_nodes, nb_nodes), -1)

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

            Hs, Hp = calculate_Hs_Hp(i, j, Hs, Hp, S, J, triangles)
            H = calculate_H(nb_nodes, triangles, Hp, Hs)
            
            # Calculate the difference in H
            delta_H = H - H_old

            p = random.random()
            
            if not H < H_old and not (H == H_old and p < 0.5):
                # Revert to the old configuration
                J = J_old
                S = S_old
                H = H_old
        #        print(f"Iteration {n+1}, H: {H} (Reverted), ΔH: {delta_H}")

            # Check for global or local minimum
            if H == -1:
                break
            elif n == n_iterations - 1:
                n_jammed += 1
    
    return n_jammed / n_calculations
                


for i in range(len(rho0)):
    for j in range(len(alpha)):
        print(f"Calculating probability for rho0 = {rho0[i]}, alpha = {alpha[j]}")  
        proba[i, j] = probability(rho0[i], alpha[j])
        
# Transpose the proba array
proba = proba.T

plt.figure(figsize=(6,4))

# Define the levels and custom colormap
levels = np.linspace(0, 1, 11)
colors = ['darkblue', 'blue', 'deepskyblue', 'cyan', 'green', 'yellow', 'orange', 'red', "firebrick", 'darkred']
cmap = plt.cm.colors.ListedColormap(colors)

contour = plt.contourf(RHO0, ALPHA, proba, levels=levels, cmap=cmap)

cbar = plt.colorbar(contour)

plt.xlabel(r'$\rho_0$', fontsize=14)
plt.ylabel(r'$\alpha$', fontsize=14)

plt.show()