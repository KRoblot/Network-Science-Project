import igraph as ig 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap, BoundaryNorm
import math as m
import random
from calculate import initialize_J, initialize_S, initialize_Hs_Hp, calculate_Hs_Hp, calculate_H

nb_nodes = 20
n_iterations = 200
n_calculations = 10

rho0 = np.linspace(0.05, 0.5, 100)
alpha = np.linspace(0.1, 1, 100)
RHO0, ALPHA = np.meshgrid(rho0, alpha)


proba = np.zeros((len(rho0), len(alpha)))

def probability(rho0, alpha):
    
    rho_inf = 0
    
    for i in range(n_calculations):
        
        n_same = 0
        
        proportion_ones = np.random.uniform(0.5, 1)
        
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
                
                if p < 0.5:
                    n_same += 1
            
            elif H < H_old or (H == H_old and p >= 0.5):
                n_same = 0
                
#        print(f"Iteration {n+1}, H: {H} (Reverted), ΔH: {delta_H}")

            # Check for global or local minimum
            if H == -1 or n == n_iterations - 1 or n_same == 250:
                break
#            elif n == n_iterations - 1:
#                break
#            elif n_same == 250:
#                print("Jammed state")
#                break
    
        rho_inf += np.count_nonzero(S == -1) / nb_nodes
            
    return rho_inf / n_calculations
                


for i in range(len(rho0)):
    for j in range(len(alpha)):
        print(f"Calculating probability for rho0 = {rho0[i]}, alpha = {alpha[j]}")  
        proba[i, j] = probability(rho0[i], alpha[j]) - rho0[i]
        proba[i, j] = m.ceil(proba[i, j] / 0.05) * 0.05
        print(proba[i, j])
        
# Transpose the proba array
proba = proba.T

# Sort rho0 and alpha arrays
rho0 = np.sort(rho0)
alpha = np.sort(alpha)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Create meshgrid for rho0 and alpha
RHO0, ALPHA = np.meshgrid(rho0, alpha)

# Define custom colormap
colors = ['#000080', '#0000CD', '#0000FF', '#ADD8E6', '#50C878', '#90EE90', '#FFFF00', '#FFA500', '#FF0000', '#8B0000']
cmap = ListedColormap(colors)
bounds = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
norm = BoundaryNorm(bounds, cmap.N)

# Plot the surface
surf = ax.plot_surface(RHO0, ALPHA, proba, cmap=cmap, norm=norm)

fig.colorbar(surf, boundaries=bounds, ticks=bounds)

ax.set_xlabel('rho0')
ax.set_ylabel('alpha')
ax.set_zlabel('rho_inf - rho0')

# Inversion de l’axe x (rho0)
ax.invert_xaxis()

plt.show()