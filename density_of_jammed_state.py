import igraph as ig 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import random
from calculate import initialize_J, initialize_S, initialize_Hs_Hp, calculate_Hs_Hp, calculate_H

nb_nodes = 5
n_iterations = 100
n_calculations = 100
proportion_ones = 0.75
epsilon = 0.01


rho0 = np.linspace(0, 1, 15)
alpha = np.linspace(0, 1, 15)
RHO0, ALPHA = np.meshgrid(rho0, alpha)

proba = np.zeros((len(rho0), len(alpha)))

def count_connected_components(J, nb_nodes):
    """Compte le nombre de composantes connexes en ne gardant que les arêtes J(i,j) = +1"""
    g = ig.Graph(directed=False)
    g.add_vertices(nb_nodes)
    for i in range(nb_nodes):
        for j in range(i + 1, nb_nodes):
            if J[i, j] == 1:
                g.add_edge(i, j)
    return len(g.components().subgraphs())

def probability(rho0, alpha):
    
    n_jammed = 0

    for i in range(n_calculations):
        
        J = initialize_J(nb_nodes, proportion_ones)
        S = initialize_S(nb_nodes, rho0)
        g = ig.Graph.Erdos_Renyi(n=nb_nodes, p=1)
        triangles = g.cliques(min=3, max=3)
        Hs, Hp = initialize_Hs_Hp(triangles, J, S)
        H = calculate_H(nb_nodes, triangles, Hp, Hs)

        # Perform the simulated annealing
        for n in range(n_iterations):
            J_old = J.copy()
            S_old = S.copy()
            H_old = H

            i, j = random.sample(range(nb_nodes), 2)
            if g.are_adjacent(i, j):
                if S[i] == S[j]:
                    J[i, j] = -J[i, j]
                    J[j, i] = J[i, j]
                elif S[i] != S[j] and J[i, j] == 1:
                    if random.random() < alpha:
                        S[i] = S[j] = -1
                    J[i, j] = -J[i, j]
                    J[j, i] = J[i, j]

            Hs, Hp = calculate_Hs_Hp(i, j, Hs, Hp, S, J, triangles)
            H = calculate_H(nb_nodes, triangles, Hp, Hs)
            
            delta_H = H - H_old
            p = random.random()
            
            if H < H_old:
                pass
            elif (abs(H-H_old)<=epsilon):

                if (p < 0.5):
                    pass
            else:
                J = J_old
                S = S_old
                H = H_old

            # Nouvelle condition d'arrêt : vérification des composantes connexes
            if count_connected_components(J, nb_nodes) >= 2:
                n_jammed += 1
                break
            if H == -1:
                break

    
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
