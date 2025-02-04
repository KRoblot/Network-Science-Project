import igraph as ig 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import random
from calculate import initialize_J, initialize_S, initialize_Hs_Hp, calculate_Hs_Hp, calculate_H

alpha = 0.5
rho0 = 0.35
nb_nodes = 20

proportion_ones = np.linspace(0, 1, 100)

n_iterations = 10000
n_calculations = 20

proportion = np.zeros_like(proportion_ones)

def calculate_proportion(proportion_ones):
    
    rho_inf = 0
    
    for i in range(n_calculations):
        
        n_same = 0
        
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

for i in range(len(proportion_ones)):
    print(f"Calculating proportion for proportion_ones = {proportion_ones[i]}")  
    proportion[i] = calculate_proportion(proportion_ones[i]) - rho0
    print(proportion[i])
        
plt.plot(proportion_ones, proportion)
plt.xlabel("Initial proportion of infected nodes (%)")
plt.ylabel("Fraction of nodes that become infected at the steady state (%)")
plt.title("Fraction of nodes that become infected at the steady state as a function of the initial proportion of infected nodes")
plt.show()