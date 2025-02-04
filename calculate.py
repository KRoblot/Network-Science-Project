import math as m
import random
import numpy as np

def initialize_J(nb_nodes, proportion_ones):
    J = np.full((nb_nodes, nb_nodes), -1)

    num_ones = m.ceil(proportion_ones * nb_nodes * (nb_nodes - 1) / 2)
    indices = [(i, j) for i in range(nb_nodes) for j in range(i + 1, nb_nodes)]
    selected_indices = random.sample(indices, num_ones)
    
    for i, j in selected_indices:
        J[i, j] = 1
        J[j, i] = 1
    
    return J

def initialize_S(nb_nodes, rho0):
    num_neg_ones = m.ceil(nb_nodes * rho0)
    num_ones = nb_nodes - num_neg_ones
    
    S = np.array([-1] * num_neg_ones + [1] * num_ones)
    np.random.shuffle(S)
    
    return S

def initialize_Hs_Hp(triangles, J, S):
    Hs, Hp = {}, {}
    
    for triangle in triangles:
        i, j, k = triangle
        Hs[tuple(triangle)] = - J[i, j] * J[j, k] * J[k, i]
        Hp[tuple(triangle)] = -(
            (S[i] * S[j])**((3 - J[i, j])/2) * (S[j] * S[k])**((3 - J[j, k])/2) +
            (S[j] * S[k])**((3 - J[j, k])/2) * (S[k] * S[i])**((3 - J[k, i])/2) +
            (S[k] * S[i])**((3 - J[k, i])/2) * (S[i] * S[j])**((3 - J[j, i])/2) - 1
        ) / 2
    
    return Hs, Hp

def calculate_Hs_Hp(i, j, Hs, Hp, S, J, triangles):
    for triangle in triangles:
        if i in triangle and j in triangle:
            i_temp, j_temp, k = triangle
            Hs[tuple(triangle)] = - J[i_temp, j_temp] * J[j_temp, k] * J[k, i_temp]
            Hp[tuple(triangle)] = -(
                (S[i_temp] * S[j_temp])**((3 - J[i_temp, j_temp])/2) * (S[j_temp] * S[k])**((3 - J[j_temp, k])/2) +
                (S[j_temp] * S[k])**((3 - J[j_temp, k])/2) * (S[k] * S[i_temp])**((3 - J[k, i_temp])/2) +
                (S[k] * S[i_temp])**((3 - J[k, i_temp])/2) * (S[i_temp] * S[j_temp])**((3 - J[j_temp, i_temp])/2) - 1
            ) / 2
    return Hs, Hp
            
def calculate_H(nb_nodes, triangles, Hp, Hs):
    return -(1 / m.comb(nb_nodes, 3)) * sum(
        Hp[triangle] * (Hs[triangle] - 1) - (Hs[triangle] + 1) for triangle in triangles
    ) / 2