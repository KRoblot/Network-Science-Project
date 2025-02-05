import math as m
import random
import numpy as np

def initialize_J(nb_nodes, proportion_ones):
    
    """Initialize the matrix J by creating a numpy array full of -1 and changing randomly a proportion of the edges to 1

    Args:
        nb_nodes (int): number of nodes
        proportion_ones (float): proportion of friendly edges

    Returns:
        np.array: matrix J
    """
    
    J = np.full((nb_nodes, nb_nodes), -1)                                        # initialize the matrix J by creating a numpy array full of -1 

    num_ones = m.ceil(proportion_ones * nb_nodes * (nb_nodes - 1) / 2)           # calculate the number of friendly edges
    indices = [(i, j) for i in range(nb_nodes) for j in range(i + 1, nb_nodes)]  
    selected_indices = random.sample(indices, num_ones)                          # select randomly num_ones edges to be friendly
    
    for i, j in selected_indices:                                                # change the matrix J
        J[i, j] = 1
        J[j, i] = 1
    
    return J

def initialize_S(nb_nodes, rho0):
    
    """Initialize the matrix J by creating a numpy array full of -1 and changing randomly a proportion of the edges to 1

    Args:
        nb_nodes (int): number of nodes
        rho0 (float): proportion of infected nodes

    Returns:
        np.array: matrix S
    """
    
    num_neg_ones = m.ceil(nb_nodes * rho0)                                      # calculate the initial number of infected nodes
    num_ones = nb_nodes - num_neg_ones                                          # calculate the initial number of susceptible nodes
    
    S = np.array([-1] * num_neg_ones + [1] * num_ones)                          # initialize the matrix S 
    np.random.shuffle(S)                                                        # choose randomly the infected nodes
    
    return S

def initialize_Hs_Hp(triangles, J, S):     
    
    """Initialize the dictionaries Hs and Hp

    Args:
        triangles (list): list of triangles
        J (np.array): matrix J
        S (np.array): matrix S

    Returns:
        (dictionary, dictionary): Hs and Hp
    """
    
    #  initialize two dictionaries giving Hs and Hp for each triple
    Hs, Hp = {}, {}
    
    for triangle in triangles:
        i, j, k = triangle
        Hs[tuple(triangle)] = - J[i, j] * J[j, k] * J[k, i]                     # calculate energy Hs for the triple
        Hp[tuple(triangle)] = -(                                                # calculate energy Hp for the triple
            (S[i] * S[j])**((3 - J[i, j])/2) * (S[j] * S[k])**((3 - J[j, k])/2) +
            (S[j] * S[k])**((3 - J[j, k])/2) * (S[k] * S[i])**((3 - J[k, i])/2) +
            (S[k] * S[i])**((3 - J[k, i])/2) * (S[i] * S[j])**((3 - J[j, i])/2) - 1
        ) / 2
    
    return Hs, Hp

def calculate_Hs_Hp(i, j, Hs, Hp, S, J, triangles):
    
    """Calculate the energy Hs and Hp for the triples containing the nodes i and j

    Args:
        i (int): node i
        j (int): node j
        Hs (dictionary): dictionary of Hs
        Hp (dictionary): dictionary of Hp
        S (np.array): matrix S
        J (np.array): matrix J
        triangles (list): list of triangles

    Returns:
        (dictionary, dictionary): updated Hs and Hp
    """
    
    for triangle in triangles:
        if i in triangle and j in triangle:                                     # consider each triple with i and j 
            i_temp, j_temp, k = triangle
            Hs[tuple(triangle)] = - J[i_temp, j_temp] * J[j_temp, k] * J[k, i_temp]         # change Hs and Hp for these triples
            Hp[tuple(triangle)] = -(
                (S[i_temp] * S[j_temp])**((3 - J[i_temp, j_temp])/2) * (S[j_temp] * S[k])**((3 - J[j_temp, k])/2) +
                (S[j_temp] * S[k])**((3 - J[j_temp, k])/2) * (S[k] * S[i_temp])**((3 - J[k, i_temp])/2) +
                (S[k] * S[i_temp])**((3 - J[k, i_temp])/2) * (S[i_temp] * S[j_temp])**((3 - J[j_temp, i_temp])/2) - 1
            ) / 2
    return Hs, Hp
            
def calculate_H(nb_nodes, triangles, Hp, Hs):    
     
    """Calculate the energy H

    Args:
        nb_nodes (int): number of nodes
        triangles (list): list of triangles
        Hp (dictionary): dictionary of Hp
        Hs (dictionary): dictionary of Hs

    Returns:
        float: energy H
    """
    
    return -(1 / m.comb(nb_nodes, 3)) * sum(
        Hp[triangle] * (Hs[triangle] - 1) - (Hs[triangle] + 1) for triangle in triangles
    ) / 2
