import math as m

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