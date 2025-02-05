import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import igraph as ig
from calculate import initialize_J, initialize_S, initialize_Hs_Hp, calculate_Hs_Hp, calculate_H
from multiprocessing import Pool, cpu_count

# Définir la fonction pour calculer rho_inf - rho0
def calculate_rho_inf(nb_nodes, S):
    l = []
    for i in range(nb_nodes):
        if S[i] == -1:  # Si le nœud est infecté
            l.append(i)
    rho_inf = len(l) / nb_nodes  # Proportion de nœuds infectés
    return rho_inf

# Fonction de simulation pour une combinaison donnée de rho et alpha
def simulate_combination(params):
    rho_value, alpha_value = params
    nb_nodes = 15
    n_iterations = 9000  # Réduire les itérations à 1000

    # Initialisation des variables pour l'algorithme
    J = initialize_J(nb_nodes, 0.5)
    S = initialize_S(nb_nodes, rho_value)
    g = ig.Graph.Erdos_Renyi(n=nb_nodes, p=1)
    triangles = g.cliques(min=3, max=3)
    Hs, Hp = initialize_Hs_Hp(triangles, J, S)
    H = calculate_H(nb_nodes, triangles, Hp, Hs)

    # Perform the simulated annealing
    for n in range(n_iterations):
        J_old = J.copy()
        S_old = S.copy()
        H_old = H
        probability = np.random.uniform(0, 1.0)
        i, j = random.sample(range(nb_nodes), 2)

        if S[i] == S[j]:
            J[i, j] = -J[i, j]
            J[j, i] = J[i, j]
        elif S[i] != S[j] and J[i, j] == 1:
            if probability < alpha_value:
                S[i] = S[j] = -1
            J[i, j] = -J[i, j]
            J[j, i] = J[i, j]

        Hs, Hp = calculate_Hs_Hp(i, j, Hs, Hp, S, J, triangles)
        H = calculate_H(nb_nodes, triangles, Hp, Hs)

        if H == -1:
            break

    # Calculer rho_inf - rho0
    rho_inf = calculate_rho_inf(nb_nodes, S)
    rho0 = rho_value
    return (rho_value, alpha_value, rho_inf - rho0)  # Retourne le résultat

def main():
    # Paramètres
    rhos = np.linspace(0.0, 0.5, 10)  # Différents rho0 de 0.0 à 0.5 (réduit)
    alphas = np.linspace(1.0, 0.0, 10)  # Différents alpha de 0.0 à 1.0 (réduit)
    alphas = alphas[::-1]  # Inverser l'ordre de alpha

    # Créer une grille de (rho, alpha)
    params = [(rho_value, alpha_value) for rho_value in rhos for alpha_value in alphas]

    # Utiliser multiprocessing pour paralléliser le calcul
    with Pool(cpu_count()) as pool:  # Utilisation de tous les cœurs disponibles
        results = pool.map(simulate_combination, params)

    # Organiser les résultats dans une structure de données
    rho_diff = np.zeros((len(alphas), len(rhos)))
    for result in results:
        rho_value, alpha_value, diff = result
        i = np.where(rhos == rho_value)[0][0]
        j = np.where(alphas == alpha_value)[0][0]
        rho_diff[j, i] = diff

    # Création de la figure 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Meshgrid pour les valeurs de rho et alpha
    RHO, ALPHA = np.meshgrid(rhos, alphas)

    # Tracer la surface 3D
    ax.plot_surface(RHO, ALPHA, rho_diff, cmap='viridis')

    # Ajouter des étiquettes et un titre
    ax.set_xlabel('Rho')
    ax.set_ylabel('Alpha')
    ax.set_zlabel('Rho_inf - Rho_0')
    ax.set_title('Surface 3D de Rho_inf - Rho_0 en fonction de Rho et Alpha')

    # Afficher le graphique
    plt.show()

# Ajoute cette ligne pour éviter des erreurs sous Windows
if __name__ == '__main__':
    main()
