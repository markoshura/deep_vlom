from math import pi
from Tabular_values import Na, a_0
from Changing_parameters import Temperature_system, rho_system
from working_progonka import progonka


# Число итераций
p = 5

mixture = []
elem_1 = {'Atom_weight': 1, 'Z': 1, 'mass': 1}
elem_2 = {'Atom_weight': 2, 'Z': 2, 'mass': 2}
mixture.append(elem_1)
mixture.append(elem_2)
amount_of_elements = len(mixture)

small_const_b = (3 / 4 / pi / Na)**(1 / 3) / a_0

for i in range(len(mixture)):
    big_const_b = small_const_b**3 * (rho_system**(- 1) * mixture[i]["mass"])





def mixture_calculation(Temperature_system, rho_system, mixture):

# 0-e приближение r_0**3
    def x_0(Atom_weight):
        return Atom_weight * rho_system**(- 1)* small_const_b**3
    def mu_0(Atom_weight, z):
        PHI = progonka(Temperature_system, rho_system, Atom_weight, z)
        return PHI[len(PHI) - 1]

    def mu_shtrih_0(Atom_weight, z):
        x0 = x_0(Atom_weight)
        x1 = x_0(Atom_weight + Atom_weight * 0.01)
        PHI_0 = progonka(Temperature_system, z / 4 / 3 / pi / x0, Atom_weight, z)
        mu0 = PHI_0[len(PHI_0) - 1]
        PHI_1 = progonka(Temperature_system, z / 4 / 3 / pi / x1, Atom_weight + Atom_weight * 0.01, z)
        mu1 = PHI_1[len(PHI_1) - 1]
        return (mu1 - mu0) / (x1 - x0)

    def mu_shtrih(i, p):
        return (mu_elems_iterations[i][p] - mu_elems_iterations[i][p - 1]) / (
                    x_elems_iterations[i][p] - x_elems_iterations[i][p - 1])

    def mu(x, Atom_weight, z):
        PHI = progonka(Temperature_system, 1, Atom_weight, z)
        return PHI[len(PHI) - 1]

    def Y(p):
        promezh_summ_1 = 0
        promezh_summ_2 = 0
        for i in range(len(mixture)):
            promezh_summ_1 += mixture[i]["mass"] / mixture[i]["Atom_weight"] / mu_shtrih_elems_iterations[i][p]
            print('i = ', i, 'promezh_summ_1 = ', promezh_summ_1, 'p = ', p)
        for i in range(len(mixture)):
            promezh_summ_2 += mixture[i]["mass"] / mixture[i]["Atom_weight"] * (mu_elems_iterations[i][p] / mu_shtrih_elems_iterations[i][p] - x_elems_iterations[i][p]) * (promezh_summ_1)**(- 1)

        return big_const_b + promezh_summ_2

    # Хим потенциалы элемента i на итерации p


    mu_elems_iterations = [[0 for i in range(p + 1)] for j in range(len(mixture))]

    # Искомый куб r_0 элемента i на итерации p
    x_elems_iterations = [[0 for i in range(p + 1)] for j in range(len(mixture))]

    #Производная хим потенциала
    mu_shtrih_elems_iterations = [[0 for i in range(p + 1)] for j in range(len(mixture))]

    for i in range(len(mixture)):
        mu_elems_iterations[i][0] = mu_0(mixture[i]["Atom_weight"], mixture[i]["Z"])
        x_elems_iterations[i][0] = x_0(mixture[i]['Atom_weight'])
        mu_shtrih_elems_iterations[i][0] = mu_shtrih_0(mixture[i]['Atom_weight'], mixture[i]['Z'])

    p_current = 1

    while p_current < p:

        for i in range(len(mixture)):
            x_elems_iterations[i][p_current] = (Y(p_current - 1) - mu_elems_iterations[i][p_current - 1]) / mu_shtrih_elems_iterations[i][p_current - 1] + x_elems_iterations[i][p_current - 1]
            mu_elems_iterations[i][p_current] = mu(x_elems_iterations[i][p_current], mixture[i]["Atom_weight"], mixture[i]["Z"])
            mu_shtrih_elems_iterations[i][p_current] = mu_shtrih(i, p_current)

        p_current += 1

    return x_elems_iterations[p]

print(mixture_calculation(Temperature_system, rho_system, mixture))
#print(len(mixture))



