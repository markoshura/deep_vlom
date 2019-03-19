from math import pi
from Tabular_values import Na, a_0
from Changing_parameters import Temperature_system, rho_system

# Элементов в смеси
amount_of_elements = 5
# Число итераций
p = 5

mixture = []
elem_1 = {'Atom_weight': 1, 'Z': 1, 'mass': 1}
elem_2 = {'Atom_weight': 2, 'Z': 3, 'mass': 4}
mixture.append(elem_1)
mixture.append(elem_2)

small_const_b = (3 / 4 / pi / Na)**(1 / 3) / a_0

whole_mass = 0
for i in range(len(mixture)):
    big_const_b = small_const_b**3 * (rho_system**(- 1) * mixture[i]["mass"])

#0-e приближение r_0**3



def mixture_calculation(T, rho, mixture):

    def x_0(Atom_weight, rho_system):
        return (1.388) ** 3 * (Atom_weight / rho_system)

    def mu_shtrhih(i, p):
        return (mu_elems_iterations[p][i] - mu_elems_iterations[p - 1][i]) / (
                    x_elems_iterstions[p][i] - x_elems_iterstions[p - 1][i])

    def Y(p):
        promezh_summ_1 = 0
        promezh_summ_2 = 0
        for i in range(len(mixture)):
            promezh_summ_1 += mixture[i]["mass"] / mixture[i]["Atom_weight"] / mu_shtrih_elems_iterations[p][i]
            promezh_summ_2 += mixture[i]["mass"] / mixture[i]["Atom_weight"] * (mu_elems_iterations[p][i] / mu_shtrih_elems_iterations[p][i] - x_elems_iterstions[p][i]) * (promezh_summ_1)**(- 1)

        return big_const_b + promezh_summ_2

    # Хим потенциалы элемента i на итерации p
    mu_elems_iterations = [[0 for i in range(amount_of_elements + 1)] for j in range(p + 1)]

    # Искомый куб r_0 элемента i на итерации p
    x_elems_iterstions = [[0 for i in range(amount_of_elements + 1)] for j in range(p + 1)]

    #Производная хим потенциала
    mu_shtrih_elems_iterations = [[0 for i in range(amount_of_elements + 1)] for j in range(p + 1)]

    for i in range(len(mixture)):
        mu_elems_iterations[0][i] =
        x_elems_iterstions[0][i] = x_0(mixture[i]['Atom_weight'])
        mu_shtrih_elems_iterations[0][i] =

    p_current = 0

    while p_current < p:

        for i in range(len(mixture)):
            x_elems_iterstions[p_current+1][i] = (Y(p_current) - mu_elems_iterations[p_current][i]) / mu_shtrih_elems_iterations[p_current][i] + x_elems_iterstions[p_current][i]


