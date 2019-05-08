from math import pi
from Tabular_values import Na, a_0
#from Changing_parameters import T, rho
from working_progonka import progonka


mixture = []
elem_1_si = {'Atom_weight': 28, 'Z': 14, 'mass': 28, 'amount': 1}
elem_2_o = {'Atom_weight': 16, 'Z': 8, 'mass': 16, 'amount': 2}
elem_3_o = {'Atom_weight': 16, 'Z': 8, 'mass': 16, 'amount': 2}



mixture.append(elem_1_si)
mixture.append(elem_2_o)
mixture.append(elem_3_o)
def mixture_calculation(T, rho, mixture):
    # Число итераций
    p = 5
    amount_of_elements = len(mixture)

    small_const_b = (3 / 4 / pi / Na)**(1 / 3) / a_0

    big_const_b = 0
    for i in range(len(mixture)):
        big_const_b += small_const_b**3 * (rho**(- 1) * mixture[i]["mass"])

# 0-e приближение r_0**3
    def x_0(Atom_weight):
        return Atom_weight * rho**(- 1)* small_const_b**3
    def mu_0(Atom_weight, z):
        PHI = progonka(T, rho, Atom_weight, z)
        return PHI[len(PHI) - 1]

    def mu_shtrih_0(Atom_weight, z):
        x0 = x_0(Atom_weight)
        x1 = x_0(Atom_weight + Atom_weight * 0.1)
        PHI_0 = progonka(T, 3 * Atom_weight / x0 / 4 / pi / Na / a_0**3, Atom_weight, z)
        mu0 = PHI_0[len(PHI_0) - 1]
        PHI_1 = progonka(T, 3 * Atom_weight / x1 / 4 / pi / Na / a_0**3, Atom_weight,  z)
        mu1 = PHI_1[len(PHI_1) - 1]
        return (mu1 - mu0) / (x1 - x0)

    def mu_shtrih(i, p):
        return (mu_elems_iterations[i][p] - mu_elems_iterations[i][p - 1]) / (
                    x_elems_iterations[i][p] - x_elems_iterations[i][p - 1])

    def mu(x, Atom_weight, z):
        PHI = progonka(T, 3 * Atom_weight / x / 4 / pi / Na / a_0**3, Atom_weight, z)
        return PHI[len(PHI) - 1]
    def Y(p):
        promezh_summ_1 = 0
        promezh_summ_2 = 0
        for i in range(len(mixture)):
            promezh_summ_1 += mixture[i]["mass"] / mixture[i]["Atom_weight"] / mu_shtrih_elems_iterations[i][p]
            #print('i = ', i, 'promezh_summ_1 = ', promezh_summ_1, 'p = ', p)
        for i in range(len(mixture)):
            promezh_summ_2 += mixture[i]["mass"] / mixture[i]["Atom_weight"] * (mu_elems_iterations[i][p] / mu_shtrih_elems_iterations[i][p] - x_elems_iterations[i][p])

        return (big_const_b + promezh_summ_2) * (promezh_summ_1)**(- 1)

    # Хим потенциалы элемента i на итерации p


    mu_elems_iterations = [[0 for i in range(p + 1)] for j in range(len(mixture))]

    # Искомый куб r_0 элемента i на итерации p
    x_elems_iterations = [[0 for i in range(p + 1)] for j in range(len(mixture))]

    #Производная хим потенциала
    mu_shtrih_elems_iterations = [[0 for i in range(p + 1)] for j in range(len(mixture))]

    for i in range(len(mixture)):
        mu_elems_iterations[i][0] = mu_0(mixture[i]["Atom_weight"], mixture[i]["Z"])
        #print("mu_elems_iterations[i][0] = ", mu_elems_iterations[i][0])
        x_elems_iterations[i][0] = x_0(mixture[i]['Atom_weight'])
        #print("x_elems_iterations[i][0] = ", x_elems_iterations[i][0])
        mu_shtrih_elems_iterations[i][0] = mu_shtrih_0(mixture[i]['Atom_weight'], mixture[i]['Z'])
        #print("mu_shtrih_elems_iterations[i][0] = ", mu_shtrih_elems_iterations[i][0])

    p_current = 1

    while p_current <= p:

        for i in range(len(mixture)):
            x_elems_iterations[i][p_current] = (Y(p_current - 1) - mu_elems_iterations[i][p_current - 1]) / mu_shtrih_elems_iterations[i][p_current - 1] + x_elems_iterations[i][p_current - 1]
            #print("x_elems_iterations[i][p_current] = ", x_elems_iterations[i][p_current])
            mu_elems_iterations[i][p_current] = mu(x_elems_iterations[i][p_current], mixture[i]["Atom_weight"], mixture[i]["Z"])
            #print("mu_elems_iterations[i][p_current] = ", mu_elems_iterations[i][p_current])
            mu_shtrih_elems_iterations[i][p_current] = mu_shtrih(i, p_current)


        p_current += 1

    rho_elems_iterations = x_elems_iterations
    for i in range(amount_of_elements):
        for j in range(p + 1):
            rho_elems_iterations[i][j] = rho_elems_iterations[i][j]**(1/3)

    found_rho = [0 for i in range(amount_of_elements)]
    for i in range(amount_of_elements):
        found_rho[i] = rho_elems_iterations[i][p]
    return found_rho
    #return rho_elems_iterations

#print(mixture_calculation(10, 1, mixture))
#print(len(mixture))



