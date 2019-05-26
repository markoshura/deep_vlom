#from bilinear_interpolation_energy import bilinear_interpolation_energy
#from ionic_contribution import ionic_contribution_energy
#from cold_curve import cold_curve_energy
#from Cell import volume
from Internal_energy import Energy
from math import log10



#n = 100
#cold_curve = [[0 for value in range(n + 1)] for j in range(2)]
#Интерполируем тепловой вклад

#bilinear_interpolation(T, rho)

#Из функции делаем ключ-значение

#def func_cold_curve():
#    cold_curve = [[0 for value in range(n + 1)] for j in range(2)]
#
#
#
#
#    return cold_curve
#
#
##Интерполируем холодную кривую
#
#def interp_cold_curve(rho):
#    cold_curve = func_cold_curve()
#    result = 0
#    for i in range(n + 1):
#        for j in range(1, n + 1):
#            if i == j:
#                result += cold_curve[1][i]
#            else:
#                result += cold_curve[1][i] * (rho - cold_curve[0][j]) / (cold_curve[0][i] - cold_curve[0][j])
#


#Прибавляем ионный вклад

def state_function_energy(T, rho, Atom_weight, z):

    #return 2.626 * 10**3 / Atom_weight * bilinear_interpolation_energy(T, rho) + ionic_contribution_energy(T, volume(rho, Atom_weight)) + cold_curve_energy(volume(rho, Atom_weight))
    return 2.626 * 10**3 / Atom_weight * Energy(T, rho, z, Atom_weight)


print("S_f_e = ", log10(state_function_energy(1, 1, 1, 1)))