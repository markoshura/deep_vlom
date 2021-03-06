from bilinear_interpolation_pressure import bilinear_interpolation_pressure
from ionic_contribution import  ionic_contribution_pressure
from cold_curve import cold_curve_pressure
from Cell import volume

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


def state_function_pressure(T, rho, Atom_weight):
    #bilinear_interpolation_pressure(T, rho)

    return 2.942 * 10**4 *bilinear_interpolation_pressure(T, rho) + ionic_contribution_pressure(T, volume(rho, Atom_weight)) + cold_curve_pressure(volume(rho, Atom_weight))


for i in range(2, 502):
    T_now = 2 + (i - 1) * (-0.02)
    print(state_function_pressure(10**T_now, 10, 55.84))
    #print(T_now)
