from math import pi
from scipy import integrate
from Atom_parameters import E_0, Atom_weight, z
from Changing_parameters import N , Temperature_system, rho_system
from Dirak_functions import integral_3_2
from Cell import volume, theta, r_0
from working_progonka import X, progonka, eta
#from sveryaem_phi import excel_phi as PHI, excel_x as X
#from Pressure import P_e
#from Changing_parameters import Temperature_system, rho_system


Z = [(i / N) for i in range(N + 1)]
# ВСПОМ. ИНТЕГРАЛ
def E_sub_int(T, rho):

    PHI = progonka(T, rho)

    max_i = 1
    while PHI[max_i] / X[max_i] >= 10 ** 6:
        max_i += 1

    if max_i % 2 != 0:
        max_i -= 1

    integr_int_1 = [Z[i] for i in range(max_i + 1)]
    integr_int_2 = [Z[i] for i in range(max_i, N + 1)]
    # integr_int_2 = [X[i] for i in range(max_i, N + 1)]

    integr_func_1 = [4 / 5 * PHI[i] ** 2.5 for i in range(max_i + 1)]
    integr_func_2 = [2 * Z[i] ** 5 * integral_3_2(PHI[i] / Z[i] ** 2) for i in range(max_i, N + 1)]
    print('f_1 = ', integr_func_1)
    print("f_2 = ", integr_func_2)

    print(integr_func_1)


    integr_res_1 = integrate.simps(integr_func_1, integr_int_1)
    print('res_1 = ', integr_res_1)

    integr_res_2 = integrate.simps(integr_func_2, integr_int_2)
    print('res_2 = ', integr_res_2)


PHI = progonka(0.2, 0.001)
max_i = 1
while PHI[max_i] / X[max_i] >= 10 ** 6:
    max_i += 1

if max_i % 2 != 0:
    max_i -= 1

#integr_func_2 = [2 * Z[i] ** 5 * integral_3_2(PHI[i] / Z[i] ** 2) for i in range(max_i, N + 1)]

print(max_i)