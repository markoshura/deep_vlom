from math import pi, log10
from working_progonka import progonka, X
from Dirak_functions import integral_3_2
from Changing_parameters import N
from Cell import volume, theta, r_0
from Atom_parameters import z
from Entrophy import S
from Atom_parameters import E_0
from scipy import integrate
from Tabular_values import E_h
Z = [(i / N) for i in range(N + 1)]



def E_sub_int(T, rho):
    max_i = 2
    while PHI[max_i] / X[max_i] >= 10 ** 6:
        max_i += 1

    if max_i % 2 != 0:
        max_i -= 1

    integr_int_1 = [Z[i] for i in range(max_i + 1)]
    integr_int_2 = [Z[i] for i in range(max_i, N + 1)]
    # integr_int_2 = [X[i] for i in range(max_i, N + 1)]

    integr_func_1 = [4 / 5 * PHI[i] ** 2.5 for i in range(max_i + 1)]
    integr_func_2 = [2 * Z[i] ** 5 * integral_3_2(PHI[i] / Z[i] ** 2) for i in range(max_i, N + 1)]
    # integr_func_2 = [X[i] ** 2 * integral_3_2(PHI[i] / X[i]) for i in range(max_i, N + 1)]

    # i = 0
    # integr_res_1 = 0
    # while i <= (max_i - 2):
    #
    #    integr_res_1 += 1/6 * (integr_func_1[i] + 4 * integr_func_1[i+1] + integr_func_1[i+2]) * (integr_int_1[i+2] - integr_int_1[i])
    #    i += 2
    #
    # integr_res_2 = 0
    # i = max_i
    # while i <= N - 2:
    #
    #    integr_res_2 += 1/6 * (integr_func_2[i - max_i] + 4 * integr_func_2[i+1 - max_i] + integr_func_2[i+2-max_i]) * (integr_int_2[i+2-max_i] - integr_int_2[i-max_i])
    #    i += 2

    integr_res_1 = integrate.simps(integr_func_1, integr_int_1)

    integr_res_2 = integrate.simps(integr_func_2, integr_int_2)

    return integr_res_1 + integr_res_2

T_now = 10 ** (-2) / 36.7493224786
rho_now = 1.0 / (8.923608963522 * 0.01 * 10 ** (8))

a = 2**(7 / 6) * 3** (2 / 3) * pi** (-5 / 3) * theta(T_now)**(1 / 2) * volume(rho_now)** (2 / 3)

PHI = progonka(T_now, rho_now)
S1 = S(T_now, rho_now)
#E2 = 2 * (2 / pi)**(2/3) * a**(1 / 3) * PHI[0]**(-7/3) * (2 / 3 * integral_3_2(PHI[N]) - E_sub_int(T_now, rho_now)) - E_0

print(log10(S1))

