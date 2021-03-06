# ВСЁ ДЛЯ ЛИНЕАРИЗАЦИИ И ПРОГОНКИ С ИТЕРАЦИЯМИ

import matplotlib.pyplot as plt

from math import pi
from Changing_parameters import N, Temperature_system, rho_system
from Dirak_functions import integral_1_2, integral_minus_1_2, igrek, integral_3_2_1
#from Atom_parameters import z
from Cell import r_0, theta, eta_0

# Число итераций
s = 9

# СЕТКА
X = [(i / N) ** 2 for i in range(N + 1)]



# ШАГ СЕТКИ
h_N = 1 - X[N - 1]

# КОЭФФИЦИЕНТЫ
A = [0] + [1 / (X[i] - X[i - 1]) for i in range(1, N)]
C = [0] + [1 / (X[i + 1] - X[i]) for i in range(1, N)]


def progonka(T, rho, Atom_weight, z):
    # константа а в уравнении
    const_a = 4 * (2 * theta(T)) ** 0.5 / pi * (r_0(rho, Atom_weight)) ** 2

    PHI_S = [[0 for i in range(N + 1)] for j in range(s + 1)]
    RESULT = [[0 for i in range(N + 1)] for j in range(s)]


    for i in range(N + 1):
        if i == 0:
            PHI_S[0][0] = z / (theta(T) * r_0(rho, Atom_weight))
        else:
            PHI_S[0][i] = z / (theta(T) * r_0(rho, Atom_weight)) * (1 - 3 / 2 * X[i] + 1 / 2 * (X[i]) ** 3) - eta_0(T, rho) * X[i]

    s_current = 0

    while s_current < s:
        B = [0] + [-A[i] - C[i] - const_a / 4 * (X[i + 1] - X[i - 1]) * integral_minus_1_2(PHI_S[s_current][i] / X[i])
                   for i in range(1, N)]
        D = [0] + [const_a * (X[i + 1] - X[i - 1]) * ((X[i] / 2 * integral_1_2(PHI_S[s_current][i]/X[i])) -
                                                      PHI_S[s_current][i] / 4 * integral_minus_1_2(PHI_S[s_current][i] /
                    X[i])) for i in range(1, N)]

        ALPHA = [0] * (N + 1)
        BETA = [0] * (N + 1)
        for i in range(N - 1, -1, -1):
            if i == N - 1:
                ALPHA[i] = 1 / (1 - h_N + h_N**2 / 4 * const_a * integral_minus_1_2(PHI_S[s_current][N]))
                BETA[i] = -h_N**2 / 2 * const_a * (
                        integral_1_2(PHI_S[s_current][N]) - PHI_S[s_current][i] / 2 * integral_minus_1_2(
                    PHI_S[s_current][N])) * ALPHA[N - 1]
            else:
                ALPHA[i] = -A[i + 1] / (B[i + 1] + C[i + 1] * ALPHA[i + 1])
                BETA[i] = (D[i + 1] - C[i + 1] * BETA[i + 1]) / (B[i + 1] + C[i + 1] * ALPHA[i + 1])

        # ИСКОМОЕ УРАВНЕНИЕ
        Y = [0] * (N + 1)
        for i in range(N + 1):
            if i == 0:
                Y[i] = z / (theta(T) * r_0(rho, Atom_weight))
            else:
                Y[i] = ALPHA[i - 1] * Y[i - 1] + BETA[i - 1]

        for i in range(N + 1):
            RESULT[s_current][i] = PHI_S[s_current][i] * theta(T)

        for i in range(N + 1):
            PHI_S[s_current+1][i] = Y[i]

        s_current += 1
    PHI = [RESULT[8][i] / theta(T) for i in range(N + 1)]

    ##СТРОИМ ГРАФИК
    #fig = plt.figure()
    #for i in range(s):
    #    graph1 = plt.plot(X, PHI)
    #plt.title("T = 0 keV")
    #plt.grid(True)
    #
    #plt.xlabel('x')
    #plt.ylabel('y(s)')
    #plt.savefig('progonka1')
    #plt.show()

    return PHI


#progonka(Temperature_system, rho_system)

def eta(T, rho, Atom_weight, z):
    PHI = progonka(T, rho, Atom_weight, z)
    return - PHI[N]









