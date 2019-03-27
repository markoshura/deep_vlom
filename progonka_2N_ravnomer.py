# ВСЁ ДЛЯ ЛИНЕАРИЗАЦИИ И ПРОГОНКИ С ИТЕРАЦИЯМИ

import matplotlib.pyplot as plt

from math import pi
from Changing_parameters import N, Temperature_system, rho_system
from Dirak_functions import integral_1_2, integral_minus_1_2
from Atom_parameters import z
from Cell import r_0, theta, eta_0

N = 2 * N

# Число итераций
s = 9

# СЕТКА
U = [(i / N) for i in range(N + 1)]

X = [U[i]**2 for i in range(N + 1)]

# ШАГ СЕТКИ
h_N = 1 / N

# КОЭФФИЦИЕНТЫ
A = [0] + [2 * U[i] + h_N for i in range(1, N)]
C = [0] + [2 * U[i] - h_N for i in range(1, N)]


def progonka_2N(T, rho):
    # константа а в уравнении
    const_a = 4 * (2 * theta(T)) ** 0.5 / pi * (r_0(rho)) ** 2

    PHI_S = [[0 for i in range(N + 1)] for j in range(s + 1)]
    RESULT = [[0 for i in range(N + 1)] for j in range(s)]


    for i in range(N + 1):
        if i == 0:
            PHI_S[0][0] = z / (theta(T) * r_0(rho))
        else:
            PHI_S[0][i] = z / (theta(T) * r_0(rho)) * (1 - 3 / 2 * U[i] + 1 / 2 * (U[i]) ** 3) - eta_0(T, rho) * U[i]

    s_current = 0

    while s_current < s:
        B = [0] + [-4 * U[i]*(1 + const_a * h_N**2 * U[i]**2 * integral_minus_1_2(PHI_S[s_current][i] / U[i]**2))
                   for i in range(1, N)]
        D = [0] + [4 * const_a * h_N**2 * U[i]**3 * (2 * U[i]**2 * integral_1_2(PHI_S[s_current][i] /
                    U[i]**2) - PHI_S[s_current][i] * integral_minus_1_2(PHI_S[s_current][i] /
                    U[i]**2)) for i in range(1, N)]

        ALPHA = [0] * (N + 1)
        BETA = [0] * (N + 1)
        for i in range(N - 1, -1, -1):
            if i == N - 1:
                ALPHA[i] = 1 / (1 - 2 * h_N + h_N**2 * (1 +  const_a * integral_minus_1_2(PHI_S[s_current][N])))
                BETA[i] = -h_N**2 * const_a * (2 * integral_1_2(PHI_S[s_current][N]) - PHI_S[s_current][i] * integral_minus_1_2(
                    PHI_S[s_current][N])) * ALPHA[N - 1]
            else:
                ALPHA[i] = -A[i + 1] / (B[i + 1] + C[i + 1] * ALPHA[i + 1])
                BETA[i] = (D[i + 1] - C[i + 1] * BETA[i + 1]) / (B[i + 1] + C[i + 1] * ALPHA[i + 1])

        # ИСКОМОЕ УРАВНЕНИЕ
        Y = [0] * (N + 1)
        for i in range(N + 1):
            if i == 0:
                Y[i] = z / (theta(T) * r_0(rho))
            else:
                Y[i] = ALPHA[i - 1] * Y[i - 1] + BETA[i - 1]

        for i in range(N + 1):
            RESULT[s_current][i] = PHI_S[s_current][i] * theta(T)

        for i in range(N + 1):
            PHI_S[s_current+1][i] = Y[i]

        s_current += 1
    PHI = [RESULT[8][i] / theta(T) for i in range(N + 1)]

    return PHI

    ##СТРОИМ ГРАФИК
    #fig = plt.figure()
    #for i in range(s):
    #    graph1 = plt.plot(X, PHI)
    #plt.title("T = 0 keV")
    #plt.grid(True)
    #
    #plt.xlabel('U')
    #plt.ylabel('y(s)')
    #plt.show()
















