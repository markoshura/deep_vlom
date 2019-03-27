from Tabular_values import V_a_e, Na
from math import  pi, log10

from Internal_energy import Energy

m = 50
n = 100


TABLE_E = [[0 for i in range(m + 1)] for j in range(n + 1)]
#f = open('pressure_output.txt', 'w')
#LG_T = [3.000, 2.833, 2.667, 2.500, 2.333, 2.167, 2.000, 1.833, 1.667, 1.500, 1.333, 1.167, 1.000, 0.833, 0.667, 0.500, 0.333, 0.167, 0.000]
#LG_V = [-4.000]
#TA
TABLE_E[0][1] = 10**(-2) #плотность по 100 точек на порядок
TABLE_E[1][0] = 10**(2) #температура кэВ по 50 точек на порядок
i = 0
while i < 4:

    for j in range(2, (m + 1) // 4 * (i + 1)):
        TABLE_E[0][j] = TABLE_E[0][1] + (j - 1) * 0.001

    i += 1

j = 0
while j < 4:

    for i in range(2, (n + 1) //4 * (j + 1)):
        TABLE_E[i][0] = TABLE_E[1][0] - (i - 1) / n * 0.02
    j += 1
#
#for i in range(1, n + 1):
#    for j in range(1, m + 1):
#        T_now = 10 ** TABLE_E[i][0] / 36.7493224786
#        rho_now = 1.0 / (8.923608963522 * 0.01 * 10 ** TABLE_E[0][j])
#        # print(T_now, rho_now)
#        # print(P(T_now, rho_now))                       к
#        TABLE_E[i][j] = log10(Energy(T_now, rho_now))
#        print("T = ", TABLE_E[i][0], "rh0 = ", TABLE_E[0][j], "LG E = ", TABLE_E[i][j])
#T_now = 10 ** (-0.167) /36.749
#rho_now = 1.0 / (8.923608963522 * 0.01 * 10 ** 4)
#print(P(T_now, rho_now))
#print(log10(P(T_now,rho_now)))






for i in range(len(TABLE_E)):
    for j in range(len(TABLE_E[i])):
        print(TABLE_E[i][j], end=' ')
    print()