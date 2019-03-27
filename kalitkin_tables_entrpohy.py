from Tabular_values import V_a_e, Na
from math import  pi, log10

from Entrophy import S

m = 39
n = 48


TABLE_S = [[0 for i in range(m + 1)] for j in range(n + 1)]
#f = open('pressure_output.txt', 'w')
#LG_T = [3.000, 2.833, 2.667, 2.500, 2.333, 2.167, 2.000, 1.833, 1.667, 1.500, 1.333, 1.167, 1.000, 0.833, 0.667, 0.500, 0.333, 0.167, 0.000]
#LG_V = [-4.000]
#TA
TABLE_S[0][1] = - 4.000
TABLE_S[1][0] = 3.000
for j in range(2, m + 1):
    TABLE_S[0][j] = TABLE_S[0][1] + (j - 1) / m * 13

for i in range(2, n + 1):
    TABLE_S[i][0] = TABLE_S[1][0] + (i - 1) / n * (-8)

for i in range(1, n + 1):
    for j in range(1, m + 1):
        T_now = 10 ** TABLE_S[i][0] / 36.7493224786
        rho_now = 1.0 / (8.923608963522 * 0.01 * 10 ** TABLE_S[0][j])
        # print(T_now, rho_now)
        # print(P(T_now, rho_now))
        TABLE_S[i][j] = log10(S(T_now, rho_now))
        print("T = ", TABLE_S[i][0], "rh0 = ", TABLE_S[0][j], "LG S = ", TABLE_S[i][j])
#T_now = 10 ** (-0.167) /36.749
#rho_now = 1.0 / (8.923608963522 * 0.01 * 10 ** 4)
#print(P(T_now, rho_now))
#print(log10(P(T_now,rho_now)))






for i in range(len(TABLE_S)):
    for j in range(len(TABLE_S[i])):
        print(TABLE_S[i][j], end=' ')
    print()
#
#
#        T_now = 10 ** LG_T[j] / 36.7493224786
#        rho_now = 1.0 / (8.923608963522 * 0.01 * 10 ** LG_V[i])
#        f.write(str(log10(P(T_now, rho_now))) + '\n')
#        #print('lg v = ', LG_V[i], 'lg T = ', LG_T[j], 'lg P = ', log10(P(T_now, rho_now)))
#
#T_now = 10 ** (3)/ 36.7493224786
#rho_now = 1.0/(8.923608963522 * 0.01 * 10 ** 0.333)
#
#print("P_e = ", log10(P_e(T_now, rho_now)))
#
#
#print('theta = ', theta(T_now))
#print('eta = ', eta(T_now, rho_now))
#print('Atom_weight = ', Atom_weight)
#print('z = ', z)
#print(log10(10))
#

