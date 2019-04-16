from Tabular_values import V_a_e, Na
from math import  pi, log10

from Internal_energy import Energy

m = 401
n = 501

z = 13
Atom_weight = 26.981

TABLE_E = [[0 for i in range(m + 1)] for j in range(n + 1)]
TABLE_E[0][1] = - 2.000
TABLE_E[1][0] = 2.000
for j in range(2, m + 1):
    TABLE_E[0][j] = TABLE_E[0][1] + (j - 1) / m * 4

for i in range(2, n + 1):
    TABLE_E[i][0] = TABLE_E[1][0] + (i - 1) / n * (-10)

for i in range(1, n + 1):
    for j in range(1, m + 1):
        T_now = 10 ** TABLE_E[i][0]
        rho_now = 10 ** TABLE_E[0][j]
        print(T_now, rho_now)
        TABLE_E[i][j] = Energy(T_now, rho_now, Atom_weight, z)
        print("T = ", TABLE_E[i][0], "rh0 = ", TABLE_E[0][j], "E = ", TABLE_E[i][j])







for i in range(len(TABLE_E)):
    for j in range(len(TABLE_E[i])):
        print(TABLE_E[i][j], end=' ')
    print()