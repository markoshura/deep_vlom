from Tabular_values import V_a_e, Na
from math import  pi, log10

from working_progonka import progonka
from Pressure import P
from corrections_P import delta_P

m = 39
n = 48


TABLE_P = [[0 for i in range(m + 1)] for j in range(n + 1)]

TABLE_P[0][1] = - 4.000
TABLE_P[1][0] = 3.000
for j in range(2, m + 1):
    TABLE_P[0][j] = TABLE_P[0][1] + (j - 1) / m * 13



for i in range(2, n + 1):
    TABLE_P[i][0] = TABLE_P[1][0] + (i - 1) / n * (-8)

for i in range (1, n + 1):
    for j in range (1, m + 1):
        T_now = 10 ** TABLE_P[i][0] / 36.7493224786
        rho_now = 1.0 / (8.923608963522 * 0.01 * 10 ** TABLE_P[0][j])
        PHI = progonka(T_now, rho_now, 1, 1)
        TABLE_P[i][j] = log10(abs(delta_P(T_now, rho_now)))
        print ("T = ", TABLE_P[i][0], "rh0 = ", TABLE_P[0][j], "LG P = ", TABLE_P[i][j])






for i in range(len(TABLE_P)):
    for j in range(len(TABLE_P[i])):
        print(TABLE_P[i][j], end=' ')
    print()
