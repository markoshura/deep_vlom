from Tabular_values import  Na


from Internal_energy import Energy
from Pressure import P
from working_progonka import progonka
m = 401
n = 501

z = 26
Atom_weight = 55.845

TABLE_E = [[0 for i in range(m + 1)] for j in range(n + 1)]
TABLE_E[0][1] = - 2.000
TABLE_E[1][0] = 2.000
for j in range(2, m + 1):
    TABLE_E[0][j] = TABLE_E[0][1] + (j - 1) * 0.01

for i in range(2, n + 1):
    TABLE_E[i][0] = TABLE_E[1][0] + (i - 1) * (-0.02)


TABLE_P = [[0 for i in range(m + 1)] for j in range(n + 1)]
TABLE_P[0][1] = - 2.000
TABLE_P[1][0] = 2.000
for j in range(2, m + 1):
    TABLE_P[0][j] = TABLE_P[0][1] + (j - 1) * 0.01

for i in range(2, n + 1):
    TABLE_P[i][0] = TABLE_P[1][0] + (i - 1) * (-0.02)

for i in range(167, n + 1):
    for j in range(1, m + 1):
        T_h = 10 ** TABLE_E[i][0] / z**(4 / 3)
        rho_h = 10 ** TABLE_E[0][j] * Na * 11.19 * 1.4818 * 10**(-25) / Atom_weight / z
        PHI = progonka(T_h, rho_h, 1, 1)
        TABLE_E[i][j] = Energy(T_h, rho_h, 1, 1, PHI) * z ** (7 / 3)
        TABLE_P[i][j] = P(T_h, rho_h, PHI) * z ** (10 / 3)
        print(i, j)

f = open('energy_ferrum.txt', 'w')
for i in range(len(TABLE_E)):
    for j in range(len(TABLE_E[i])):
        f.write(str(TABLE_E[i][j]))
        f.write(" ")
    f.write("\n")

f = open('pressure_ferrum.txt', 'w')
for i in range(len(TABLE_P)):
    for j in range(len(TABLE_P[i])):
        f.write(str(TABLE_P[i][j]))
        f.write(" ")
    f.write("\n")





