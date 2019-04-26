from Tabular_values import V_a_e, Na
from Pressure import P

m = 401
n = 501

z = 26
Atom_weight = 55.84

TABLE_E = [[0 for i in range(m + 1)] for j in range(n + 1)]
TABLE_E[0][1] = - 2.000
TABLE_E[1][0] = 2.000
for j in range(2, m + 1):
    TABLE_E[0][j] = TABLE_E[0][1] + (j - 1) * 0.01

for i in range(2, n + 1):
    TABLE_E[i][0] = TABLE_E[1][0] + (i - 1) * (-0.02)

for i in range(1, n + 1):
    for j in range(1, m + 1):
        T_h = 10 ** TABLE_E[i][0] / z**(4 / 3)
        rho_h = 10 ** TABLE_E[0][j] * Na * 11.19 * 1.4818 * 10**(-25) / Atom_weight / z
        TABLE_E[i][j] = P(T_h, rho_h) * z**(10 / 3)
        print(i, j)

f = open('text1.txt', 'w')
for i in range(len(TABLE_E)):
    for j in range(len(TABLE_E[i])):
        f.write(str(TABLE_E[i][j]))
        f.write(" ")
    f.write("\n")
