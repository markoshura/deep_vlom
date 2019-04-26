from Tabular_values import V_a_e, Na
from math import  pi, log10
from Cell import theta, r_0   , eta_0

from Internal_energy import Energy

m = 401
n = 1

z = 13
Atom_weight = 26

TABLE_E = [[0 for i in range(m + 1)] for j in range(n + 1)]
TABLE_E[0][1] = - 2.000
TABLE_E[1][0] = -9.000
for j in range(2, m + 1):
    TABLE_E[0][j] = TABLE_E[0][1] + (j - 1) * 0.01




for j in range(1, m + 1):
    T_h = 10 ** TABLE_E[1][0] / z**(4 / 3)
    rho_h = 10 ** TABLE_E[0][j] * Na * 11.19 * 1.4818 * 10**(-25) / Atom_weight / z
    TABLE_E[1][j] = Energy(T_h, rho_h, 1, 1) * z**(7 / 3)
    q = 3 * pi / 2 * 1 / 2 / theta(T_h)**(3 / 2) / r_0(rho_h, 1)**3
    #print("q = ", q)
    #print(- (3 / 2 * q) ** (2 / 3))
    print(TABLE_E[1][j])
    #print(eta_0(T_h, rho_h))

#f = open('text4.txt', 'w')
#for i in range(len(TABLE_E)):
#    for j in range(len(TABLE_E[i])):
#        f.write(str(TABLE_E[i][j]))
#        f.write(" ")
#    f.write("\n")
