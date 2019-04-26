from mixture import mixture_calculation
from Pressure import P
from Internal_energy import Energy

from Tabular_values import Na
from working_progonka import progonka

mixture = []
elem_1_al = {'Atom_weight': 26, 'Z': 13, 'mass': 26, 'amount': 2}


elem_2_o = {'Atom_weight': 16, 'Z': 8, 'mass': 16, 'amount': 3}



nuclear_amount = 5

Weight = 100



mixture.append(elem_1_al)

mixture.append(elem_2_o)




m = 401
n = 501



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





for i in range(1, n + 1):
    for j in range(1, m + 1):
        #print(10 ** TABLE_E[i][0], 10 ** TABLE_E[0][j])
        found_rho = mixture_calculation(10 ** TABLE_E[i][0], 10 ** TABLE_E[0][j], mixture)
        for k in range(len(mixture)):
            T_h = 10 ** TABLE_E[i][0] / mixture[k]["Z"]**(4 / 3)
            rho_h = 10 ** TABLE_E[0][j] * Na * 11.19 * 1.4818 * 10**(-25) / mixture[k]["Atom_weight"] / mixture[k]["Z"]
            PHI = progonka(T_h, rho_h, 1, 1)
            TABLE_P[i][j] += P(T_h, rho_h, PHI) * mixture[k]["Z"]**(10 / 3) * mixture[k]["amount"] / nuclear_amount
            TABLE_E[i][j] += Energy(T_h, rho_h, 1, 1, PHI) * mixture[k]["Z"]**(7 / 3) * mixture[k]["amount"] / Weight
        print(i, j)

f = open('energy_al2o3.txt', 'w')
for i in range(len(TABLE_E)):
    for j in range(len(TABLE_E[i])):
        f.write(str(TABLE_E[i][j]))
        f.write(" ")
    f.write("\n")

f = open('pressure_al2o3.txt', 'w')
for i in range(len(TABLE_P)):
    for j in range(len(TABLE_P[i])):
        f.write(str(TABLE_P[i][j]))
        f.write(" ")
    f.write("\n")

