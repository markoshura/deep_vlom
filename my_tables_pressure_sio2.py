from mixture import mixture_calculation
from Pressure import P

from Tabular_values import Na

mixture = []
elem_1_si = {'Atom_weight': , 'Z': , 'mass': , 'amount': 1}
elem_2_o = {'Atom_weight': , 'Z': , 'mass': , 'amount' : 2}

nuclear_amount = 0
for elem in range(len(mixture)):
    nuclear_amount += mixture[elem]["amount"]

mixture.append(elem_1_si)
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




for i in range(1, n + 1):
    for j in range(1, m + 1):
        found_rho = mixture_calculation(10 ** TABLE_E[i][0], 10 ** TABLE_E[0][j], mixture)
        for k in range(len(mixture)):
            T_h = 10 ** TABLE_E[i][0] / mixture[i]["Z"]**(4 / 3)
            rho_h = 10 ** TABLE_E[0][j] * Na * 11.19 * 1.4818 * 10**(-25) / mixture[i]["Atom_weight"] / mixture[i]["Z"]
            TABLE_E[i][j] += P(T_h, rho_h) * #mixture[i]["Z"]**(7 / 3) #* mixture[i]["amount"] / nuclear_amount
        print(i, j)

f = open('text.txt', 'w')
for i in range(len(TABLE_E)):
    for j in range(len(TABLE_E[i])):
        f.write(str(TABLE_E[i][j]))
        f.write(" ")
    f.write("\n")
