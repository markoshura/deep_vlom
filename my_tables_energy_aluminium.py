from math import  pi, log10

from Internal_energy import Energy

m1 = -2
m2 = 2
n1 = -8
n2 = -2


TABLE_E = [[0 for i in range((m2 - m1 + 1) * 90)] for j in range((n2 - n1 + 1) * 90)]
#f = open('pressure_output.txt', 'w')
#LG_T = [3.000, 2.833, 2.667, 2.500, 2.333, 2.167, 2.000, 1.833, 1.667, 1.500, 1.333, 1.167, 1.000, 0.833, 0.667, 0.500, 0.333, 0.167, 0.000]
#LG_V = [-4.000]
#TA
TABLE_E[0][1] = 10 ** m1 #плотность
TABLE_E[0][(m2 - m1 + 1) * 91 - 1] = 10 ** m2
TABLE_E[1][0] = 10 ** n1 #температура
TABLE_E[(n2 - n1 + 1) * 91 - 1][0] = 10 ** n2


#k = TABLE_E[0][1]
#while (k < TABLE_E[0][(m2 - m1 + 1) * 100 - 1]):
#    iter = k
#    zapomni = k
#    for i in range(100):
#        if iter < TABLE_E[0][(m2 - m1 + 1) * 100 - 1]:
#            iter = zapomni + zapomni * i
#            print(iter)
#            j += 1
#    k = k * 10

#number = 0
#while number != m2 - m1 + 1:
#    k = TABLE_E[0][1]
#    iter = k
#    zapomni = k
#    for i in range(99):
#        if iter < k * 10:
#            iter = k + k * i / 10
#            print(iter)
#    number += 1

d = 0
j = 2
while d < m2 - m1:
    i = 10** (m1 + d)
    z = i
    for k in range(1, 91):
        if i < 10**(m1 + d + 1):
            i = 10 ** (m1 + d) + 10**(m1 + d - 1) * k
            TABLE_E[0][j] = i
            j += 1

    d += 1

for i in range(len(TABLE_E)):
    for j in range(len(TABLE_E[i])):
        print(TABLE_E[i][j], end=' ')
    print()












