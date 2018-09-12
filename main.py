# deep_vlom


import math
import numpy as np
import matplotlib.pyplot as plt

r_0 = 8.08
eta = [2.78, 4.75, 6.99, 10.1]
T = [0.01, 0.1, 1, 10]
# Z_0 = [3.80, 17.1, 57.8, 77.9]
# r_* = [2.63, 2.58, 0.83, 0.11]
# theta =

# I_1_2 = 1/2*pow(math.pi, 1/2)*pow(math.e, x)
# I_3_2 = 3/4*pow(math.pi, 1/2)*pow(math.e, x)
# I_k = pow(x, k+1)/(k+1)
A = 196.96657
rho = 1

# rho_e = pow((2*theta),(3/2))*(I(-eta)_(1/2))/(2*pow(math.pi,2))
# rho_e = Z_0/((4/3)*math.pi*pow(r_0,3))
# Z_0 = 4/3*math.pi*pow(r_0,3)*pow((2*theta),(3/2))/(2*pow(math.pi,2))*(I(-eta)_(1/2))

# если eta >> 1
Z_0 = 317.5 * A * pow(T[2], (3 / 2)) * math.exp(-eta[2]) / (rho)
y1 = []
r = 0.01
r1 = []
while r < r_0-0.01:
    V = Z_0 / r * (1 - 3 / 2 * r / r_0 + 1 / 2 * pow((r / r_0), 3))
    y = r * V
    y1 += y
    r1 += r
    r += 0.01


