
import math
import numpy as np
import scipy
from scipy import integrate
from math import cos, pi
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N, T, rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2

from Atom_parameters import Atom_weight, z
from Cell import z_0, r_0, volume, theta
from State_functions import eta, rho_e
from working_progonka import RESULT3, X, PHI

const = 3 * 2 ** (1 / 2) * volume(rho) * theta(T) ** (5 / 2) / math.pi ** 2
Z = [0] * (N + 1)
for i in range(N + 1):
    Z[i] = X[i] ** (1 / 2)


# КИНЕТИЧЕСКАЯ ЭНЕРГИЯ
def E_k(T, rho):
    subfunc2 = []
    subfunc22 = []

    max_i = 0

    for i in range(1, N + 1):
        if PHI[i]/X[i] >= 10 ** 6:
            max_i = i

    subfunc1 = [0] * (max_i + 1)
    subx = [0] * (max_i + 1)
    for i in range(max_i + 1):
        subx[i] = Z[i]
    nadx = []
    for i in range(max_i + 1, N + 1):
        nadx.append(Z[i])

    for i in range(max_i + 1):
        subfunc1[i] = 4 / 5 * (PHI[i]) ** (5 / 2)

    for i in range(max_i + 1, N + 1):
        subfunc2.append(X[i]**2*integral_3_2(PHI[i]/X[i]))




    a1 = scipy.integrate.simps(subfunc1, subx)
    a2 = scipy.integrate.simps(subfunc2, nadx)
    print("a1=", a1)
    print('a2=', a2)


    return const * (a1 + a2)


# ПОТЕНЦИАЛЬНАЯ ЭНЕРГИЯ
def E_p(T, rho):
    const = 3 * 2 ** (1 / 2) * volume(rho) * theta(T) ** (5 / 2) / math.pi ** 2
    return const * integral_3_2(-eta(T, rho)) - 2 * E_k(T, rho)


# ВНУТРЕННЯЯ ЭНЕРГИЯ ЭЛЕКТРОНОВ
def E_e(T, rho):
    return E_k(T, rho) + E_p(T, rho)


# ПОЛНАЯ ЭНЕРГИЯ
def E(T, rho):
    return 2.626 * 10 ** 3 / Atom_weight * (E_e(T, rho) + 0.76874512 * z ** (7 / 3) + 3 / 2 * theta(T))


#ENERGY_ISOHORE_RHO = [[], [], [], [], [],[]]
#RHO = [[], [], [], [], [],[]]
#k = -3
#while k < 3:
#    rho_is = 0.001
#    while rho_is < 100:
#        RHO[k + 3].append(rho_is)
#        ENERGY_ISOHORE_RHO[k + 3].append(E(10 ** k, rho_is))
#        rho_is += 0.1
#    k += 1
#for i in range(5):
#    plt.plot(RHO[i], ENERGY_ISOHORE_RHO[i])
##plt.xscale('log')
##plt.yscale('log')
#plt.xlabel("rho")
#plt.ylabel('E')
#plt.title('Energy isohore')
#plt.show()
#
#ENERGY_ISOTERM_T = [[], [], [], [], []]
#TT = [[], [], [], [], []]
#k = -3
#while k < 2:
#    T_is = 0.001
#    while T_is < 10:
#        TT[k + 3].append(T_is)
#        ENERGY_ISOTERM_T[k + 3].append(E(T_is, 10 ** k))
#        T_is += 0.1
#    k +=1
#for i in range(5):
#    plt.plot(TT[i], ENERGY_ISOTERM_T[i])
##plt.xscale('log')
##plt.yscale('log')
#plt.xlabel("T")
#plt.ylabel('E')
#plt.title('Energy isoterm')
#plt.show()

print(E_e(0.0000001,0.00000001))
#должны быть равны

print("E0 = ", -0.76874512*z**(7/3))


