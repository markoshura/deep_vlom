import math
import numpy as np
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N,T,rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2, igrek

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
import scipy
from scipy import integrate

from Phi_to_V_transformation import mu

from working_progonka import PHI, X

from hi_function import HI, Diff_HI

from Internal_energy import E

#ПОПРАВКИ К МОДЕЛИ ТОМАСА-ФЕРМИ

Z = [0]*(N+1)
for i in range(N+1):
    Z[i] = X[i]**(1/2)

a1 = 0
a2 = 0
b = 0

def E_integrals(T, rho):


    max_i = 0
    for i in range(1, N + 1):
        if PHI[i]/X[i] >= 10 ** 6:
            max_i = i
    subfunc1 = [0]*(max_i+1)
    subfunc2 = []
    subx =[0]*(max_i+1)
    for i in range(max_i+1):
        subx[i] = X[i]
    nadx = []
    for i in range(max_i+1,N+1):
        nadx.append(X[i])
    for i in range(max_i+1):
        subfunc1[i] = 4/3*(11*PHI[i]**2*Z[i]+PHI[i]**(3/2)*HI[i])
    for i in range(max_i+1,N+1):
        subfunc2.append(X[i]*HI[i]*integral_1_2(PHI[i]/X[i])+2*X[i]**2*igrek(PHI[i]/X[i]))
    return scipy.integrate.trapz(subfunc1, subx) + scipy.integrate.trapz(subfunc2, nadx)


def S_integrals(T,rho):
    max_i = 0
    for i in range(1, N + 1):
        if PHI[i] / X[i] >= 10 ** 6:
            max_i = i
    subfunc1 = [0] * (max_i + 1)
    subfunc2 = []
    subx = [0] * (max_i + 1)
    for i in range(max_i + 1):
        subx[i] = X[i]
    nadx = []
    for i in range(max_i + 1, N + 1):
        nadx.append(X[i])
    for i in range(max_i + 1):
        subfunc1[i] = 4 / 3 * (22 * PHI[i] ** 2 * Z[i] + PHI[i] ** (3 / 2) * HI[i])
    for i in range(max_i + 1, N + 1):
        subfunc2.append(X[i]*HI[i]*integral_1_2(PHI[i]/X[i])+4*Z[i]**2*igrek(PHI[i]/X[i]))
    return scipy.integrate.trapz(subfunc1, subx) + scipy.integrate.trapz(subfunc2, nadx)




#ДАВЛЕНИЕ
def delta_P(T,rho):
    return theta(T)**2/(3*math.pi**3)*(HI[N]*integral_1_2(PHI[N]) + igrek(PHI[N]))



def delta_E(T,rho):
    return 2*theta(T)**2/(3*math.pi**2)*r_0(rho)**3*E_integrals(T,rho) + 0.2690017*z**(5/3)

def delta_S(T,rho):
    return 2*theta(T)/(3*math.pi**2)*r_0(rho)**2*S_integrals(T,rho) + 2**(1/2)*z*Diff_HI[0]/(6*math.pi*theta(T)**(1/2))


def delta_mu(T, rho):
    return (2*theta(T))**(1/2)/(6*math.pi)*(1/2*integral_minus_1_2(PHI[N]) + HI[N])


DELTA_E_ISOTERM = [[],[],[],[],[]]

RHO_ISOTERM = [[],[],[],[],[]]




k = -3
while k<2:

    rho_is = 0.001
    while rho_is<100:
        DELTA_E_ISOTERM[k+3].append(delta_E(10**k, rho_is))
        RHO_ISOTERM[k+3].append(rho_is)
        rho_is += 0.1
    k+=1

for i in range(5):
    plt.plot(RHO_ISOTERM[i], DELTA_E_ISOTERM[i])


plt.xscale('log')
plt.yscale('log')
plt.xlabel("rho")
plt.ylabel('E')
plt.title('Delta Energy isoterm')
plt.show()








