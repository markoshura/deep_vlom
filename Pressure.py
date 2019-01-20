import math
import numpy as np
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N,T,rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
from State_functions import eta, rho_e

#Давление электронов
def P_e(T,rho):
    return (2*theta(T))**(5/2)/(6*math.pi**2)*integral_3_2(-eta(T,rho))

#Аппроксимация давления электронного газа
def P_e_approximation(T, rho):
    return rho_e(T, rho)*((theta(T))**3+3.36*rho_e(T, rho)*(theta(T))**(3/2)+9/125*math.pi**4*(rho_e(T, rho))**2)**(1/3)

#Полное давление (ГПа)
def P(T, rho):
    return 2.942*10**4*(P_e(T,rho)+theta(T)/volume(rho))


PRESSURE_ISOTERM_RHO =[[],[],[],[],[]]
RHO =[[],[],[],[],[]]
k = -3
while k<2:
    rho_is = 0.001
    while rho_is<100:
        RHO[k+3].append(rho_is)
        PRESSURE_ISOTERM_RHO[k+3].append(P(10**k,rho_is))
        rho_is+=1
    k+=1
for i in range(5):
    plt.plot(RHO[i], PRESSURE_ISOTERM_RHO[i])
plt.xscale('log')
plt.yscale('log')
plt.xlabel("rho")
plt.ylabel('P')
plt.title('PRESSURE isohore')
plt.grid('true')
plt.savefig('preshore')
plt.show()

PRESSURE_ISOTERM_T =[[],[],[],[],[]]
TT =[[],[],[],[],[]]
k = -3
while k<2:
    T_is = 0.001
    while T_is<10:
        TT[k+3].append(T_is)
        PRESSURE_ISOTERM_T[k+3].append(P(T_is,10**k))
        T_is+=1
    k+=1
for i in range(5):
    plt.plot(TT[i], PRESSURE_ISOTERM_T[i])
plt.xscale('log')
plt.yscale('log')
plt.xlabel("T")
plt.ylabel('P')
plt.grid('true')
plt.title('PRESSURE isoterm')
plt.savefig('presterm')
plt.show()