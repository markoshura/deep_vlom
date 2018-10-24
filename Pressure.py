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

Y_P = [[],[],[],[],[],[]]
X_rho = [[],[],[],[],[],[]]

i = -3
while i<3:
    T_current = T**i
    Y_P[0].append(math.log(P(T_current,rho),10))
    X_rho[0].append(i)
    i+=0.2

fig = plt.figure()

graph1 = plt.plot(X_rho[0], Y_P[0])
plt.title('Plot: '+str(1)+' '+'T= '+ str(10))


plt.grid(True)
plt.show()
