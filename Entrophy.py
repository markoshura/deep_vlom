import math
import numpy as np
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N,T,rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
import scipy
from scipy import integrate
from State_functions import eta, rho_e

from working_progonka import  X, PHI

Z = [0]*(N+1)
for i in range(N+1):
    Z[i] = X[i]**(1/2)

#ЭЛЕКТРОННАЯ ЭНТРОПИЯ
def S_e(T,rho):
    const =  4*2**(1/2)*theta(T)**(3/2)*r_0(rho)**3/math.pi

    subfunc2 = []

    max_i = 0
    for i in range(1, N + 1):
        if PHI[i] / (theta(T) * X[i]) >= 10 ** 6:
            max_i = i
    subfunc1 = [0]*(max_i+1)
    subx =[0]*(max_i+1)
    for i in range(max_i+1):
        subx[i] = X[i]
    nadx = []
    for i in range(max_i+1,N+1):
        nadx.append(X[i])
    for i in range(max_i+1):
        subfunc1[i] = 1/3*math.pi**2*(PHI[i])**(1/2)*X[i]**(3/2)
    for i in range(max_i+1,N+1):
        subfunc2.append(((5/3*integral_3_2(PHI[i]/X[i])) - (PHI[i]/X[i])*integral_1_2(PHI[i]/X[i]))*X[i]**2)
    a1 = scipy.integrate.trapz(subfunc1,subx)
    a2 = scipy.integrate.trapz(subfunc2,nadx)

    return const*(a1+a2)


#ПОЛНАЯ ЭНТРОПИЯ

def S(T,rho):
    return 0.9648*10**2/Atom_weight*(S_e(T,rho)+3/2*math.log(1836*Atom_weight*theta(T)*volume(rho)**(2/3)/2/math.pi,math.e)) + 5/2

#ENTROPHY_ISOTERM_RHO =[[],[],[],[],[]]
#RHO =[[],[],[],[],[]]
#k = -3
#while k<2:
#    rho_is = 0.001
#    while rho_is<100:
#        RHO[k+3].append(rho_is)
#        ENTROPHY_ISOTERM_RHO[k+3].append(S(10**k,rho_is))
#        rho_is+=1
#    k+=1
#for i in range(5):
#    plt.plot(RHO[i], ENTROPHY_ISOTERM_RHO[i])
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel("rho")
#plt.ylabel('S')
#plt.title('Entrophy изохора')
#plt.show()
#
#ENTROPHY_ISOTERM_T =[[],[],[],[],[]]
#TT =[[],[],[],[],[]]
#k = -3
#while k<2:
#    T_is = 0.001
#    while T_is<10:
#        TT[k+3].append(T_is)
#        ENTROPHY_ISOTERM_T[k+3].append(S(T_is,10**k))
#        T_is+=1
#    k+=1
#for i in range(5):
#    plt.plot(TT[i], ENTROPHY_ISOTERM_T[i])
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel("T")
#plt.ylabel('S')
#plt.title('Entrophy isoterm')
#plt.show()