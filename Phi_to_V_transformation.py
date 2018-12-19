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

from working_progonka import PHI, X
RESULT_V = [0]*(N+1)
R = []
F = [0]*(N+1)
V = [0]*(N+1)
for i in range(N+1):
    R.append(X[i]*r_0(rho))

mu = PHI[N]
for i in range(N+1):
    if i==0:
        F[0] = z/theta(T)
    else:
        F[i] = PHI[i]/X[i]*R[i]
for i in range(1,N+1):
    if i==N:
        V[N] = 0
    else:
        V[i] = F[i]/R[i]*theta(T) - mu
for i in range(N+1):
    if i==0:
        RESULT_V[0] = z
    else:
        RESULT_V[i] = (V[i]*R[N] + eta(T, rho)*R[i])*theta(T)
RESULT_R = [0]*(N+1)
for i in range(N+1):
    RESULT_R[i] = R[i]

fig = plt.figure()
graph1 = plt.plot(RESULT_R,RESULT_V)
plt.title('V(r)*r')

plt.grid(True)
plt.show()

