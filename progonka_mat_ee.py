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

const = 21.034*((T/z)**(4/3))**(1/2)/(rho/(Atom_weight*z))**(2/3)
#СЕТКА
y = [0]*3
for i in range(3):
    y[i] = [0]*(N+1)
U =[0]*(N+1)
X = [0]*(N+1)
A = [0]*N
C = [0]*N
B = [0]*N
D = [0]*N
h = 1/N
for i in range(N+1):
    U[i] = i*h

PHI = [0]*3
for i in range(3):
    PHI[i] = [0]*(N+1)
PHI[0][0] = z/(theta(T)*r_0(rho))
for i in range(1,N+1):
    PHI[0][i] = z/(theta(T)*r_0(rho))*(1-3/2*U[i]**2+1/2*U[i]**2**3)-theta(T)*U[i]**2

for i in range(1,N):
    A[i] = 2*U[i]+h

for i in range(1,N):
    C[i] = 2*U[i] - h
for i in range(1,N):
    B[i] = -4*U[i]*(1+const*h**2*U[i]**2*integral_minus_1_2(PHI[0][i]/U[i]**2))
for i in range(1,N):
    D[i] = 4*const*h**2*U[i]**3*(2*U[i]**2*integral_1_2(PHI[0][i]/U[i]**2)-PHI[0][i]*integral_minus_1_2(PHI[0][i]/U[i]**2))

ALPHA = [0]*N
BETA = [0]*N
ALPHA[N-1] = 1/(1-2*h+h**2*(1+const*integral_minus_1_2(PHI[0][N])))
BETA[N-1] = -(const*h**2*(2*integral_1_2(PHI[0][N]) - PHI[0][N]*integral_minus_1_2(PHI[0][N])))*ALPHA[N-1]
for i in range(N-2,0,-1):
    ALPHA[i] = -A[i+1]/(B[i+1]+C[i+1]*ALPHA[i+1])
    BETA[i] = (D[i+1]-C[i+1]*BETA[i+1])/(B[i+1]+C[i+1]*ALPHA[i+1])

for i in range(N+1):
    y[0][i] = theta(T)*PHI[0][i]

for i in range(N):
    y[1][i] = ALPHA[i]*y[0][i]+BETA[i]
for i in range(N+1):
    X[i] = U[i]**2

print(U)
print(A)
print(B)
print(C)
print(D)
print(ALPHA)
print(BETA)
#graph1 = plt.plot(X,y[0])
#graph2 = plt.plot(X, y[1])
#plt.grid()
#plt.show()
