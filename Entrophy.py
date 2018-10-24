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
def phi_0(x):
    return z/(theta(T)*r_0(rho))*(1-3/2*x+1/2*x**3)-eta(T,rho)*x
def podint(x):
    print(-5/3*integral_3_2(phi_0(x)/x))
    return -5/3*integral_3_2(phi_0(x)/x)


XX=[]
YY = []
x=0.00000000001
while x<1:

    XX.append(x)
    YY.append(podint(x))
    x+=0.0001
fig = plt.figure()
graph1 = plt.plot(XX, YY)
plt.grid(True)
plt.show()