
import math
import numpy as np
from scipy import integrate
from math import cos, pi
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N,T,rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
from State_functions import eta, rho_e
from scipy import integrate
def phi_0(x):
    return z/(theta(T)*r_0(rho))*(1-3/2*x+1/2*x**3)-eta(T,rho)*x
XX=[]
YY = []
x=0.00000000001
while x<1:

    XX.append(x)
    YY.append(phi_0(x)/x)
    x+=0.0001
fig = plt.figure()
graph1 = plt.plot(XX, YY)
plt.grid(True)
plt.show()

#КИНЕТИЧЕСКАЯ ЭНЕРГИЯ
#def E_k():
#    func = lambda x: x**2*integral_3_2()
#    answer = integrate.quad(subfunc, 0, 1)
#    print(answer)
#    print(z/(theta(T)*r_0(rho)))
#E_k()












