#ПАРАМЕТРЫ ЯЧЕЙКИ

import math
import numpy as np


from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2
from Tabular_values import a_0, Na, E_h
from Atom_parameters import Atom_weight,z



#Безразмерный потенциал
def eta(T,rho):
    q = 2.795*10**(-3)*z*rho/(Atom_weight*T**(3/2))
    return 1/2*math.log(math.pi/6,math.e)-3/2*math.log((math.exp((2/3*q**2)**(1/3))-1),math.e)

#Заряд ядра
def z_0(T, rho):
    return 317.5 * Atom_weight * T**(3/2)*2/(np.pi**(1/2)) * integral_1_2(-eta(T,rho))/(rho)

#Средний радиус атомной ячейки
def r_0(rho):
    return 1/a_0*(3/(4*math.pi)*Atom_weight/(rho*Na))**(1/3)

#Объём атомной ячейки
def volume(rho):
    return 4/3*math.pi*(r_0(rho))**3


#Характеристическая температура
def theta(T):
    return T/E_h


