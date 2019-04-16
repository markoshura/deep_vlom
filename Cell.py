# ПАРАМЕТРЫ ЯЧЕЙКИ

import numpy as np

from math import pi, log
from Dirak_functions import integral_1_2
from Tabular_values import a_0, Na, E_h
from Atom_parameters import Atom_weight, z



# Безразмерный потенциал
def eta_0(T, rho, Atom_weight):
    q = 2.795 * 10**(-3) * z * rho / (Atom_weight * T**1.5)
    try:
        res = 0.5 * log(pi / 6) - 1.5 * log((np.exp((2/3 * q**2)**(1/3), dtype=np.float64)-1))
    except OverflowError:
        print('Overflow epta!')
        res = float('inf')
    return res
    #return 0.361127101



# Заряд ядра
def z_0(T, rho, Atom_weight):
    return 317.5 * Atom_weight * T**1.5 * 2 * integral_1_2(-eta(T, rho)) / pi**0.5 / rho


# Средний радиус атомной ячейки
def r_0(rho, Atom_weight):
    return (1 / a_0) * (3 / (4 * pi) * Atom_weight / (rho * Na))**(1/3)


# Объём атомной ячейки
def volume(rho, Atom_weight):
    return 4/3 * pi * r_0(rho, Atom_weight)**3


# Характеристическая температура
def theta(T):
    return T / E_h


# Число свободных электронов на атом
def rho_e(T, rho):
    return z_0(T, rho) / volume(rho)


# Потенциал в модели ППЭ
def V(r, rho):
    return z / r * (1 - 3/2 * r / r_0(rho) + 1/2 * (r / r_0(rho))**3)

# print(z_0(0.001, 1))
