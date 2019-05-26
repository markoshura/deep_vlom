# ПАРАМЕТРЫ ЯЧЕЙКИ

import numpy as np

from math import pi, log, exp
from Dirak_functions import integral_1_2
from Tabular_values import a_0, Na, E_h




# Безразмерный потенциал
def eta_0(T, rho):

    q = 3 * pi / 2 * 1 / 2 / theta(T)**(3 / 2) / r_0(rho, 1)**3

    if q > 10**2:
        return - (3 / 2 * q) **(2 / 3)
    if q < 10 ** (- 3):
        return log(pi ** (1 / 2) / 2 / q)
    else:
        return 0.5 * log(pi / 6) - 1.5 * log((exp((2/3 * q**2)**(1/3))))





# Заряд ядра
def z_0(T, rho, Atom_weight):
    return 317.5 * Atom_weight * T**1.5 * 2 * integral_1_2(-eta(T, rho)) / pi**0.5 / rho


# Средний радиус атомной ячейки
def r_0(rho, Atom_weight):
    #return (1 / a_0) * (3 / (4 * pi) * Atom_weight / (rho * Na))**(1/3)
    return 1.388 * (Atom_weight / rho) ** (1 / 3)


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
