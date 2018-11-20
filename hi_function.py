import math
import numpy as np
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N,T,rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2, igrek_sht

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
import scipy
from scipy import integrate
from State_functions import eta, rho_e

from working_progonka import RESULT3, X


const = 4*(2*theta(T))**(1/2)/math.pi*r_0(rho)**2

#СДЕЛАЛИ СЕТКУ РАВНОМЕРНОЙ
T = [0]*(N+1)
for i in range(N+1):
    T[i] = X[i]**(1/2)

#ЗАМЕНА ПЕРЕМЕННЫХ: Y - ХИ, Z - ПРОИЗВОДНАЯ ХИ
Y = [0]*(N+1)
Y[N] = 0

Z = [0]*(N+1)
Z[N] = 0

h[i] =(2*i+1)/N**2

k1[i] = h[i]*const*(1/2*integral_minus_1_2(RESULT3[i]/T[i])*Y[i]+T[i]*igrek_sht(RESULT3[i]/T[i]))

q1[i] = h[i]*Z[i]

k2[i] = h[i]*const*(1/2*integral_minus_1_2(((RESULT3[i]+RESULT3[i-1])/2)/((T[i]+T[i-1])/2))*(Y[i]+q1[i]/2)+((T[i]+T[i-1])/2)*igrek_sht((RESULT3[i]+RESULT3[i-1])/2)/((T[i]+T[i-1])/2))

q2[i] = h[i]*(Z[i]+k1[i]/2)

k3[i] = h[i]*const*(1/2*integral_minus_1_2(((RESULT3[i]+RESULT3[i-1])/2)/((T[i]+T[i-1])/2))*(Y[i]+q2[i]/2)+((T[i]+T[i-1])/2)*igrek_sht((RESULT3[i]+RESULT3[i-1])/2)/((T[i]+T[i-1])/2))

q3[i] = h[i]*(Z[i]+k2[i]/2)

k4[i] = h[i]*const*(1/2*integral_minus_1_2(RESULT3[i-1]/T[i-1])*(Y[i]+q3[i])+T[i-1]*igrek_sht(RESULT3[i-1]/T[i-1]))

q4[i] = h*(Z[i]+k3[i])

Y[i-1] = Y[i] + (q1[1]+2*q2[i]+2*q3[i]+q4[i])/6

Z[i-1] = Z[i] + (q1[1]+2*q2[i]+2*q3[i]+q4[i])/6
