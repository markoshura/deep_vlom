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


Z = [0]*(N+1)

h = [0]*(N+1)
k1 = [0]*(N+1)
q1 = [0]*(N+1)
k2 = [0]*(N+1)
q2 = [0]*(N+1)
k3 = [0]*(N+1)
q3 = [0]*(N+1)
k4 = [0]*(N+1)
q4 = [0]*(N+1)

def hi_function(sigma):

    Y[N] = sigma

    Z[N] = sigma



    for i in range(N,1,-1):


        h[i] =(2*i+1)/N**2

        k1[i] = h[i]*const*(1/2*integral_minus_1_2(RESULT3[i]/T[i])*Y[i]+T[i]*igrek_sht(RESULT3[i]/T[i]))

        q1[i] = h[i]*Z[i]

        k2[i] = h[i]*const*(1/2*integral_minus_1_2(((RESULT3[i]+RESULT3[i-1])/2)/((T[i]+T[i-1])/2))*(Y[i]+q1[i]/2)+((T[i]+T[i-1])/2)*igrek_sht((RESULT3[i]+RESULT3[i-1])/2)/((T[i]+T[i-1])/2))

        q2[i] = h[i]*(Z[i]+k1[i]/2)

        k3[i] = h[i]*const*(1/2*integral_minus_1_2(((RESULT3[i]+RESULT3[i-1])/2)/((T[i]+T[i-1])/2))*(Y[i]+q2[i]/2)+((T[i]+T[i-1])/2)*igrek_sht((RESULT3[i]+RESULT3[i-1])/2)/((T[i]+T[i-1])/2))

        q3[i] = h[i]*(Z[i]+k2[i]/2)

        k4[i] = h[i]*const*(1/2*integral_minus_1_2(RESULT3[i-1]/T[i-1])*(Y[i]+q3[i])+T[i-1]*igrek_sht(RESULT3[i-1]/T[i-1]))

        q4[i] = h[i]*(Z[i]+k3[i])

        Y[i-1] = Y[i] + (q1[i]+2*q2[i]+2*q3[i]+q4[i])/6

        Z[i-1] = Z[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6

    Z[0] = Z[1]+1/3*(Z[1]-Z[2])
    i = 1
    h[i] = (2 * i + 1) / N ** 2

    k1[i] = h[i] * const * (1 / 2 * integral_minus_1_2(RESULT3[i] / T[i]) * Y[i] + T[i] * igrek_sht(RESULT3[i] / T[i]))

    q1[i] = h[i] * Z[i]

    k2[i] = h[i] * const * (1 / 2 * integral_minus_1_2(((RESULT3[i] + RESULT3[i - 1]) / 2) / ((T[i] + T[i - 1]) / 2)) * (
                Y[i] + q1[i] / 2) + ((T[i] + T[i - 1]) / 2) * igrek_sht((RESULT3[i] + RESULT3[i - 1]) / 2) / (
                                        (T[i] + T[i - 1]) / 2))

    q2[i] = h[i] * (Z[i] + k1[i] / 2)

    k3[i] = h[i] * const * (1 / 2 * integral_minus_1_2(((RESULT3[i] + RESULT3[i - 1]) / 2) / ((T[i] + T[i - 1]) / 2)) * (
                Y[i] + q2[i] / 2) + ((T[i] + T[i - 1]) / 2) * igrek_sht((RESULT3[i] + RESULT3[i - 1]) / 2) / (
                                        (T[i] + T[i - 1]) / 2))

    q3[i] = h[i] * (Z[i] + k2[i] / 2)



    q4[i] = h[i] * (Z[i] + k3[i])

    Y[i - 1] = Y[i] + (q1[i] + 2 * q2[i] + 2 * q3[i] + q4[i]) / 6

    return Y[0]
def sigma_2(sigma_0,sigma_1):
    return sigma_1 - (sigma_1-sigma_0)*hi_function(sigma_1)/(hi_function(sigma_1) - hi_function(sigma_0))


sigma_2 = sigma_2(1,-1)


hi_function(sigma_2)
print("Y = ", Y)
print("Z = ", Z)

plt.plot(T,Z)
plt.title("Производная пси (Z)")
plt.savefig("graph1")
plt.show()



plt.plot(T,Y)
plt.title("Пси (Y)")
plt.savefig("graph2")
plt.show()


HI = [0]*(N+1)
Diff_HI = [0]*(N+1)

for i in range(N+1):
    HI[i] = Y[i]
    Diff_HI[i] = Z[i]

print(sigma_2)



