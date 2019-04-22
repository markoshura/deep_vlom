import math
import matplotlib.pyplot as plt
from Changing_parameters import N
from Dirak_functions import integral_minus_1_2, igrek_sht

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta

from working_progonka import progonka, X

Y = [0]*(N+1)
Z = [0]*(N+1)
#ЗАМЕНА ПЕРЕМЕННЫХ: Y - ХИ, Z - ПРОИЗВОДНАЯ ХИ
def hi(T, rho, Y, Z):

    PHI = progonka(T, rho, 1, 1)
    const = 4*(2*theta(X))**(1/2)/math.pi*r_0(rho)**2


    h = [0]*(N+1)
    k1 = [0]*(N+1)
    q1 = [0]*(N+1)
    k2 = [0]*(N+1)
    q2 = [0]*(N+1)
    k3 = [0]*(N+1)
    q3 = [0]*(N+1)
    k4 = [0]*(N+1)
    q4 = [0]*(N+1)

    def hi_function(sigma, Y, Z):

        Y[N] = sigma

        Z[N] = sigma
        for i in range(N,1,-1):


            h[i] =(2*i+1)/N**2

            k1[i] = h[i]*const*(1/2*integral_minus_1_2(PHI[i]/X[i])*Y[i]+X[i]*igrek_sht(PHI[i]/X[i]))

            q1[i] = h[i]*Z[i]

            k2[i] = h[i]*const*(1/2*integral_minus_1_2(((PHI[i]+PHI[i-1])/2)/((X[i]+X[i-1])/2))*(Y[i]+q1[i]/2)+((X[i]+X[i-1])/2)*igrek_sht((PHI[i]+PHI[i-1])/2)/((X[i]+X[i-1])/2))

            q2[i] = h[i]*(Z[i]+k1[i]/2)

            k3[i] = h[i]*const*(1/2*integral_minus_1_2(((PHI[i]+PHI[i-1])/2)/((X[i]+X[i-1])/2))*(Y[i]+q2[i]/2)+((X[i]+X[i-1])/2)*igrek_sht((PHI[i]+PHI[i-1])/2)/((X[i]+X[i-1])/2))

            q3[i] = h[i]*(Z[i]+k2[i]/2)

            k4[i] = h[i]*const*(1/2*integral_minus_1_2(PHI[i-1]/X[i-1])*(Y[i]+q3[i])+X[i-1]*igrek_sht(PHI[i-1]/X[i-1]))

            q4[i] = h[i]*(Z[i]+k3[i])

            Y[i-1] = Y[i] + (-q1[i]-2*q2[i]-2*q3[i]-q4[i])/6

            Z[i-1] = Z[i] + (-k1[i]-2*k2[i]-2*k3[i]-k4[i])/6

        Z[0] = Z[1]+1/3*(Z[1]-Z[2])
        i = 1
        h[i] = (2 * i + 1) / N ** 2

        k1[i] = h[i] * const * (1 / 2 * integral_minus_1_2(PHI[i] / X[i]) * Y[i] + X[i] * igrek_sht(PHI[i] / X[i]))

        q1[i] = h[i] * Z[i]

        k2[i] = h[i] * const * (1 / 2 * integral_minus_1_2(((PHI[i] + PHI[i - 1]) / 2) / ((X[i] + X[i - 1]) / 2)) * (
                    Y[i] + q1[i] / 2) + ((X[i] + X[i - 1]) / 2) * igrek_sht((PHI[i] + PHI[i - 1]) / 2) / (
                                            (X[i] + X[i - 1]) / 2))

        q2[i] = h[i] * (Z[i] + k1[i] / 2)

        k3[i] = h[i] * const * (1 / 2 * integral_minus_1_2(((PHI[i] + PHI[i - 1]) / 2) / ((X[i] + X[i - 1]) / 2)) * (
                    Y[i] + q2[i] / 2) + ((X[i] + X[i - 1]) / 2) * igrek_sht((PHI[i] + PHI[i - 1]) / 2) / (
                                            (X[i] + X[i - 1]) / 2))

        q3[i] = h[i] * (Z[i] + k2[i] / 2)



        q4[i] = h[i] * (Z[i] + k3[i])

        Y[i - 1] = Y[i] + (-q1[i] - 2 * q2[i] - 2 * q3[i] - q4[i]) / 6

    def sigma_2(sigma_0,sigma_1):
        return sigma_1 - (sigma_1-sigma_0)*hi_function(sigma_1, Y, Z)/(hi_function(sigma_1, Y, Z) - hi_function(sigma_0, Y, Z))


    sigma_2 = sigma_2(1, -1)




    return hi_function(sigma_2, Y, Z)


