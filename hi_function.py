from math import pi
import matplotlib.pyplot as plt
from Changing_parameters import N
from Dirak_functions import integral_minus_1_2, igrek_sht, integral_sht_1_2


from Cell import z_0, r_0, volume, theta

from working_progonka import progonka, X


#ЗАМЕНА ПЕРЕМЕННЫХ: Y - ХИ, Z - ПРОИЗВОДНАЯ ХИ
def hi(T, rho):


    PHI = progonka(T, rho, 1, 1)

    const = 2**(7/6) * 3**(2/3) * pi**(-5 / 3) * theta(T)**(1/2) * volume(rho, 1) ** (2/3)



    def hi_function(sigma):
        Y = [0] * (N + 1)
        Z = [0] * (N + 1)

        h = [0] * (N + 1)
        k1 = [0] * (N + 1)
        q1 = [0] * (N + 1)
        k2 = [0] * (N + 1)
        q2 = [0] * (N + 1)
        k3 = [0] * (N + 1)
        q3 = [0] * (N + 1)
        k4 = [0] * (N + 1)
        q4 = [0] * (N + 1)

        Y[N] = sigma

        Z[N] = sigma
        for i in range(N, 1, -1):


            h[i] = (-2 * i + 1) / N**2
            #print("h = ", h[i])
            #print("i = ", i, "N = ", N, "sigma = ", sigma)


            k1[i] = h[i] * const * (integral_sht_1_2(PHI[i] / X[i]) * Y[i] + X[i] * igrek_sht(PHI[i] / X[i]))

            q1[i] = h[i] * Z[i]

            k2[i] = h[i] * const * (integral_sht_1_2(((PHI[i] + PHI[i-1]) / 2)/((X[i] + X[i-1]) / 2))*(Y[i] + q1[i] / 2) + ((X[i] + X[i-1]) / 2)*igrek_sht((PHI[i] + PHI[i-1]) / 2)/((X[i] + X[i-1])/2))

            q2[i] = h[i] * (Z[i] + k1[i] / 2)

            k3[i] = h[i] * const * (integral_sht_1_2(((PHI[i] + PHI[i-1]) / 2)/((X[i] + X[i-1]) / 2))*(Y[i] + q2[i] / 2) + ((X[i] + X[i-1]) / 2) * igrek_sht((PHI[i] + PHI[i-1]) / 2) / ((X[i] + X[i-1]) / 2))

            q3[i] = h[i] * (Z[i] + k2[i] / 2)

            k4[i] = h[i] * const * (integral_sht_1_2(PHI[i-1] / X[i-1]) * (Y[i] + q3[i]) + X[i-1] * igrek_sht(PHI[i-1] / X[i-1]))

            q4[i] = h[i] * (Z[i] + k3[i])

            Y[i-1] = Y[i] + (q1[i] + 2 * q2[i] + 2 * q3[i] + q4[i]) / 6

            Z[i-1] = Z[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6

        Z[0] = Z[1] + 1 / 3 * (Z[1] - Z[2])
        i = 1
        h[i] = - (- 2 * i + 1) / N ** 2

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

        Y[i - 1] = Y[i] + (q1[i] + 2 * q2[i] + 2 * q3[i] + q4[i]) / 6

        #plt.plot(Y)
        #plt.show()
        return Y

    def sigma_2(sigma_0, sigma_1):
        HI1 = hi_function(sigma_1)
        HI0 = hi_function(sigma_0)
        print("sigma = ", sigma_1 - (sigma_1 - sigma_0) * HI1[0] / (HI1[0] - HI0[0]))
        return sigma_1 - (sigma_1 - sigma_0) * HI1[0] / (HI1[0] - HI0[0])


    sigma_22 = sigma_2(100, -100)

    return(hi_function(sigma_22))
print(hi(0.001, 1))


