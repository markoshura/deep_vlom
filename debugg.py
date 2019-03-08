from math import log10, log, e, pi
from Changing_parameters import N
from Atom_parameters import z, Atom_weight
from Cell import theta, r_0
from Dirak_functions import integral_sht_1_2, igrek_sht
from working_progonka import progonka, X
from progonka_2N_ravnomer import progonka_2N
import matplotlib.pyplot as plt
# ПСИ = Y, ПРОИЗВОДНАЯ ПСИ = Z




PHI = progonka(0.001, 1)
PHI_2N = progonka_2N(0.001, 1)


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


const = 4 * (2 * theta(0.001)) ** 0.5 / pi * (r_0(1)) ** 2

Y[N] = -6.30651686

Z[N] = -6.30651686

for i_iter in range(N, 1, -1):
    h[i_iter] = (-2 * i_iter + 1) / N**2
    k1[i_iter] = h[i_iter] * const * (integral_sht_1_2(PHI[i_iter]/X[i_iter]) * Y[i_iter] + X[i_iter] * igrek_sht(PHI[i_iter]/X[i_iter]))
    q1[i_iter] = h[i_iter] * Z[i_iter]
    #k2[i_iter] = h[i_iter] * const * (integral_sht_1_2(((PHI[i_iter] + PHI[i_iter - 1]) / 2) / ((X[i_iter] + X[i_iter - 1]) / 2)) * (Y[i_iter] + q1[i_iter] / 2) + ((X[i_iter] + X[i_iter - 1]) / 2) * igrek_sht((PHI[i_iter] + PHI[i_iter - 1]) / 2) / ((X[i_iter] + X[i_iter - 1]) / 2))
    #k2[i_iter] = h[i_iter] * const * (integral_sht_1_2(PHI_2N[2 * i_iter - 1] / X[i_iter]) * (Y[i_iter] + q1[i_iter] / 2) + X[i_iter] * igrek_sht(PHI_2N[2 *i_iter - 1] / X[i_iter]))
    k2[i_iter] = h[i_iter] * const * (
                integral_sht_1_2((PHI[i_iter - 1] + PHI[i_iter - 1]) / 2 / X[i_iter]) * (Y[i_iter] + q1[i_iter] / 2) + X[i_iter] * igrek_sht(PHI_2N[2 * i_iter - 1] / X[i_iter]))
    q2[i_iter] = h[i_iter] * (Z[i_iter] + k1[i_iter] / 2)
    #k3[i_iter] = h[i_iter] * const * (integral_sht_1_2(((PHI[i_iter] + PHI[i_iter - 1]) / 2) / ((X[i_iter] + X[i_iter - 1]) / 2)) * (Y[i_iter] + q2[i_iter] / 2) + ((X[i_iter] + X[i_iter - 1]) / 2) * igrek_sht((PHI[i_iter] + PHI[i_iter - 1]) / 2) / ((X[i_iter] + X[i_iter - 1]) / 2))
    #k3[i_iter] = h[i_iter] * const * (integral_sht_1_2((PHI[i_iter - 1] + PHI[i_iter]) / 2 / X[i_iter]) * (Y[i_iter] + q2[i_iter] / 2) + X[i_iter] * igrek_sht(PHI_2N[2 *i_iter - 1] / X[i_iter]))
    k3[i_iter] = h[i_iter] * const * (
                integral_sht_1_2(PHI_2N[2 * i_iter - 1] / X[i_iter]) * (Y[i_iter] + q2[i_iter] / 2) + X[
            i_iter] * igrek_sht(PHI_2N[2 * i_iter - 1] / X[i_iter]))
    q3[i_iter] = h[i_iter] * (Z[i_iter] + k2[i_iter] / 2)
    k4[i_iter] = h[i_iter] * const * (integral_sht_1_2(PHI[i_iter - 1] / X[i_iter - 1]) * (Y[i_iter] + q3[i_iter]) + X[i_iter - 1] * igrek_sht(PHI[i_iter - 1] / X[i_iter - 1]))
    q4[i_iter] = h[i_iter] * (Z[i_iter] + k3[i_iter])
    Y[i_iter - 1] = Y[i_iter] + (q1[i_iter] + 2 * q2[i_iter] + 2 * q3[i_iter] + q4[i_iter]) / 6
    Z[i_iter - 1] = Z[i_iter] + (k1[i_iter] + 2 * k2[i_iter] + 2 * k3[i_iter] + k4[i_iter]) / 6

Z[0] = Z[1] + 1 / 3 * (Z[1] - Z[2])

for i in range(N+1):
    print(Y[i])





