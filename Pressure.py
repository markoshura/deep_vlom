# ПОСТРОЕНИЕ ИЗОТЕРМ И ИЗОХОР ДАВЛЕНИЯ

import matplotlib.pyplot as plt

from math import pi, log10
from Dirak_functions import integral_3_2
from Cell import volume, theta,  rho_e
from Changing_parameters import N
from working_progonka import eta, progonka
from Tabular_values import Na




# Давление электронов
def P_e(T, rho):
    return (2 * theta(T))**2.5 / (6 * pi**2) * integral_3_2(-eta(T, rho))


# Аппроксимация давления электронного газа
def P_e_approximation(T, rho):
    return rho_e(T, rho) * (abs(theta(T))**3 + 3.36 * rho_e(T, rho) * theta(T)**1.5 + 9/125 * pi**4 * rho_e(T, rho)**2)**(1/3)


# Полное давление (ГПа)
def P(T, rho, Atom_weight, z):
    PHI = progonka (T, rho, Atom_weight, z)
    #return (P_e(T, rho) + theta(T) / volume(rho))
    return 32 / (3 * pi**3) * (2 / pi)**(2 / 3) * (2**(7/6) * 3**(2/3) * pi**(-5 / 3) * theta(T)**(1/2) * volume(rho, 1) **(2/3) * (PHI[0])**2)**(-5 / 3) * integral_3_2(PHI[N])


PRESSURE_ISOTHERM_RHO = [[], [], [], [], []]
RHO = [[], [], [], [], []]
k = -3

while k < 2:
    rho_is = 0.001
    for i in (0.001, 0.01, 0.1, 1., 10., 60.):
        while rho_is < i * 10:
            RHO[k + 3].append(rho_is)
            T_h = 10**k / 13 ** (4 / 3)

            rho_h = rho_is * Na * 11.19 * 1.4818 * 10**(-25) / 26 / 13
            PRESSURE_ISOTHERM_RHO[k+3].append(2.942 * 10**6 * P(T_h, rho_h, 1, 1) * 13**(10 / 3))

            rho_is += i
    k += 1

plt.plot(RHO[0], PRESSURE_ISOTHERM_RHO[0], label="k = - 3")
plt.plot(RHO[1], PRESSURE_ISOTHERM_RHO[1], label="k = - 2")
plt.plot(RHO[2], PRESSURE_ISOTHERM_RHO[2], label="k = -1")
plt.plot(RHO[3], PRESSURE_ISOTHERM_RHO[3], label="k = 0")
plt.plot(RHO[4], PRESSURE_ISOTHERM_RHO[4], label="k = 1")


plt.xscale('log')
plt.yscale('log')
plt.xlabel('rho, г / см^3')
plt.ylabel('P, ГПа')
plt.title('Изотермы давления алюминия')
plt.grid('true')
plt.legend()
plt.savefig('prestherm')
plt.show()


PRESSURE_ISOHORE_T = [[], [], [], [], []]
TT = [[], [], [], [], []]
k = -3

while k < 2:
    T_is = 0.001
    for i in (0.001, 0.01, 0.1, 1.):
        while T_is < i * 10:
            TT[k+3].append(T_is)
            T_h = T_is/ 13 ** (4 / 3)

            rho_h = 10**k * Na * 11.19 * 1.4818 * 10 ** (-25) / 26 / 13
            PRESSURE_ISOHORE_T[k + 3].append(2.942 * 10 ** 6 * P(T_h, rho_h, 1, 1) * 13 ** (10 / 3))
            T_is += i
    k += 1

plt.plot(TT[0], PRESSURE_ISOHORE_T[0], label="k = - 3")
plt.plot(TT[1], PRESSURE_ISOHORE_T[1], label="k = - 2")
plt.plot(TT[2], PRESSURE_ISOHORE_T[2], label="k = -1")
plt.plot(TT[3], PRESSURE_ISOHORE_T[3], label="k = 0")
plt.plot(TT[4], PRESSURE_ISOHORE_T[4], label="k = 1")

plt.xscale('log')
plt.yscale('log')
plt.xlabel("T, кэВ")
plt.ylabel('P, ГПа')
plt.grid('true')
plt.title('Изохоры давления алюминия')
plt.legend()
plt.savefig('preshore')
plt.show()

#print(P_e(Temperature_system, rho_system))

#print(log10(P(10**(3), 10**(-4), PHI = progonka(10**3, 10**(-4), 1, 1))))