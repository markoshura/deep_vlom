from math import pi, log, e
import matplotlib.pyplot as plt
from Changing_parameters import N
from Dirak_functions import integral_1_2, integral_3_2

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
from scipy import integrate

from working_progonka import X, progonka, eta

Z = [0]*(N+1)
for i in range(N+1):
    Z[i] = i / N




def S(T, rho):


    PHI = progonka(T, rho)

    # ВСПОМ. ИНТЕГРАЛ
    def S_sub_int(T, rho):

        max_i = 2
        while PHI[max_i] / X[max_i] >= 10**3:
            max_i += 1

        integr_int_1 = [X[i] for i in range(max_i + 1)]
        integr_int_2 = [X[i] for i in range(max_i, N + 1)]

        integr_func_1 = [1/3 * pi**2 * PHI[i]**0.5 * X[i]**1.5 for i in range(max_i + 1)]
        integr_func_2 = [(5/3 * integral_3_2(PHI[i]/X[i]) - PHI[i]/X[i] * integral_1_2(PHI[i]/X[i])) * X[i]**2 for i in range(max_i, N + 1)]

        integr_res_1 = integrate.trapz(integr_func_1, integr_int_1)
        integr_res_2 = integrate.trapz(integr_func_2, integr_int_2)

        return integr_res_1 + integr_res_2

    #ЭЛЕКТРОННАЯ ЭНТРОПИЯ
    def S_e(T,rho):
        const = 4*2**(1/2)*theta(T)**(3/2)*r_0(rho)**3/pi
        return const * S_sub_int(T, rho)



    ##ПОЛНАЯ ЭНТРОПИЯ

    def S(T,rho):
        return S_e(T, rho) + 3/2 * log(1836 * Atom_weight * theta(T) * volume(rho)**(2/3) / 2 / pi, e) + 5/2

    return S(T, rho)


ENTROPHY_ISOTHERM_RHO = [[], [], [], [], []]

RHO = [[], [], [], [], []]
k = -3
                                                                                                    
while k < 2:
    rho_is = 0.001
    for i in (0.001, 0.01, 0.1, 1., 10., 60.):
        while rho_is < i * 10:
            RHO[k + 3].append(rho_is)
            ENTROPHY_ISOTHERM_RHO[k+3].append(S(10**k, rho_is))

            print(rho_is, S(10**k, rho_is))
            rho_is += i
    k += 1

for i in range(5):
    plt.plot(RHO[i], ENTROPHY_ISOTHERM_RHO[i])

plt.xscale('log')
plt.yscale('log')
plt.xlabel('rho')
plt.ylabel('S')
plt.title('ENTROPHY isotherm')
plt.grid('true')
# plt.savefig('entrophyrtherm')
plt.show()


ENTROPHY_ISOHORE_T = [[], [], [], [], []]
TT = [[], [], [], [], []]
k = -3

while k < 2:
    T_is = 0.001
    for i in (0.001, 0.01, 0.1, 1., 10.):
        while T_is < i * 10:
            TT[k+3].append(T_is)
            ENTROPHY_ISOHORE_T[k+3].append(S(T_is, 10**k))
            print(T_is, S(T_is, 10**k))
            T_is += i

    k += 1

for i in range(5):
    plt.plot(TT[i], ENTROPHY_ISOHORE_T[i])

plt.xscale('log')
plt.yscale('log')
plt.xlabel("T")
plt.ylabel('S')
plt.grid('true')
plt.title('ENTROPHY isohore')
###plt.axis((0, 10000,0,950000))
### plt.savefig('entrophyhore')
plt.show()


