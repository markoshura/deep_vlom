# ИЗОТЕРМЫ И ИЗОХОРЫ ВНУТРЕННЕЙ ЭНЕРГИИ

import matplotlib.pyplot as plt

from math import pi
from scipy import integrate
from Atom_parameters import E_0, Atom_weight, z 
from Changing_parameters import N
from Dirak_functions import integral_3_2
from Cell import volume, theta, eta  ,r_0
from working_progonka import X, PHI
#from sveryaem_phi import kakaxa_phi as PHI, kakaxa_x as X
from Pressure import P_e
#from Changing_parameters import Temperature_system, rho_system

Z = [(i / N) for i in range(N + 1)]


# ВСПОМ. ИНТЕГРАЛ
def E_sub_int(T, rho):
    max_i = 1
    while PHI[max_i] / X[max_i] >= 10**3:
        max_i += 1

    integr_int_1 = [Z[i] for i in range(max_i + 1)]
    integr_int_2 = [Z[i] for i in range(max_i, N + 1)]
    #integr_int_2 = [X[i] for i in range(max_i + 1, N + 1)]

    integr_func_1 = [4/5 * PHI[i]**2.5 for i in range(max_i + 1)]
    integr_func_2 = [2 * Z[i]**5 * integral_3_2(PHI[i] / Z[i]**2) for i in range(max_i, N + 1)]
    #integr_func_2 = [X[i] ** 2 * integral_3_2(PHI[i] / X[i]) for i in range(max_i + 1, N + 1)]

    integr_res_1 = integrate.simps(integr_func_1, integr_int_1)
    integr_res_2 = integrate.simps(integr_func_2, integr_int_2)

    # Очень нечестно
    # integr_int = [X[i] for i in range(1, N + 1)]
    # integr_func = [X[i]**2 * integral_3_2(PHI[i] / X[i]) for i in range(1, N + 1)]
    # integr_res = integrate.cumtrapz(integr_func, integr_int)[-1]
   # k1 = [i for i in range(max_i + 1)]
   # k2 = [i for i in range(max_i, N + 1)]
   #
   # plt.plot(k1,integr_func_1)
   # plt.plot(k2, integr_func_2)
   # plt.axis((0, 2048,0,950000))
   # plt.show()



    return integr_res_1 + integr_res_2


const = 2**(0.5) / (pi**2) * (theta(T)**2.5) * (4/3) * pi * r_0(rho)**3
# КИНЕТИЧЕСКАЯ ЭНЕРГИЯ
def E_k(T, rho):
    #return (3 * 2**0.5 / pi**2) * volume(rho) * theta(T)**2.5 * E_sub_int(T, rho)
    return 3*const * E_sub_int(T,rho)
#print('E_k_old=', (3 * 2**0.5 / pi**2) * volume(rho) * theta(T)**2.5 * E_sub_int(T, rho))
#print ("E_k_new=", 3*const*E_sub_int(T,rho))
#print("E_p_old=", (2 * 2**0.5 / pi**2) * volume(rho) * theta(T)**2.5 * (integral_3_2(-eta(T, rho)) - 3 * E_sub_int(T, rho)) )
#print("E_p_new=",2*const*(integral_3_2(-eta(T,rho)) - 3*E_sub_int(T,rho) )  )

# ПОТЕНЦИАЛЬНАЯ ЭНЕРГИЯ
def E_p(T, rho):
    #return (2 * 2**0.5 / pi**2) * volume(rho) * theta(T)**2.5 * (integral_3_2(-eta(T, rho)) - 3 * E_sub_int(T, rho))
    return 2 * const * (integral_3_2(-eta(T,rho)) - 3 * E_sub_int(T,rho))

# Проверяем virial theorem
print(2 * E_k(T, rho) + E_p(T, rho))
print(3 * P_e(T, rho) * volume(rho))



# ВНУТРЕННЯЯ ЭНЕРГИЯ ЭЛЕКТРОНОВ
def E_e(T, rho):
    return E_k(T, rho) + E_p(T, rho)



# ПОЛНАЯ ЭНЕРГИЯ
def E(T, rho):
    #return (2.626 * 10**3 / Atom_weight) * (E_e(T, rho) - E_0) ##3/2 * theta(T)
    return const * ( 2 * integral_3_2(-eta(T,rho)) - 3*E_sub_int(T,rho)) #+ 0.76874512421364*z**(7/3)
#print("E_old=",((E_e(T, rho) - E_0) + 3/2 * theta(T) ) )
#print("E_old_1=",((E_e(T, rho) - E_0) )  )
#print("E_new=", const*(2*integral_3_2(-eta(T,rho))-3*E_sub_int(T,rho))+0.76874512421364*z**(7/3))
#print("E_new_1=", E_p(T,rho)+E_k(T,rho)+0.76874512421364*z**(7/3))

ENERGY_ISOTHERM_RHO = [[], [], [], [], []]
RHO = [[], [], [], [], []]
k = -3

while k < 2:
    rho_is = 0.001
    for i in (0.001, 0.01, 0.1, 1., 10., 60.):
        while rho_is < i * 10:
            RHO[k + 3].append(rho_is)
            ENERGY_ISOTHERM_RHO[k+3].append(E(10**k, rho_is))
            rho_is += i
    k += 1

for i in range(5):
    plt.plot(RHO[i], ENERGY_ISOTHERM_RHO[i])

plt.xscale('log')
plt.yscale('log')
plt.xlabel('rho')
plt.ylabel('E')
plt.title('ENERGY isotherm')
plt.grid('true')
# plt.savefig('enertherm')
plt.show()


ENERGY_ISOHORE_T = [[], [], [], [], []]
TT = [[], [], [], [], []]
k = -3

while k < 2:
    T_is = 0.001
    for i in (0.001, 0.01, 0.1, 1.,10.):
        while T_is < i * 10:
            TT[k+3].append(T_is)
            ENERGY_ISOHORE_T[k+3].append(E(T_is, 10**k))
            #print(T_is, E(T_is,rho))
            T_is += i

    k += 1

for i in range(5):
    plt.plot(TT[i], ENERGY_ISOHORE_T[i])
E_sub_int(T, rho)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("T")
plt.ylabel('E')
plt.grid('true')
plt.title('ENERGY isohore')
#plt.axis((0, 10000,0,950000))
# plt.savefig('enerhore')
plt.show()

print('E=', E(10,rho))



