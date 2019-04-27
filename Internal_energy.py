# ИЗОТЕРМЫ И ИЗОХОРЫ ВНУТРЕННЕЙ ЭНЕРГИИ

import matplotlib.pyplot as plt

from math import pi
from scipy import integrate
from Atom_parameters import E_0
from Changing_parameters import N , Temperature_system, rho_system
from Dirak_functions import integral_3_2
from Cell import volume, theta, r_0
from working_progonka import X, progonka, eta
#from sveryaem_phi import excel_phi as PHI, excel_x as X
#from Pressure import P_e
#from Changing_parameters import Temperature_system, rho_system
from Tabular_values import a_0, Na

Z = [(i / N) for i in range(N + 1)]




def Energy(T, rho, z, Atom_weight, PHI):

    #PHI = progonka(T, rho, 1, 1)

    # ВСПОМ. ИНТЕГРАЛ
    def E_sub_int(T, rho):

        max_i = 2
        while (PHI[max_i] / X[max_i] >= 10**6) and (max_i < N):
            max_i += 1

        if max_i % 2 != 0:
            max_i -= 1


        integr_int_1 = [Z[i] for i in range(max_i + 1)]
        integr_int_2 = [Z[i] for i in range(max_i, N + 1)]
        #integr_int_2 = [X[i] for i in range(max_i, N + 1)]

        integr_func_1 = [4/5 * PHI[i]**2.5 for i in range(max_i + 1)]
        integr_func_2 = [2 * Z[i]**5 * integral_3_2(PHI[i] / Z[i]**2) for i in range(max_i, N + 1)]
        #integr_func_2 = [X[i] ** 2 * integral_3_2(PHI[i] / X[i]) for i in range(max_i, N + 1)]


        #i = 0
        #integr_res_1 = 0
        #while i <= (max_i - 2):
        #
        #    integr_res_1 += 1/6 * (integr_func_1[i] + 4 * integr_func_1[i+1] + integr_func_1[i+2]) * (integr_int_1[i+2] - integr_int_1[i])
        #    i += 2
        #
        #integr_res_2 = 0
        #i = max_i
        #while i <= N - 2:
        #
        #    integr_res_2 += 1/6 * (integr_func_2[i - max_i] + 4 * integr_func_2[i+1 - max_i] + integr_func_2[i+2-max_i]) * (integr_int_2[i+2-max_i] - integr_int_2[i-max_i])
        #    i += 2



        integr_res_1 = integrate.simps(integr_func_1, integr_int_1)


        integr_res_2 = integrate.simps(integr_func_2, integr_int_2)

        return integr_res_1 + integr_res_2


    const = 2**(0.5) / (pi**2) * (theta(T)**2.5) * (4/3) * pi * r_0(rho, 1)**3

    # КИНЕТИЧЕСКАЯ ЭНЕРГИЯ
   #def E_k(T, rho):
   #    #return (3 * 2**0.5 / pi**2) * volume(rho) * theta(T)**2.5 * E_sub_int(T, rho)
   #
   #    return 3 * const * E_sub_int(T, rho)
   ## ПОТЕНЦИАЛЬНАЯ ЭНЕРГИЯ
   #def E_p(T, rho):
   #    #return (2 * 2**0.5 / pi**2) * volume(rho) * theta(T)**2.5 * (integral_3_2(-eta(T, rho)) - 3 * E_sub_int(T, rho))
   #    return 2 * const * (integral_3_2(-eta(T, rho, 1, 1)) - 3 * E_sub_int(T, rho))



# Проверяем virial theorem
#print(2 * E_k(Temperature_system, rho_system) + E_p(Temperature_system, rho_system))
#print(3 * P_e(Temperature_system, rho_system) * volume(rho_system))



    ## ВНУТРЕННЯЯ ЭНЕРГИЯ ЭЛЕКТРОНОВ
    #def E_e(T, rho):
    #   return E_k(T, rho) + E_p(T, rho)
    #

    # ПОЛНАЯ ЭНЕРГИЯ
    def E(T, rho):

        #return E_e(T, rho) - E_0 + 3/2 * theta(T)
        
        return const * (2 * integral_3_2(-eta(T, rho, 1, 1)) - 3*E_sub_int(T, rho)) + 0.76874512421364*1**(7/3)

    return E(T , rho)
    #print("E_sub = ", E_sub_int(T, rho))
    #print ("E_k = ", E_k(T, rho))
    #print("E_p = ", E_p(T, rho))
    #print("E_e = ", E_e(T, rho))
    #print("E = ", E(T, rho))

#print(Energy(36.4, 100, 13, 26))
#z = 13
#A = 26
#T = 36.4
#rho = 100
#T_h = T / z**(4 / 3)
#rho_h = rho * Na * 11.19 * 1.4818 * 10**(-25) / A / z
#print("T_h = ", T / z**(4 / 3))
#print("rho_h = ", rho * Na * 11.19 * 1.4818 * 10**(-25) / A / z )
#print("E = ", Energy(T_h, rho_h, 1, 1) * z**(7 / 3))
#ENERGY_ISOTHERM_RHO = [[], [], [], [], []]
#RHO = [[], [], [], [], []]
#k = -3
#
#while k < 2:
#    rho_is = 0.001
#    for i in (0.001, 0.01, 0.1, 1., 10.):
#        while rho_is < i * 10:
#            RHO[k + 3].append(rho_is)
#            ENERGY_ISOTHERM_RHO[k+3].append(Energy(10**k, rho_is))
#            rho_is += i
#    k += 1
#
#for i in range(5):
#    plt.plot(RHO[i], ENERGY_ISOTHERM_RHO[i])
#
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('rho')
#plt.ylabel('E')
#plt.title('ENERGY isotherm')
#plt.grid('true')
## plt.savefig('enertherm')
#plt.show()
#
#
#ENERGY_ISOHORE_T = [[], [], [], [], []]
#TT = [[], [], [], [], []]
#k = -3
#
#while k < 2:
#    T_is = 0.001
#    for i in (0.001, 0.01, 0.1, 1., 10.):
#        while T_is < i * 10:
#            TT[k+3].append(T_is)
#            ENERGY_ISOHORE_T[k+3].append(Energy(T_is, 10**k))
#            T_is += i
#
#    k += 1
#
#for i in range(5):
#    plt.plot(TT[i], ENERGY_ISOHORE_T[i])
#
#
###E_sub_int(T, rho)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel("T")
#plt.ylabel('E')
#plt.grid('true')
#plt.title('ENERGY isohore')
####plt.axis((0, 10000,0,950000))
#### plt.savefig('enerhore')
#plt.show()