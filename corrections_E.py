from hi_function import hi, Y, G
from working_progonka import progonka, X
from Dirak_functions import integral_1_2, igrek
from Cell import theta
from scipy import integrate
from math import pi
from Changing_parameters import N

def delta_E(T, rho, z):
    G = [(i / N) for i in range(N + 1)]

    hi(T, rho, Y, G)

    HI = Y
    DIFF_HI = G
    PHI = progonka(T, rho, 1, 1)







    def E_integrals(T, rho):
        max_i = 2
        while PHI[max_i] / X[max_i] >= 10**6:
            max_i += 1

        if max_i % 2 != 0:
            max_i -= 1


        integr_int_1 = [G[i] for i in range(max_i + 1)]
        integr_int_2 = [G[i] for i in range(max_i, N + 1)]
        #integr_int_2 = [X[i] for i in range(max_i, N + 1)]

        integr_func_1 = [4 / 3 * (11 * PHI[i]**2 * G[i] + PHI[i]**(3/2) * HI[i]) for i in range(max_i + 1)]
        integr_func_2 = [X[i] * HI[i] * integral_1_2(PHI[i] / X[i])+2 * X[i]**2 * igrek(PHI[i] / X[i]) for i in range(max_i, N + 1)]




        integr_res_1 = integrate.simps(integr_func_1, integr_int_1)


        integr_res_2 = integrate.simps(integr_func_2, integr_int_2)

        return integr_res_1 + integr_res_2

#       max_i = 0
#       for i in range(1, N + 1):
#           if PHI[i] / X[i] >= 10 ** 6:
#               max_i = i
#       subfunc1 = [0] * (max_i + 1)
#       subfunc2 = []
#       subx = [0] * (max_i + 1)
#
#for i in range(max_i+1):
#           subx[i] = X[i]
#       nadx = []
#       for i in range(max_i+1,N+1):
#           nadx.append(X[i])
#       for i in range(max_i+1):
#           subfunc1[i] = 4 / 3 * (11 * PHI[i]**2 * G[i] + PHI[i]**(3/2) * HI[i])
#       for i in range(max_i+1,N+1):
#           subfunc2.append(X[i] * HI[i] * integral_1_2(PHI[i] / X[i])+2 * X[i]**2 * igrek(PHI[i] / X[i]))
#       return integrate.trapz(subfunc1, subx) + integrate.trapz(subfunc2, nadx)

    return  E_integrals(T,rho) + 0.2690017 * z **(5/3) + z * (2 * theta(T)**(1 / 2) / 6 / pi * DIFF_HI[0])

