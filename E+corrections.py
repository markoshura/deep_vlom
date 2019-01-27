import math
import numpy as np
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N,T,rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2, igrek

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
import scipy
from scipy import integrate

from Phi_to_V_transformation import mu

from working_progonka import PHI, X, RESULT3

from hi_function import HI, Diff_HI

from Internal_energy import E
from corrections import delta_E
DELTA_E_ISOTERM = [[],[],[],[],[]]

RHO_ISOTERM = [[],[],[],[],[]]

DELTA_E_plus_E_ISOTERM = [[],[],[],[],[]]


k = -3
while k<2:

    rho_is = 0.001
    while rho_is<100:
        DELTA_E_plus_E_ISOTERM.append(delta_E(10**k, rho_is)+ E(10**k, rho_is))

        DELTA_E_ISOTERM[k+3].append(delta_E(10**k, rho_is))
        RHO_ISOTERM[k+3].append(rho_is)
        rho_is += 0.1
    k+=1

for i in range(5):
    plt.plot(RHO_ISOTERM[i], DELTA_E_ISOTERM[i])
    plt.plot(RHO_ISOTERM[i], DELTA_E_plus_E_ISOTERM[i])

plt.xscale('log')
plt.yscale('log')
plt.xlabel("rho")
plt.ylabel('E')
plt.title('Delta Energy isoterm')
plt.show()