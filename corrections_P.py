#ДАВЛЕНИЕ

from working_progonka import progonka
from math import pi
from hi_function import hi
from Cell import theta
from Dirak_functions import integral_1_2, igrek
from Changing_parameters import N
from Cell import volume
def delta_P(T, rho):
    HI = hi(T, rho)
    PHI = progonka(T, rho, 1, 1)
    print("Hi_n = ", HI[N], "len = ", len(HI))
    return 8 / (3 * pi**4) * (2 / pi)**(1 / 3) * (2**(7/6) * 3**(2/3) * pi**(-5 / 3) * theta(T)**(1/2) * volume(rho, 1)**(2/3) * PHI[0]**2)**(- 4 / 3) * (HI[N] * integral_1_2(PHI[N]) + igrek(PHI[N]))


