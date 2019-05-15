from math import log
from Cell import volume
a = [0 for i in range(4)]
b = [0 for i in range(3)]
a[0] = 3785.043
a[1] = 2034.798
a[2] = -97.789
a[3] = 317.442
b[1] = -2766.977
b[2] = 797.079
V_0_c = 0.5668

def cold_curve_energy(V):
    sigma_c = V_0_c / V
    summ_1 = 0
    summ_2 = 0
    for i in range(1, 4):
        summ_1 += a[i] / i * (sigma_c **( - i / 3) - 1)
    for i in range(1 , 3):
        summ_2 += b[i] / i * (sigma_c **(i / 3) - 1)


    return a[0] * V_0_c * log(sigma_c) - 3 * V_0_c * summ_1 + 3 * V_0_c * summ_2


def cold_curve_pressure(V):
    sigma_c = V_0_c / V
    summ_1 = 0
    summ_2 = 0
    for i in range(1, 4):
        summ_1 += a[i] / i * (i / 3 * (V_0_c / V)**( -i / 3 - 1) * V_0_c / V**2)

    for i in range (1, 3):
        summ_2 += b[i] / i * (i / 3 * (V_0_c / V)**(i / 3 - 1) * (-V_0_c / V**2))

    return - (a[0] * V_0_c * V / V_0_c * (- V_0_c / V**2) - 3 * V_0_c * summ_1 + 3 * V_0_c * summ_2)


print(cold_curve_energy(volume(10, 56)))


