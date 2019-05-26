from math import log
from Cell import volume
a = [0 for i in range(4)]
b = [0 for i in range(3)]
a[0] = 6923.207
a[1] = - 4772.762
a[2] = 1263.499
a[3] = 50.038
b[1] = - 5061.082
b[2] = 1597.1
V_0_c = 0.3618

def cold_curve_energy(V):
    sigma_c = V_0_c / V
    summ_1 = 0
    summ_2 = 0
    for i in range(1, 4):
        summ_1 += a[i] / i * (sigma_c ** (- i / 3) - 1)
    for i in range(1 , 3):
        summ_2 += b[i] / i * (sigma_c ** (i / 3) - 1)

    print(summ_1)
    print(summ_2)
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


print("energy = ", cold_curve_energy(volume(1, 1)))


