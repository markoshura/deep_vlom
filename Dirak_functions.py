# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ

import math

# Функция Ферми-Дирака I_1/2
def integral_1_2(x):
    sqx = x ** (1 / 2)
    xsqx = x * sqx
    # Апроксимация функции дирака трёхчленной формулой при k=1/2
    if (x >= 176868.709):
        return xsqx / 1.5 + 0.82246703342411321823620758332301 / sqx + 0.71027462212293443818237742585514 / (x * xsqx);
    elif (x < -17.9):
        return 0.88622692545275801364908374167057 * math.exp(x)
    else:
        if (x > 30.0):
            i0 = x
        else:
            i0 = math.log(1.0 + math.exp(x), math.e)

        return (0.8862269254528) * i0 * ((1.0 + 1.178 * i0 + (0.1812102675) * (i0 ** 3.0)) ** (1.0 / 6.0))

    # Функция Ферми-Дирака I_3/2


def integral_3_2(x):
    return 3 / 10 * integral_1_2(x) * (125 + 60 * integral_1_2(x) + 18 * (integral_1_2(x)) ** 2) ** (1 / 3)


def integral_minus_1_2(x):
    if (x > 500000.0):
        return 2 * x ** (1 / 2)
    else:
        if (x < -17.9):
            return 2 * 0.88622692545275801364908374167057 * math.exp(x)
        else:
            if (x > 30.0):
                i0 = x
            else:
                i0 = math.log(1.0 + math.exp(x))

            return 2 * (0.8862269254528) * i0 * ((1.0 + 1.614 * i0 + (0.4844730731296) * (i0 ** 3.0)) ** (-1.0 / 6.0))

def igrek_sht(x):
    dx = max(0.14, 0.001 * abs(x))
    return (-1/2*integral_minus_1_2(x + dx) + 1/2*integral_minus_1_2(x - dx) + 168 * 1/2*integral_minus_1_2(x + dx / 2.0) - 168 * 1/2*integral_minus_1_2(x - dx / 2.0) - 5376 * 1/2*integral_minus_1_2(x + dx / 4.0) + 5376 * 1/2*integral_minus_1_2(x - dx / 4.0) + 32768 * 1/2*integral_minus_1_2(x + dx / 8.0) - 32768 * 1/2*integral_minus_1_2(x - dx / 8.0)) / (5670 * dx)


