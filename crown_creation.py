from bilinear_interpolation import bilinear_interpolation

T = 100
rho = 0.01
n = 100
cold_curve = [[0 for value in range(n + 1)] for j in range(2)]
#Интерполируем тепловой вклад

bilinear_interpolation(T, rho)

#Из функции делаем ключ-значение

def func_cold_curve():
    cold_curve = [[0 for value in range(n + 1)] for j in range(2)]




    return cold_curve


#Интерполируем холодную кривую

def interp_cold_curve(rho):
    cold_curve = func_cold_curve()
    result = 0
    for i in range(n + 1):
        for j in range(1, n + 1):
            if i == j:
                result += cold_curve[1][i]
            else:
                result += cold_curve[1][i] * (rho - cold_curve[0][j]) / (cold_curve[0][i] - cold_curve[0][j])



#Прибавляем ионный вклад

def state_function(T, rho):

    return bilinear_interpolation(T, rho) + interp_cold_curve(rho) + ion(T, rho)