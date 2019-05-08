#Интерполируем тепловой вклад


НА ВХОД T и rho


billinear_approximation

#Интерполируем холодную кривую

n = 100

cold_curve = [[0 for value in range(n + 1)] for j in range(2)]

def interp_cold_curve(rho):
    result = 0
    for i in range(n + 1):
        for j in range(1, n + 1):
            if i == j:
                result += cold_curve[1][i]
            else:
                result += cold_curve[1][i] * (rho - cold_curve[0][j]) / (cold_curve[0][i] - cold_curve[0][j])



#Прибавляем ионный вклад
