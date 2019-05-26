import matplotlib.pyplot as plt

X_1 = [1.06, 1.08, 1.13, 1.18, 1.32, 1.44, 1.57, 1.7, 1.83, 1.9, 2.1, 2.15, 2.27, 2.4, 2.6, 2.7]
Y_1 = [3, 4, 8, 14, 22, 30, 36, 47, 100, 110, 200, 205, 300, 450, 630, 800]

X_2 = [2.16, 2.3, 2.47, 2.53, 3.18]
Y_2 = [260, 400, 630, 700, 1800]

X_3 = [1.34, 1.39, 1.57, 1.9]
Y_3 = [18, 22, 30, 70]

#график ударной адиабаты кварца

plt.scatter(X_1, Y_1, color = 'green', marker='P')

plt.scatter(X_2, Y_2, color = 'red',  marker='+')
plt.scatter(X_3, Y_3, color = 'blue', marker='X')
#plt.scatter(X_2, Y_2,  marker='*')
#plt.scatter(X_2, Y_2,  marker='h')
#plt.scatter(X_2, Y_2,  marker='X')
#plt.scatter(X_2, Y_2,  marker='D')
#plt.scatter(X_2, Y_2,  marker='_')
#plt.scatter(X_2, Y_2,  marker='v')
#plt.scatter(X_2, Y_2,  marker='P')

#диапазон осей
plt.xlim(1,4)
plt.ylim(1,10000)
plt.yscale('log')
#заголовки
plt.title('Экспериментальные данные по ударному сжатию кварца')
#метки x и y
plt.xlabel('sigma')
plt.ylabel('P, GPa')
#отображение графика
plt.show()