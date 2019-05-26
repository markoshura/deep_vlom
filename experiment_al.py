import matplotlib.pyplot as plt

X_1 = [4.4, 4.45, 4.47, 4.48, 4.55, 3.95, 3.94, 3.6]
Y_1 = [10000, 12000, 16000, 15500, 30000, 4000, 3700, 3200]
X_2 = [1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2]
Y_2 = [4, 5, 5, 8, 7, 8, 9, 10, 11, 12, 13, 14.5, 16, 17, 18, 20]

X_3 = [1.22, 1.23, 1.24, 1.27, 1.28, 1.32, 1.37, 1.4, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5]
Y_3 = [23, 25, 27, 30, 33, 35, 45, 45, 50, 55, 65, 65, 70, 74, 77, 76, 78, 74]

X_4 = [1.52, 1.54, 1.56, 1.58, 1.6, 4.4, 4.45, 4, 3.43, 3.35, 3.3, 3.25, 3.07, 2.7]
Y_4 = [80, 84, 88, 92, 96, 10000, 12000, 3000, 1800, 2500, 2000, 1700, 1000, 1000]

X_5 = [2.2, 2.22, 2.25]
Y_5 = [800, 900, 900]

X_6 = [2.3, 2.4, 2.6, 2.6, 2.65, 2.64]
Y_6 = [1000, 1000, 800, 700, 800, 1200]

X_7 = [2.4, 2.45]
Y_7 = [950, 1000]

X_8 = [1.61, 1.62, 1.63, 1.64, 1.65, 1.66, 1.67, 1.68, 1.69, 1.7]
Y_8 = [110, 110, 115, 120, 120, 130, 135, 140, 145, 145]

X_9 = [1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8]
Y_9 = [155, 165, 175, 170, 175, 180, 190, 196, 200, 215]

X_10 = [1.81, 1.82, 1.83, 1.93, 1.94, 1.95]
Y_10 = [225, 240, 270, 350, 400, 410]

X_11 = [1.84, 1.85, 1.86, 1.87]
Y_11 = [200, 280, 220, 300]

X_12 = [1.88, 1.89, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95]
Y_12 = [350, 250, 400, 400, 410, 350, 400, 410]

X_13 = [1.96, 1.97, 1.98, 2.1, 2.15]
Y_13 = [440, 480, 500, 540, 700]

X_14 = [1.99, 2]
Y_14 = [520, 520]
#график ударной адиабаты алюминия

plt.scatter(X_1, Y_1,  marker='s', facecolors='none', edgecolors='b')

plt.scatter(X_2, Y_2,  marker='+', facecolors='none', edgecolors='b')
plt.scatter(X_3, Y_3,  marker='o', facecolors='none', edgecolors='b')
plt.scatter(X_4, Y_4,  marker='*', facecolors='none', edgecolors='b')
plt.scatter(X_5, Y_5,  marker='h', facecolors='none', edgecolors='b')
plt.scatter(X_6, Y_6,  marker='X', facecolors='none', edgecolors='b')
plt.scatter(X_7, Y_7,  marker='D', facecolors='none', edgecolors='b')
plt.scatter(X_8, Y_8,  marker='_')
plt.scatter(X_9, Y_9,  marker='v', facecolors='none', edgecolors='b')
plt.scatter(X_10, Y_10,  marker='D', facecolors='none', edgecolors='b')
plt.scatter(X_11, Y_11,  marker='>', facecolors='none', edgecolors='b')
plt.scatter(X_12, Y_12,  marker='<', facecolors='none', edgecolors='b')
plt.scatter(X_13, Y_13,  marker='x', facecolors='none', edgecolors='b')
plt.scatter(X_14, Y_14,  marker='H', facecolors='none', edgecolors='b')

#диапазон осей
plt.xlim(1,5)
plt.ylim(1,10000000)
plt.yscale('log')
#заголовки
plt.title('Экспериментальные данные по ударному сжатию алюминия')
#метки x и y
plt.xlabel('sigma')
plt.ylabel('P, GPa')
plt.savefig('experiment_al')
#отображение графика
plt.show()