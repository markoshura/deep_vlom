import math
import numpy as np
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N,T,rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
from State_functions import eta, rho_e


# ВСЁ ДЛЯ ЛИНЕАРИЗАЦИИ И ПРОГОНКИ С ИТЕРАЦИЯМИ

s = 5

s_current = 0

PHI_S = [[0] * (N + 1)] * 5
X = [0] * (N + 1)
A = [0] * (N + 1)
B = [0] * (N + 1)
C = [0] * (N + 1)
D = [0] * (N + 1)
ITERATIONS = [[0] * (N + 1)] * 5
RESULT0 = [0] * (N + 1)
RESULT1 = [0] * (N + 1)
RESULT2 = [0] * (N + 1)
RESULT3 = [0] * (N + 1)
RESULT4 = [0] * (N + 1)
ALPHA = [0] * (N + 1)
BETA = [0] * (N + 1)
Y = [0] * (N + 1)
# константа а в уравнении
const = 4 * (2 * theta(T)) ** 2 / math.pi * (r_0(rho)) ** 2


# СЕТКА
def x():
    for i in range(N + 1):
        X[i] = (i / N) ** 2


def phi(s=s_current, T=T, rho=rho):
    for i in range(N + 1):

        if i == 0:
            PHI_S[s_current][0] = z / (theta(T) * r_0(rho))
        else:
            PHI_S[s_current][i] = z / (theta(T) * r_0(rho)) * (1 - 3 / 2 * X[i] + 1 / 2 * X[i] ** 3) - eta(T,rho) * X[i]


# ШАГ СЕТКИ

def h(N):
    x()
    return 1 - X[N - 1]


# КОЭФФИЦИЕНТЫ
def a():
    x()
    for i in range(1, N):
        A[i] = 1 / (X[i] - X[i - 1])


def c():
    for i in range(1, N):
        C[i] = 1 / (X[i + 1] - X[i])


def b(s=s_current):
    for i in range(1, N):
        B[i] = -A[i] - C[i] - const / 4 * (X[i + 1] - X[i - 1]) * integral_minus_1_2(PHI_S[s_current][i] / X[i])


def d(s=s_current):
    for i in range(1, N):
        D[i] = const * (X[i + 1] - X[i - 1]) * (
                    (X[i] / 2 * integral_1_2(PHI_S[s_current][i]/X[i])) - PHI_S[s_current][i] / 4 * integral_minus_1_2(
                PHI_S[s_current][i] / X[i]))


def alpha(s=s_current):
    for i in range(N - 1, -1, -1):
        if i == N - 1:
            ALPHA[i] = 1 / (1 - h(N) + h(N) ** 2 / 4 * const * integral_minus_1_2(PHI_S[s_current][N]))
        else:
            ALPHA[i] = -A[i + 1] / (B[i + 1] + C[i + 1] * ALPHA[i + 1])


def beta(s=s_current):
    for i in range(N - 1, -1, -1):
        if i == N - 1:
            BETA[i] = -h(N) ** 2 / 2 * const * (
                        integral_1_2(PHI_S[s_current][N]) - PHI_S[s_current][i] / 2 * integral_minus_1_2(
                    PHI_S[s_current][N])) * ALPHA[N - 1]
        else:
            BETA[i] = (D[i + 1] - C[i + 1] * BETA[i + 1]) / (B[i + 1] + C[i + 1] * ALPHA[i + 1])


# ИСКОМОЕ УРАВНЕНИЕ
def y(T=T):
    for i in range(N + 1):
        if i == 0:
            Y[i] = z / (theta(T) * r_0(rho))
        else:
            Y[i] = ALPHA[i - 1] * Y[i - 1] + BETA[i - 1]



h(N)

x()
# print("X=",X)
phi()
# print(phi(0))
# print("PHI_S_0 =", PHI_S[s_current])
a()
# print("A=", A)
c()
# print("C=",C)
b()
# print("B=", B)
d()
# print("D=", D)
alpha()
# print("ALPHA=", ALPHA)
beta()
# print("BETA=", BETA)
y()
# print("Y=", Y)
for i in range(N + 1):
    RESULT0[i] = PHI_S[s_current][i] * theta(T)
s_current = s_current + 1

for i in range(N + 1):
    PHI_S[s_current][i] = Y[i]
    # print(PHI_S[s_current][i])

a()
# print("A=", A)
c()
# print("C=",C)
b()
# print("B=", B)
d()
# print("D=", D)
alpha()
# print("ALPHA=", ALPHA)
beta()
# print("BETA=", BETA)
y()
# print("Y=", Y)
for i in range(N + 1):
    RESULT1[i] = Y[i] * theta(T)
s_current = s_current + 1

for i in range(N + 1):
    PHI_S[s_current][i] = Y[i]
    # print(PHI_S[s_current][i])

a()
# print("A=", A)
c()
# print("C=",C)
b()
# print("B=", B)
d()
# print("D=", D)
alpha()
# print("ALPHA=", ALPHA)
beta()
# print("BETA=", BETA)
y()
# print("Y=", Y)
for i in range(N + 1):
    RESULT2[i] = Y[i] * theta(T)
s_current = s_current + 1

for i in range(N + 1):
    PHI_S[s_current][i] = Y[i]
    # print(PHI_S[s_current][i])

a()
# print("A=", A)
c()
# print("C=",C)
b()
# print("B=", B)
d()
# print("D=", D)
alpha()
# print("ALPHA=", ALPHA)
beta()
# print("BETA=", BETA)
y()
# print("Y=", Y)
for i in range(N + 1):
    RESULT3[i] = Y[i] * theta(T)
s_current = s_current + 1

for i in range(N + 1):
    PHI_S[s_current][i] = Y[i]
    # print(PHI_S[s_current][i])


PHI = [0]*(N+1)
for i in range(N+1):
    PHI[i] = RESULT3[i]/theta(T)



#fig = plt.figure()
#graph1 = plt.plot(X, RESULT0)
#graph2 = plt.plot(X, RESULT1)
#graph3 = plt.plot(X, RESULT2)
#graph4 = plt.plot(X, RESULT3)

#plt.grid(True)
#plt.show()

#plt.plot(X,RESULT3)
#plt.show()