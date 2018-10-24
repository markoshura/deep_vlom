import math
import matplotlib.pyplot as plt
from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
from Changing_parameters import N,T,rho



#ВАЖНЫЕ ФУНКЦИИ СОСТОЯНИЯ

#Потенциал в модели ППЭ
def V(r, T, rho):
    return z/r * (1 - 3 / 2 * r / r_0(rho) + 1 / 2 * pow((r / (r_0(rho))), 3))


#Число свободных электронов на атом
def rho_e(T, rho):
    return z_0(T, rho)/volume(rho)

#Безразмерный потенциал
def eta(T,rho):
    q = 2.795*10**(-3)*z*rho/(Atom_weight*T**(3/2))
    return 1/2*math.log(math.pi/6,math.e)-3/2*math.log((math.exp((2/3*q**2)**(1/3))-1),math.e)
Y_V_r = [[],[],[],[],[],[]]
X_r = [[],[],[],[],[],[]]
for i in range(0,6):
    T = 0.01 + i*2
    r = 0.001
    r_0(rho)
    while r < r_0(rho):
        eta(T,rho)
        z_0(T,rho)
        V(r,T,rho)

        y = V(r,T,rho)*r
        Y_V_r[i].append(y)
        X_r[i].append(r)
        r += 0.01
#fig = plt.figure()
#for i in range(0,6):
#    graph1 = plt.plot(X_r[i], Y_V_r[i])
#    plt.title('Plot: ' + str(i+1) + ' ' + 'T= ' +  str(0.01 + i*2))
#plt.grid(True)
#plt.show()
