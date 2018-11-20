import scipy
from scipy import integrate
from working_progonka import X
from Changing_parameters import N,T,rho
func = []
x = []
for i in range(N+1):
    x.append(i/N)
for i in range(N+1):
    func.append(x[i]**2)
print(scipy.integrate.trapz(func,x))
