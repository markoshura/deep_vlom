def eta(i):
    q = 2.795*10**(-3)*z*rho/(A*T[i]**(3/2))
    print(q)
    return 1/2*math.log(math.pi/6,math.e)-3/2*math.log((math.exp((2/3*q**(2))**(1/3))-1),math.e)

T = [0.01, 0.1, 1, 10]
z = 79
a_0 =  0.52917721067*10**(-8)
Na = 6.022* 10**(23)
r_0 = 1/a_0*(3/(4*math.pi)*A/(rho*Na))**(1/3)
A = 196.96657
rho = 1

eta(0)