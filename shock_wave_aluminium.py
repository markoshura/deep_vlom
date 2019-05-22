from Pressure import P
from Internal_energy import Energy
from numpy import dot
from numpy import transpose
from working_progonka import progonka
rho_00 = 2.712
E_0 = Energy(293 * 10**(-3) / (1.160 * 10**4), rho_00, 13, 26)
P_0 = P(293 * 10**3 / (1.160 * 10**4), rho_00, 26, 13)



def f1(T, rho):
    return P(T, rho, 26, 13) - P_0


def f2(T, rho):
    return Energy(T, rho, 13, 26) - E_0 - 1 / 2 * (P(T, rho, 26, 13) + P_0) * (1 / rho_00 - 1 / rho)

def newton_solver(P, T_pribl, rho_pribl):
    delta = [[0 for i in range(1)] for j in range(2)]
    delta[0][0] = 10**(- 8)
    delta[1][0] = 10**(- 8)
    epsilon = 10**(- 8)
    mod_delta = 100


    W = [[0 for i in range(2)] for j in range(2)]
    F = [[0 for i in range(1)] for j in range(2)]
    priblizhenie = [[0 for i in range(1)] for j in range(2)]

    priblizhenie[0][0] = T_pribl
    priblizhenie[1][0] = rho_pribl

    next_priblizhenie = [[0 for i in range(1)] for j in range(2)]
    while mod_delta > epsilon:
        print("mod_delta = ", mod_delta)



        W[0][0] = (f1(T_pribl + 0.1 * T_pribl, rho_pribl) - f1(T_pribl - 0.1 * T_pribl, rho_pribl) )/ (2 * T_pribl * 0.1)
        W[0][1] = (f1(T_pribl, rho_pribl + 0.1 * rho_pribl) - f1(T_pribl, rho_pribl - 0.1 * rho_pribl)) / (2 * rho_pribl * 0.1)
        W[1][0] = (f2(T_pribl + 0.1 * T_pribl, rho_pribl) - f2(T_pribl - 0.1 * T_pribl, rho_pribl) )/ (2 * T_pribl * 0.1)
        W[1][1] = (f2(T_pribl, rho_pribl + 0.1 * rho_pribl) - f2(T_pribl, rho_pribl - 0.1 * rho_pribl)) / (2 * rho_pribl * 0.1)

        F[0][0] = f1(T_pribl, rho_pribl)
        F[1][0] = f2(T_pribl, rho_pribl)


        W_1 = transpose(W)

        next_priblizhenie[0][0] = priblizhenie[0][0] - (W_1 @ F)[0][0]
        next_priblizhenie[1][0] = priblizhenie[1][0] - (W_1 @ F)[1][0]

        delta[0][0] = next_priblizhenie[0][0] - priblizhenie[0][0]
        delta[1][0] = next_priblizhenie[1][0] - priblizhenie[1][0]

        print("next_pribl = ", next_priblizhenie)
        print("pribl = ", priblizhenie)

        priblizhenie[0][0] = next_priblizhenie[0][0]
        priblizhenie[1][0] = next_priblizhenie[1][0]

        T_pribl = priblizhenie[0][0]
        rho_pribl = priblizhenie[1][0]
        print("T = ", T_pribl)
        print("rho = ", rho_pribl)

        mod_delta = (delta[0][0]**2 + delta[1][0]**2)**(1 / 2)

    return next_priblizhenie


print(newton_solver(100, 293 * 10**(-3) / (1.160 * 10**4) * 2 , 6))

#delta = [[0 for i in range(1)] for j in range(2)]
#delta[0][0] = 10**(- 8)
#delta[1][0] = 10**(- 8)
#print(delta)
#








