#ДАВЛЕНИЕ
def delta_P(T,rho):
    HI = hi_function(T, rho)
    PHI = progonka(T, rho)
    return theta(T)**2/(3*math.pi**3)*(HI[N]*integral_1_2(PHI[N]) + igrek(PHI[N]))