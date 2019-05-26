from math import exp, log
from scipy import arctan
from Cell import volume
V_0 = 0.3687
T_a = 0.00089443 * 10**4 * 1.1604521 # KK??????? Спросить Максима. 1 эв = 1.1604521 * 10**4 К

teta_0 = 0.2 * 10**4 * 1.1604521 # kK спросить
gamma_0 = 1.95
B = 0.5
D = 0.35667
R = 8.31
def ionic_contribution_energy(T, V):
    #F - T * dF / dT
    sigma = V_0 / V
    f1 = - exp(-teta(V) / T - (T_a * sigma**(2 / 3) / T)**(1 / 2))
    f2 = -teta(V) / T
    f3 = - (T_a * sigma**(2 / 3) / T)**(1 / 2)
    f4 =  exp((gamma_0 - 2 / 3) * (B**2 + D**2) / B * arctan(B * log(sigma) / (B**2 + D *(log(sigma) + D))))
    f5 = B * log(sigma) / (B**2 + D *(log(sigma) + D))
    f6 = B * log(sigma)
    f7 = B**2 + D *(log(sigma) + D)
    f8 = (T_a * sigma**(2 / 3) / T)**(1 / 2)
    df8_T = - T_a * sigma**(2 / 3) / T**2
    df3_T = - 1 / 2 * f8**(- 1 / 2) * df8_T
    df2_T = teta(V) / T**2
    df1_T = f1 * (df2_T + df3_T)
    dF_T = 3 * R * log(1 + f1) + 3 * R * T / (1 + f1) * df1_T
    return 3 * R * T * log(1 + f1) - T * dF_T

def ionic_contribution_pressure(T, V):
    # - dF / dV

    sigma = V_0 / V
    f1 = - exp(-teta(V) / T - (T_a * sigma**(2 / 3) / T)**(1 / 2))
    f2 = -teta(V) / T
    f3 = - (T_a * sigma**(2 / 3) / T)**(1 / 2)
    f4 = exp((gamma_0 - 2 / 3) * (B**2 + D**2) / B * arctan(B * log(sigma) / (B**2 + D *(log(sigma) + D))))
    f5 = B * log(sigma) / (B**2 + D *(log(sigma) + D))
    f6 = B * log(sigma)
    f7 = B**2 + D * (log(sigma) + D)
    f8 = (T_a * sigma**(2 / 3) / T)**(1 / 2)
    df8_V = T_a / T * 2 / 3 * sigma**(- 1 / 3) * (-V_0 / V**2)
    df3_V = - 1 / 2 * (f8)**(- 1 / 2) * df8_V
    df5_V = (B * V / V_0 * (-V_0 / V**2) * f7 - f6 * D * V / V_0 * (-V_0 / V**2)) / (f7**2)
    darctanf5_V = 1 / (1 + f5**2) * df5_V
    df4_V = f4 * (gamma_0 - 2 / 3) * (B**2 + D**2) / B * darctanf5_V
    df2_V = -teta_0 / T * (2 / 3 * (V_0 / V) ** (-1 / 3) * (-V_0 / V ** 2) * f4 + sigma ** (2 / 3) * df4_V)
    df1_V = f1 * (df2_V + df3_V)

    return - 3 * R * T / (1 + f1) * df1_V



def teta(V):
    sigma = V_0 / V
    return teta_0 * sigma**(2 / 3) * exp((gamma_0 - 2 / 3) * (B**2 + D**2) / B * arctan(B * log(sigma) / (B**2 + D * (log(sigma) + D))))

print(ionic_contribution_pressure(1, volume(1, 1)))