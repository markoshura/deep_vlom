
from Pressure import P
from Tabular_values import Na
import matplotlib.pyplot as plt
ENERGY_ISOTHERM_RHO = [[], [], [], [], []]
RHO = [[], [], [], [], []]
k = 0

while k < 3:
    rho_is = 0.001
    for i in (0.001, 0.01, 0.1, 1., 10.):
        while rho_is < i * 10:
            RHO[k].append(rho_is)
            T_h = 3.98 / (10**3) * 10**k / 13 **(4 / 3)
            rho_h = rho_is * Na * 11.19 * 1.4818 * 10**(-25) / 26 / 13
            ENERGY_ISOTHERM_RHO[k].append(2.942 * 10**6 * P(T_h, rho_h, 1, 1) * 13**(10 / 3))
            rho_is += i
    k += 1

plt.plot(RHO[0], ENERGY_ISOTHERM_RHO[0], label="3.98 эв")
plt.plot(RHO[1], ENERGY_ISOTHERM_RHO[1], label="39.8 эв")
plt.plot(RHO[2], ENERGY_ISOTHERM_RHO[2], label="398 эв")
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.001, 100)
plt.ylim(0, 10000000)
plt.xlabel('rho, г / см^3')
plt.ylabel('P, ГПа')
plt.title('Электронное давление алюминия по различным моделям')
plt.grid('true')
plt.legend()
plt.savefig('compare_pressure')
plt.show()



