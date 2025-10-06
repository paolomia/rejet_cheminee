import pint
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
import math

def m_dot(Q_F, f_m1, f_m2, sigma_CO2):
    r"""
    \dot{m} = \left(\frac{f_{m1}}{\sigma(\text{CO}_2)} + f_{m2}\right) \cdot Q_F
    """
    res = (f_m1 / sigma_CO2 + f_m2) * Q_F.to("kW").magnitude
    return Q_(res/1000, "kg/s")


def sigma_CO2_f(Q_N, f_x1, f_x2, f_x3):
    r"""
    \sigma(\text{CO}_2) =
    \begin{cases}
    \frac{f_{x1}}{1 - f_{x2} \cdot \log_{10} Q_N} & \text{pour } Q_N \leq 100 \, \text{kW} \\
    f_{x3} & \text{pour } Q_N > 100 \, \text{kW}
    \end{cases}
    """
    # pg 101 B2
    if Q_N <= 100 * ureg.kW:
        return f_x1 / (1 - f_x2 * math.log10(Q_N.magnitude))
    else:
        return f_x3

# Q est la puissance utile
# Q_N est la puissance utile nominale
# Q_F est le debit calorifique

# Nous devons valider les exigencer pour la puissance utile nomintale de l'appareil à combustion déclarée par le fabricant
# ous sinon la puissance utile nominale et la plus petite valeur de la plage de puissance utile declarée par le fabricant
# Pour nous le discours est un peu compliqué, car nous avons un bruleur + foyer
# d'après les essais de Frank, le rendement du foyer est de 86.5 %
# la puissance du bruleur est de 48/79 - 195 kW mais je n'ai pas de nominale
# je vais donc prendre le nominal le plus bas (condition la plus défavorable) -> 48 kW

Q = 140 * ureg.kW  # puissance utile
Q_F = Q / 0.86async  # debit calorifique
# et donc on va prendre le nominal le plus bas
Q_N = Q

f_m1 = 3.735
f_m2 = 0.0535
# Annexe B3, valeurs pour gaz et tirage forcé
f_x1 = 8.6
f_x2 = 0.078
f_x3 = 10.2

sigma_CO2 = sigma_CO2_f(Q_F, f_x1, f_x2, f_x3)
# sigma_CO2 = 10.0  # valeur mesurée par Frank
print(f"Fraction massique de CO2 dans les fumées: {sigma_CO2:.1f} %")
mdot = m_dot(Q_F, f_m1, f_m2, sigma_CO2)
print(f"Débit massique des fumées: {mdot:.5f}")