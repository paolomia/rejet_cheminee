import pandas as pd
from calculs_condensation_maille import *




def p_H_f(H, rho_L, rho_m):
    r"""
    p_H = H \cdot g \cdot (\rho_L - \rho_m)
    """
    # (31) pg 34
    # Tirage théorique du au tirage naturel
    g = 9.81 * ureg.m / ureg.s ** 2
    return H * g * (rho_L - rho_m)

def main():
    # import csv
    df = pd.read_csv('calculs_condensation_maille.csv', sep=';')

    # Prenons la valeur moyenne de rho_r
    rho_m_mean = df['rho_m'].mean() * ureg.kg / ureg.m**3
    print(f"rho_m moyen: {rho_m_mean}")

    T_L = Q_(-15, 'degC')
    R_L = 288 * ureg.joule / (ureg.kg * ureg.kelvin)
    z = 41 * ureg.meter
    p_L = p_L_f(T_L, z=z, R_L=R_L)
    rho_L = p_L / (R_L * T_L.to('kelvin'))

    # Tirage naturel théorique
    H = 7.5 * ureg.meter
    p_H = p_H_f(H, rho_L, rho_m_mean)
    print(f"Tirage naturel théorique p_H: {p_H.to('pascal')}")

    # Difference de pression du au changement de vitesse
    w_1 = df['w_m'].loc[0] * ureg.m / ureg.s
    rho_1 = df['rho_m'].loc[0] * ureg.kg / ureg.m**3
    w_2 = df['w_m'].loc[len(df)-1] * ureg.m / ureg.s
    rho_2 = df['rho_m'].loc[len(df)-1] * ureg.kg / ureg.m**3

    P_G = 0.5 * (rho_2 * w_2**2 - rho_1 * w_1**2)
    print(f"Différence de pression due au changement de vitesse P_G: {P_G.to('pascal')}")

    # Perte de charge du au frottement
    L_v = 0.2 * ureg.m
    L_fum = 6.3 * ureg.m + 1.2 * ureg.m
    L_tot = L_v + L_fum

    D_h = 0.2 * ureg.m
    psi_mean = df['psi'].mean()
    w_mean = df['w_m'].mean() * ureg.m / ureg.s

    P_R_frottement = (psi_mean * L_tot/D_h * (rho_m_mean * (w_mean**2) / 2)).to('pascal')
    print(f"Perte de charge due au frottement P_R_frottement: {P_R_frottement.to('pascal')}")
    S_E = 1.2 # pag 28, 5.7.8 dans notre cas
    S_EG = S_E if P_G.magnitude >= 0 else 1

    # Perte de charge du au raccordement
    # pg 106 avec m_dot_1 = 0
    Zeta_racc = 1.2

    # Nous avons une seule changement de direction/raccord, entre le circuit de raccordement et le circuit de fumée
    sum_Zeta_i = Zeta_racc

    P_R_raccord = (sum_Zeta_i * (rho_m_mean * (w_mean**2) / 2)).to('pascal')
    P_E = P_R_frottement + P_R_raccord
    P_R = S_E * P_E + S_EG * P_G
    print(f"Perte de charge totale P_R: {P_R.to('pascal')}")

    # Pression de la vitesse du vent
    P_L = 25 * ureg.pascal # pg 36, 5.10.4

    P_ZO = P_R - p_H + P_L
    P_ZO_min = P_R - p_H
    print(f'{P_ZO=}, {P_ZO_min=}')



if __name__ == "__main__":
    main()