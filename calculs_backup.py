from report import *
import math
import json

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


def eta_w(Q_N):
    r"""
    \eta_w =
    \begin{cases}
    85.0 + 1.0 \cdot \log_{10} Q_N & \text{pour } Q_N \leq 1000 \,\text{kW} \\
    88.0 & \text{pour } Q_N > 1000 \,\text{kW}
    \end{cases}
    """
    # pg 101 B2
    if Q_N <= 1000 * ureg.kW:
        return 85.0 + 1.0 * math.log10(Q_N.magnitude)
    else:
        return 88.0


def sigma_H2O_f(f_w, sigma_CO2):
    r"""
    \sigma(H_2O) =
    \frac{100}{1 + \frac{f_w}{\sigma(CO_2)}} + 1.1
    """
    # pg 99
    return 100 / (1 + f_w / sigma_CO2) + 1.1


def P_L_f(T_L, z, R_L):
    r"""
    P_L = 97000 \cdot e^{\frac{-g \cdot z}{R_L \cdot T_L }}
    """
    # (12) pg 27
    g = 9.81 * ureg.m / ureg.s ** 2
    return 97000 * ureg.pascal * math.exp(-g * z / (R_L * (T_L.to('kelvin'))))


def eta_A_f(t_m):
    r"""
    \eta_A = 15 \cdot 10^{-6} + 47 \cdot 10^{-9} \cdot t_m - 20 \cdot 10^{-12} \cdot t_m^2
    \quad\text{avec } t_m \text{ en °C}
    """
    # pg 100
    if hasattr(t_m, "to"):
        t_m = t_m.to(ureg.degC).magnitude
    value = 15e-6 + 47e-9 * t_m - 20e-12 * (t_m ** 2)
    return value * ureg.pascal * ureg.second


def lambda_A_f(t_m):
    r"""
    \lambda_A = 0.0223 + 0.000065 \cdot t_m
    \quad\text{avec } t_m \text{ en °C}
    """
    # B.9 pg 99
    if hasattr(t_m, "to"):
        t_m = t_m.to(ureg.degC).magnitude

    value = 0.0223 + 0.000065 * t_m
    return value * ureg.watt / (ureg.meter * ureg.kelvin)


def R_f(R_L, f_R, sigma_CO2):
    r"""
    R = R_L \cdot \left( 1 + f_R \cdot \sigma(CO_2) \right)
    """
    return R_L * (1 + f_R * sigma_CO2)


def c_p_f(t_m, sigma_CO2, f_c0, f_c1, f_c2, f_c3):
    r"""
    c_p = \frac{1011 + 0.05 \, t_m + 0.0003 \, t_m^2 + (f_{c0} + f_{c1} \, t_m + f_{c2} \, t_m^2)\,\sigma(CO_2)}{1 + f_{c3}\,\sigma(CO_2)}
    \quad\text{avec } t_m \text{ en °C}
    """
    # conversion en °C si pint.Quantity
    if hasattr(t_m, "to"):
        t_m = t_m.to(ureg.degC).magnitude

    numerator = 1011 + 0.05 * t_m + 0.0003 * (t_m ** 2) + (f_c0 + f_c1 * t_m + f_c2 * (t_m ** 2)) * sigma_CO2
    denominator = 1 + f_c3 * sigma_CO2

    return (numerator / denominator) * ureg.joule / (ureg.kg * ureg.kelvin)


def Nu_f(psi, psi_smooth, Re, Pr, D_h, L_tot):
    r"""
    Nu = \left(\frac{\psi}{\psi_{\text{smooth}}}\right)^{0.67}
         \cdot 0.0214 \cdot \left(Re^{0.8} - 100\right) \cdot Pr^{0.4}
         \cdot \left[1 + \left(\frac{D_h}{L_{tot}}\right)^{0.67}\right]
    """
    res = ((psi / psi_smooth) ** 0.67 * 0.0214 * (Re ** 0.8 - 100) * (Pr ** 0.4) * (1 + (D_h / L_tot) ** 0.67))
    # Ensure dimensionless result
    return res.magnitude


def colebrook_psi(Re, r, D_h, tol=1e-6, max_iter=20):
    r"""
    \frac{1}{\sqrt{\Psi}} = -2 \cdot \log_{10}\!\left(\frac{2{,}51}{Re \cdot \sqrt{\Psi}} + \frac{r}{3{,}71 \cdot D_h}\right)
    """
    # Ensure dimensionless inputs where required
    Re_val = Re
    rel_rough = (r / D_h)

    # Swamee–Jain explicit approximation as initial guess (robust start)
    # f ≈ 0.25 / [log10( (ε/3.7D) + 5.74/Re^0.9 )]^2
    denom = rel_rough / 3.7 + 5.74 / (Re_val ** 0.9)
    psi = 0.25 / (math.log10(denom) ** 2)

    # Fixed-point iterations on 1/sqrt(psi)
    for _ in range(max_iter):
        rhs = -2.0 * math.log10(2.51 / (Re_val * math.sqrt(psi)) + rel_rough / 3.71)
        new_psi = 1.0 / (rhs * rhs)
        if abs(new_psi - psi) <= tol * psi:
            psi = new_psi
            break
        psi = new_psi

    return psi


def solve_conduit(T_e, T_u, T_u0, L, dm, D_h, Res_therm):


    #################################Ss
    # VALEUR PROVISOIRS CALCULS
    ##################################

    # Température moyenne des fumées
    T_m = Q_(223, 'degC')
    T_m_prev = None

    ###################################

    while (T_m_prev is None or abs((T_m - T_m_prev)).magnitude > 1e-1):
        T_m_prev = T_m
        # Circonférence intérieure du conduit de fumée
        U = 3.1416 * D_h  # m

        # Section cheminée/appareil à combustion
        A = 3.1416 * (D_h / 2) ** 2  # m^2
        show_value(A)

        # Surface totale raccordement
        A_v = U * L_v  # m^2

        # Surface totale conduit de fumée
        A_t = U * L  # m^2
        show_value(A_t)

        # Coefficient de calcul de la teneur en vapeur d'eau des fumées
        # Annexe B3, valeurs pour gaz et tirage forcé
        f_w = 56
        # [56, 57]  # %

        f_m1 = [3.72, 3.75]
        f_m2 = [0.053, 0.054]

        # Puissance thermique nominale
        Q_F = 48 * ureg.kW
        # 48 - 195

        # Annexe B3
        f_x1 = 8.6
        f_x2 = 0.078
        f_x3 = 10.2

        # Puissance utile nominale
        Q_N = Q_F * 0.85

        # Coefficient de calcul de la teneur en CO2 des fumées
        sigma_CO2 = sigma_CO2_f(Q_N=Q_N, f_x1=f_x1, f_x2=f_x2, f_x3=f_x3)
        show_latex(sigma_CO2_f.__doc__)
        show_value(sigma_CO2)

        # rendement de l'appareil à combustion (%)
        eta = eta_w(Q_N=Q_N)
        show_latex(eta_w.__doc__)
        show_value(eta)

        # Teneur en vapeur d'eau des fumées (%)
        sigma_H2O = sigma_H2O_f(f_w=f_w, sigma_CO2=sigma_CO2)
        show_latex(sigma_H2O_f.__doc__)

        # Hauteur au dessus du niveau de la mer
        z = 10 * ureg.m

        # Température de l'air extérieur
        T_L = Q_(-15, 'degC')
        # Prendre +15 °C pour tirage minimal ou pressio positive maximale ou niveau de l'admisson des fumées dans le conduits

        # Constant de gaz de l'air
        R_L = 288 * ureg.joule / (ureg.kg * ureg.kelvin)

        # Pression de l'air extérieur
        p_L = P_L_f(T_L=T_L, z=z, R_L=R_L)
        show_latex(P_L_f.__doc__)
        show_value(p_L)



        # Viscosité dynamique des fumées
        eta_A = eta_A_f(t_m=T_m)
        show_latex(eta_A_f.__doc__)
        show_value(eta_A, 'eta_A')

        # Coefficient de conductivité thermique des fumées
        lambda_A = lambda_A_f(t_m=T_m)
        show_latex(lambda_A_f.__doc__)
        show_value(lambda_A)

        # coefficient de calcul de la constante des gaz des fumÈes, en 1/%
        # pg 96
        f_R = 0.0033
        # ou 0.0032

        R = R_f(f_R=f_R, R_L=R_L, sigma_CO2=sigma_CO2)
        show_value(R.to('joule/(kg*kelvin)'))

        # Densité moyenne des fumées
        rho_m = p_L / (R * T_m.to('kelvin'))
        show_value(rho_m.to('kg/m^3'))

        # Vitesse moyenne des fumées
        w_m = dm / (rho_m * A)
        show_value(w_m.to('m/s'))

        # Nombre de Reynolds
        # Cette formule est équivalente à Re = (rho_m * w_m * D_h) / eta_A
        Re = (dm * D_h) / (A * eta_A)
        # si < 2300 alors Re = 2300
        Re = max(Re, 2300)
        show_value(Re, 'Re')

        f_c0 = 23  # 23.5
        f_c1 = 0.015
        f_c2 = -0.000007
        f_c3 = 0.0142  # 0.0144

        # Capacité calorifique specifique des fumées (pg 99)
        c_p = c_p_f(t_m=T_m, sigma_CO2=sigma_CO2, f_c0=f_c0, f_c1=f_c1, f_c2=f_c2, f_c3=f_c3)
        show_latex(c_p_f.__doc__)
        show_value(c_p, 'c_p')

        # Nombre de Prandtl
        Pr = (c_p * eta_A) / lambda_A
        show_value(Pr, 'Pr')
        if Pr < 0.6 or Pr > 1.5:
            show_latex("**Attention**: Pr n'est pas entre 0.6 et 1.5, les calculs peuvent être incorrects.")

        # valeur moyenne de rugosité tab. B.4 pg 102
        r = 0.001 * ureg.m

        # Coefficient de frottement
        psi = colebrook_psi(Re=Re, r=r, D_h=D_h)
        psi_smooth = colebrook_psi(Re=Re, r=0 * ureg.m, D_h=D_h)

        show_value(psi, 'psi')
        show_value(psi_smooth, 'psi_smooth')

        # Nombre de Nusselt
        Nu = Nu_f(psi=psi, psi_smooth=psi_smooth, Re=Re, Pr=Pr, D_h=D_h, L_tot=L + L_v)
        show_value(Nu, 'Nu')

        # Coefficient intener de transfert de chaleur par convection
        alpha_i = (Nu * lambda_A) / D_h
        show_value(alpha_i, 'alpha_i')



        # Coefficient externe de transfert thermique exterieur
        alpha_a_ext = 25 * ureg.watt / (ureg.meter ** 2 * ureg.kelvin)
        alpha_a_int = 8 * ureg.watt / (ureg.meter ** 2 * ureg.kelvin)
        # on prend une interpolation
        alpha_e = 10 * ureg.watt / (ureg.meter ** 2 * ureg.kelvin)

        # Coefficient de transfert thermique température NON équilibrée
        k = 1 / (1 / alpha_i + 0.5 * (Res_therm + 1 / alpha_e))
        show_value(k, 'k')

        # Coefficient de transfert thermique température équilibrée
        k_b = 1 / (1 / alpha_i + Res_therm + 1 / alpha_e)
        show_value(k_b, 'k_b')

        # Coefficient de refroidissement température NON équilibrée
        K = (k * U * L) / (dm * c_p)
        show_value(K, 'K')
        # Coefficient de refroidissement température équilibrée
        K_b = (k_b * U * L) / (dm * c_p)

        T_u = Q_(15, 'degC')

        # Température des fumées à la sortie  équilibrée
        T_ob = T_u + (T_e - T_u) * math.exp(-K_b)
        show_value(T_ob, 'T_ob')

        # Température moyenne des fumées
        T_m = T_u + (T_e - T_u) * (1 - math.exp(-K)) / K
        show_value(T_m, 'T_m')

        # Température de la paroi intérieure à la sortie
        T_iob = T_ob + (T_uo - T_ob) * (k_b / alpha_i)
        show_value(T_iob, 'T_iob')

    return {
        'A': A,
        'A_v': A_v,
        'A_t': A_t,
        'sigma_CO2': sigma_CO2,
        'eta': eta,
        'sigma_H2O': sigma_H2O,
        'p_L': p_L,
        'T_m': T_m,
        'eta_A': eta_A,
        'lambda_A': lambda_A,
        'R': R,
        'rho_m': rho_m,
        'w_m': w_m.to('m/s'),
        'Re': Re,
        'Pr': Pr,
        'c_p': c_p,
        'psi': psi,
        'psi_smooth': psi_smooth,
        'Nu': Nu,
        'alpha_i': alpha_i,
        'Res_therm': Res_therm,
        'k': k,
        'k_b': k_b,
        'K': K,
        'K_b': K_b,
        'T_ob': T_ob,
        'T_iob': T_iob
        }


# Température fumée sortie foyer)
# Pg 20
# la valeur doit être connue, ou en alternative il faut prendre les resultats des formules de l'annexe B
T_e = Q_(300, 'degC')


# Diamètre cheminée/appareil à combustion
D_h = 200 * ureg.mm

# Longueur du conduit de raccordement
L_v = 0.2 * ureg.m

# Longueur du conduit de fumée
L = 7.4 * ureg.m

# Débit massique (fourni)
dm = 0.01 * ureg.kg / ureg.s  # kg/s

# Débit massique calculé
# Pg 99
# Je divise par 1000 pour réobtenir le résultat en kg/s
m_prime_calc = [0.019646085, 0.084213628]  # kg/s

# Température air ambiant à la sortie du conduit de fumée
# 5.7.1.3 pg 25, 0 °C (sans condensation)
T_uo = Q_(0, 'degC')
# ou -15 °C (avec condensation)
# T_uo = Q_(-15, 'degC')

# Température air ambiant dans la salle des chaudières
# 15 °C
T_ub = Q_(15, 'degC')

# Température air ambiant dans les zones chauffées
# 20 °C
T_uh = Q_(20, 'degC')

# Température air ambiant extérieure bâtiment
# égale à T_uo
T_ul = T_uo

# Température air ambiant non chauffées intérieur du bâtiment
# 0 °C
T_uu = Q_(0, 'degC')

# Température pour calcul d'exigeance de température de fumée

# raccord, on va dire 15 °C
T_uv = T_ub

# conduit de fumée, etant donné que la partie de cheminée à l'extérieur est < 1/4 de la longueur totale, on peut prendre 15 °C (Note 1 pg 26)
T_ut = Q_(15, 'degC')

# resistance thermique
Res_therm = 0 / (ureg.watt / (ureg.meter ** 2 * ureg.kelvin))


res = solve_conduit(T_e=T_e, T_u=T_ut, T_u0=T_uo, L=L, dm=dm, D_h=D_h, Res_therm=Res_therm)

print(json.dumps(res, indent=2, default=str))