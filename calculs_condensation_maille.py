import math

from calculs import alpha_i
from report import *
import json
from joblib import Memory



def m_dot(Q_F, f_m1, f_m2, sigma_CO2):
    r"""
    \dot{m} = \left(\frac{f_{m1}}{\sigma(\text{CO}_2)} + f_{m2}\right) \cdot Q_F
    """
    return (f_m1 / sigma_CO2 + f_m2) * Q_F



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


def p_L_f(T_L, z, R_L):
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


def R_f(R_L, sigma_CO2, avec_condensation=False):
    r"""
    R = R_L \cdot \left( 1 + f_R \cdot \sigma(CO_2) \right)
    """
    f_R = 0.00325 if not avec_condensation else 0.00025

    return R_L * (1 + f_R * sigma_CO2)


def c_p_f(t_m, sigma_CO2, f_c0, f_c1, f_c2, f_c3):
    r"""
    c_p = \frac{1011 + 0.05 \, t_m + 0.0003 \, t_m^2 + (f_{c0} + f_{c1} \, t_m + f_{c2} \, t_m^2)\,\sigma(CO_2)}{1 + f_{c3}\,\sigma(CO_2)}
    \quad\text{avec } t_m \text{ en °C}
    """
    # conversion en °C
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

def f_K_f(T_e):
    r"""
    f_K = 132.7 - 2.6 \cdot t_e + 0.0133 \cdot t_e^2
    \quad\text{avec } t_e \text{ en °C}
    """
    # pag 89. Normalement la formule est donnée uniquement pour T_e entre 40 et 80 °C
    if hasattr(T_e, "to"):
        T_e = T_e.to(ureg.degC).magnitude

    # La formuls est donnée uniquement pour T_e entre 40 et 80 °C. Comme la fonction est décroissante, on limite la valeur
    if T_e < 40 or T_e > 80:
        print("**Attention**: T_e en dehors de l'intervalle 40-80 °C, la formule pour f_K peut être incorrecte.")
    T_e = max(40, min(80, T_e))

    return 132.7 - 2.6 * T_e + 0.0133 * (T_e ** 2)

def p_Do_j_f(T_iob_j):
    r"""
    P_{D0,j} = \exp\!\left( 23.6448 - \frac{4077.9}{T_{iob,j} - 36.48} \right)
    \quad\text{avec } T_{iob,j} \text{ en K, } P_{D0,j} \text{ en Pa}
    """
    # pg 88
    T_iob_j = T_iob_j.to(ureg.kelvin).magnitude
    val = math.exp(23.6448 - 4077.9 / (T_iob_j - 36.48))
    return val * ureg.pascal

def delta_m_D_j_f(m_dot, R, p_D, p_L, p_Do_prev, p_Do_j, f_K):
    r"""
    \Delta \dot{m}_{D,j} =
    \dot{m} \cdot \frac{R}{R_D} \cdot \left( 1 - \frac{p_D}{p_L} \right)
    \cdot \left( \frac{p_{Do,j-1}}{p_L - p_{Do,j-1}} - \frac{p_{Do,j}}{p_L - p_{Do,j}} \right)
    \cdot \frac{f_K}{100}
    """
    """
    Paramètres :
    ------------
    m_dot : débit massique des fumées avant condensation (kg/s)
    R : constante des gaz des fumées (J/(kg·K))
    p_D : pression partielle de la vapeur d'eau avant condensation (Pa)
    p_L : pression de l'air extérieur (Pa)
    p_Do_prev : pression partielle vapeur eau liée à T_iob du segment précédent (Pa)
    p_Do_j : pression partielle vapeur eau liée à T_iob du segment courant (Pa)
    f_K : facteur de condensation (%)
    """
    R_D = 496 * ureg.joule / (ureg.kg * ureg.kelvin)
    term1 = m_dot * (R / R_D)
    term2 = (1 - p_D / p_L)
    term3 = (p_Do_prev / (p_L - p_Do_prev)) - (p_Do_j / (p_L - p_Do_j))
    val = term1 * term2 * term3 * (f_K / 100.0)
    val = max(0 * ureg.kg / ureg.s, val)  # ensure non-negative
    return val.to('kg/s')

def T_ob_j_f(T_ob_prev, K_b_j, q_K_j, m_dot_o_j, m_dot_o_prev, c_po_j, c_po_j_prev, T_u_j):
    r"""
    T_{ob,j} =
    \frac{\left(\frac{\dot m_{o,j-1}\,c_{po,j-1}}{\dot m_{o,j}\,c_{po,j}} - \frac{K_{b,j}}{2}\right) T_{ob,j-1}
          + \frac{q_{K,j}}{\dot m_{o,j}\,c_{po,j}} + K_{b,j}\,T_{u,j}}
         {1 + \frac{K_{b,j}}{2}}
    \quad\text{(T en K, }K_{b,j}\text{ sans dimension)}
    """
    # EN 13384-1, eq. (165)
    T_ob_prev = T_ob_prev.to(ureg.kelvin)
    T_u_j = T_u_j.to(ureg.kelvin)
    K_b_j = K_b_j.to('dimensionless')

    ratio = (m_dot_o_prev * c_po_j_prev) / (m_dot_o_j * c_po_j)  # dimensionless
    denominator = (1.0 + K_b_j/2.0)

    terme_refroidissement = (((m_dot_o_prev * c_po_j_prev) / (m_dot_o_j * c_po_j) - (K_b_j / 2)) * T_ob_prev + K_b_j * T_u_j) / denominator

    contr_condens= (q_K_j / (m_dot_o_j * c_po_j) / denominator).to('delta_degC')

    val = terme_refroidissement + contr_condens

    return val.to(ureg.kelvin), contr_condens

def Re_f(w_m, D_h, rho_m, eta_A):
    r"""
    Re = \frac{w_m \cdot D_h}{\eta_A}
    \quad\text{avec } w_m \text{ en m/s, } D_h \text{ en m, } \eta_A \text{ en Pa·s}
    """
    if w_m < 0.5 * ureg.m / ureg.s:
        # pg. 31 "Pour une vitesse moyenne des fumées w_m < 0.5 m/s, prendre un nombre de Nusselt approprrié à w_m = 0.5 m/s."
        w_m = 0.5 * ureg.m / ureg.s
    Re = (w_m * D_h * rho_m) / eta_A
    Re = Re.to('dimensionless')
    return max(Re.magnitude, 2300)  # minimum 2300

# def m_dot_theor


def solve_conduit_helper(T_m, T_in, T_u, T_uo, L, m_dot, D_h, Res_therm, alpha_a, T_L=None, T_e=None, data_j_prev=None, L_tot=None) -> dict:
    """
    Calculs communs aux segments avec et sans condensation.
    """
    if T_L is None:
        T_L = T_uo

    # Q est la puissance utile
    # Q_N est la puissance utile nominale
    # Q_F est le debit calorifique

    # Nous devons valider les exigencer pour la puissance utile nomintale de l'appareil à combustion déclarée par le fabricant
    # ous sinon la puissance utile nominale et la plus petite valeur de la plage de puissance utile declarée par le fabricant
    # Pour nous le discours est un peu compliqué, car nous avons un bruleur + foyer
    # d'après les essais de Frank, le rendement du foyer est de 86.5 %
    # la puissance du bruleur est de 48/79 - 195 kW mais je n'ai pas de nominale
    # je vais donc prendre le nominal le plus bas (condition la plus défavorable) -> 48 kW

    Q = 140 * ureg.kW
    Q_F = Q / 0.86
    # et donc on va prendre le nominal le plus bas
    Q_N = Q


    # 48 - 195
    # Hauteur au dessus du niveau de la mer
    z = 10 * ureg.m

    # Coefficient de calcul de la teneur en vapeur d'eau des fumées
    # Annexe B3, valeurs pour gaz et tirage forcé
    f_x1 = 8.6
    f_x2 = 0.078
    f_x3 = 10.2
    # B.1
    # f_w = 56.5 # pour le fioul
    f_w = 111 # gaz nat


    # Coefficient de calcul de la teneur en CO2 des fumées
    sigma_CO2 = sigma_CO2_f(Q_N=Q_N, f_x1=f_x1, f_x2=f_x2, f_x3=f_x3)
    # mesuration OMIA
    # sigma_CO2 = 14

    # Teneur en vapeur d'eau des fumées (%)
    sigma_H2O = sigma_H2O_f(f_w=f_w, sigma_CO2=sigma_CO2)
    # show_latex(sigma_H2O_f.__doc__)

    # Constant de gaz de l'air
    R_L = 288 * ureg.joule / (ureg.kg * ureg.kelvin)
    R = R_f(R_L=R_L, sigma_CO2=sigma_CO2, avec_condensation=True)

    # Pression de l'air extérieur
    p_L = p_L_f(T_L=T_L, z=z, R_L=R_L)
    # show_latex(P_L_f.__doc__)
    # show_value(p_L)

    # Pression partielle de la vapeur d'eau dans les fumées
    p_D = (sigma_H2O / 100) * p_L

    # Température de rosée
    t_p = 4077.9 / (23.6448 - math.log(p_D.magnitude)) - 236.67
    t_p = Q_(t_p, 'degC')

    # Circonférence intérieure du conduit de fumée
    U = 3.1416 * D_h  # m

    # Section cheminée/appareil à combustion
    A = 3.1416 * (D_h / 2) ** 2  # m^2
    # show_value(A)

    # Surface totale conduit de fumée
    A_t = U * L  # m^2
    # show_value(A_t)

    # rendement de l'appareil à combustion (%)
    eta = eta_w(Q_N=Q_N)
    # show_latex(eta_w.__doc__)
    # show_value(eta)

    # Viscosité dynamique des fumées
    eta_A = eta_A_f(t_m=T_m)
    # show_latex(eta_A_f.__doc__)
    # show_value(eta_A, 'eta_A')

    # Coefficient de conductivité thermique des fumées
    lambda_A = lambda_A_f(t_m=T_m)

    # Densité moyenne des fumées
    rho_m = p_L / (R * T_m.to('kelvin'))
    # show_value(rho_m.to('kg/m^3'))

    # Vitesse moyenne des fumées, il faut utiliser le nouveau debit massique apres condensation
    m_dot_segment = data_j_prev['m_dot_o_j'] if data_j_prev is not None else m_dot
    w_m = m_dot_segment / (rho_m * A)
    # show_value(w_m.to('m/s'))

    # Nombre de Reynolds
    Re = Re_f(w_m=w_m, D_h=D_h, rho_m=rho_m, eta_A=eta_A)

    # show_value(Re, 'Re')
    # B.1
    f_c0 = 23  # 23.5
    f_c1 = 0.015
    f_c2 = -0.000007
    f_c3 = 0.0142  # 0.0144

    # Capacité calorifique specifique des fumées (pg 99)
    c_p = c_p_f(t_m=T_m, sigma_CO2=sigma_CO2, f_c0=f_c0, f_c1=f_c1, f_c2=f_c2, f_c3=f_c3)

    # Nombre de Prandtl
    Pr = (c_p * eta_A) / lambda_A

    if Pr < 0.6 or Pr > 1.5:
        print("**Attention**: Pr n'est pas entre 0.6 et 1.5, les calculs peuvent être incorrects.")

    # valeur moyenne de rugosité tab. B.4 pg 102
    r = 0.001 * ureg.m

    # Coefficient de frottement
    psi = colebrook_psi(Re=Re, r=r, D_h=D_h)
    psi_smooth = colebrook_psi(Re=Re, r=0 * ureg.m, D_h=D_h)


    # Nombre de Nusselt
    Nu = Nu_f(psi=psi, psi_smooth=psi_smooth, Re=Re, Pr=Pr, D_h=D_h, L_tot=L_tot)

    # Coefficient intener de transfert de chaleur par convection
    alpha_i = (Nu * lambda_A) / D_h
    # convert to SI
    alpha_i = alpha_i.to_base_units()  # show_value(alpha_i, 'alpha_i')

    res = {
        'A': A,
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
        'alpha_i': alpha_i.to_base_units(),
        'alpha_a': alpha_a.to_base_units(),
        't_p': t_p,
        'U': U,
        'p_D': p_D,
        }
    return res


def solve_conduit_sans_condensation(T_m, T_in, T_u, T_uo, L, m_dot, D_h, Res_therm, alpha_a, T_e=None, data_j_prev=None, L_tot=None, T_L=None) -> dict:
    res = solve_conduit_helper(T_m=T_m, T_in=T_in, T_u=T_u, T_uo=T_uo, L=L, m_dot=m_dot, D_h=D_h, Res_therm=Res_therm, alpha_a=alpha_a, T_e=T_e, data_j_prev=data_j_prev, L_tot=L_tot)

    # Extract values from helper result
    alpha_i = res['alpha_i']
    t_p = res['t_p']
    c_p = res['c_p']
    U = res['U']

    # Coefficient de transfert thermique température NON équilibrée
    k = 1 / (1 / alpha_i + 0.5 * (Res_therm + 1 / alpha_a))
    # show_value(k, 'k')

    # Coefficient de transfert thermique température équilibrée
    k_b = 1 / (1 / alpha_i + Res_therm + 1 / alpha_a)
    # show_value(k_b, 'k_b')

    # Coefficient de refroidissement température NON équilibrée
    K = (k * U * L) / (m_dot * c_p)
    # check K is dimensionless
    K = K.to('dimensionless')
    # show_value(K, 'K')

    # Coefficient de refroidissement température équilibrée
    K_b = (k_b * U * L) / (m_dot * c_p)
    # check K_b is dimensionless
    K_b = K_b.to('dimensionless')

    # Température des fumées à la sortie  équilibrée
    T_ob = T_u.to('kelvin') + (T_in.to('kelvin') - T_u.to('kelvin')) * math.exp(-K_b)
    # show_value(T_ob, 'T_ob')

    # Température moyenne des fumées
    T_m_prev = T_m
    # T_m = T_u + (T_in - T_u) * (1 - math.exp(-K)) / K
    T_m = T_u.to('kelvin') + (T_in.to('kelvin') - T_u.to('kelvin')) * (1 - math.exp(-K)) / K
    # show_value(T_m, 'T_m')

    # Température de la paroi intérieure à la sortie
    T_iob = T_ob + (T_uo - T_ob) * (k_b / alpha_i)  # show_value(T_iob, 'T_iob')

    p_Do_j = p_Do_j_f(T_iob_j=T_iob)

    # Calcul point de rosée
    condensation = T_iob.to('degC') < t_p.to('degC')

    res.update({
        'type_calcul': 'segment_sans_condensation',
        'T_m': T_m,
        'alpha_i': alpha_i,
        'alpha_a': alpha_a,
        'Res_therm': Res_therm,
        'k': k,
        'k_b': k_b,
        'K': K,
        'K_b': K_b,
        'T_ob': T_ob,
        'T_iob': T_iob,
        't_p': t_p,
        'condensation': condensation,
        'p_Do_j': p_Do_j,
        'm_dot_o_j': m_dot,
        'c_po_j': c_p,
        'k_btot_j': k_b,
        'L': L,
        'T_u': T_u,
        'contr_condens': 0 * ureg.delta_degC,
        })

    return res


def solve_conduit_avec_condensation(T_iob, T_ob, T_in, T_u, T_uo, L, m_dot, D_h, Res_therm, alpha_a, T_e=None, data_j_prev=None, L_tot=None):
    T_m = (T_in.to('kelvin') + T_ob.to('kelvin')) / 2.0

    res = solve_conduit_helper(T_m=T_m, T_in=T_in, T_u=T_u, T_uo=T_uo, L=L, m_dot=m_dot, D_h=D_h, Res_therm=Res_therm, alpha_a=alpha_a, T_e=T_e, data_j_prev=data_j_prev, L_tot=L_tot)
    # Extract values from helper result
    # Extract values from helper result
    alpha_i = res['alpha_i']
    t_p = res['t_p']
    c_p = res['c_p']
    U = res['U']
    R = res['R']
    p_L = res['p_L']
    p_D = res['p_D']

    p_Do_j = p_Do_j_f(T_iob_j=T_iob)

    f_K = f_K_f(T_e=T_e)
    # show_value(f_K, 'f_K')

    if T_iob < data_j_prev['T_iob']:
        delta_m_D_j = delta_m_D_j_f(m_dot=m_dot, R=R, p_D=p_D, p_L=p_L, p_Do_prev=data_j_prev['p_Do_j'], p_Do_j=p_Do_j, f_K=f_K)
    else:
        delta_m_D_j = Q_(0, 'kg/s')

    # ENthalpie de l'eau évaporée pg 88
    r_D = 2_400_000 * ureg.joule / ureg.kg

    # Chaleur de condensation
    q_K_j = delta_m_D_j * r_D

    # l_c_j est lavProportion de surface de condensation dans le segment j
    # TODO: fix below: something is wrong with the logic
    is_first_condensing_segment = (data_j_prev['type_calcul'] != 'segment_avec_condensation')
    if not is_first_condensing_segment:
        # Ce n'est pas le premier segment avec condensation
        l_c_j = 1
    else:
        # Premier segment avec condensation
        T_iob_j_prev = data_j_prev['T_iob']
        # La température de rosée à l'entrée du premier segment (sauf que normalement elle ne dépende pas de j, donc je prends T_p)
        T_pe_1 = t_p
        # l_c_j doit être entre 0 et 1. Un crop assure que la condition soit respectée pendant les itérations
        l_c_j = max(0.1, min(1, (T_pe_1.to('degC').magnitude - T_iob.to('degC').magnitude) / (T_iob_j_prev.to('degC').magnitude - T_iob.to('degC').magnitude)))
        # L'introduction de l_c_j introduit des instabilités dans la convergence, on prend une valeur
        # de 0.5, ce qui implique un erreur maxi de 15 cm sur la longueur de la condensation, ce qui est négligeéble sur 7 m de cheminée
        l_c_j = 0.5

    # Coefficient de transfert thermique par condensation
    alpha_ioK_j = q_K_j / (l_c_j * (T_ob - T_iob) * U * L)

    # Coefficient tolla de transfert thermique, qui est donc l'effet combiné de la convection et de la condensation (pg. 86)
    alpha_iotot_j = alpha_i + alpha_ioK_j

    k_obtot_j = 1 / (1 / alpha_iotot_j + Res_therm + 1 / alpha_a)

    # Debit massique des fumées à la sortie du segment

    m_dot_o_j = data_j_prev['m_dot_o_j'] - delta_m_D_j

    ## AJOUTER ICI !!!!

    if is_first_condensing_segment:
        k_b_prev = data_j_prev['k_b']
        k_btot_j = (1 - l_c_j) * k_b_prev + l_c_j * k_obtot_j
    else:
        k_btot_prev = data_j_prev['k_btot_j']
        k_btot_j = (k_btot_prev + k_obtot_j) / 2

    # Calcul déjà fait, qui va converger
    c_po_j = c_p
    K_b = (U * k_btot_j * L) / (m_dot_o_j * c_po_j)
    K_b = K_b.to('dimensionless')

    T_ob_j, contr_condens = T_ob_j_f(T_ob_prev=data_j_prev['T_ob'], K_b_j=K_b, q_K_j=q_K_j, m_dot_o_j=m_dot_o_j, m_dot_o_prev=data_j_prev['m_dot_o_j'], c_po_j=c_po_j, c_po_j_prev=data_j_prev['c_po_j'], T_u_j=T_u)
    T_ob = T_ob_j

    T_m_prev = T_m
    T_m = (T_ob.to('kelvin') + T_in.to('kelvin')) / 2

    T_iob = T_ob + (T_uo - T_ob) * (k_obtot_j / alpha_iotot_j)

    # Nous sommes dans un segment avec condensation, donc T_iob < t_p
    T_iob = min(T_iob, t_p - Q_(0.1, 'delta_degC')).to('kelvin')

    # Calcul point de rosée
    condensation = T_iob.to('degC') < t_p.to('degC')

    q_A = data_j_prev['m_dot_o_j'] * data_j_prev['c_po_j'] * data_j_prev['T_ob'] - m_dot_o_j * c_po_j * T_ob_j + q_K_j
    q_C = alpha_iotot_j * U * L * (data_j_prev['T_ob'] - data_j_prev['T_iob'] + T_ob_j - T_iob) / 2

    res.update({
        'type_calcul': 'segment_avec_condensation',
        'p_L': p_L,
        'T_m': T_m,
        'm_dot': m_dot,
        'R': R,
        'c_p': c_p,
        'alpha_a': alpha_a,
        'Res_therm': Res_therm,
        'K_b': K_b,
        'T_ob': T_ob,
        'T_iob': T_iob,
        't_p': t_p,
        'condensation': condensation,
        'p_Do_j': p_Do_j,
        'L': L,
        'f_K': f_K,
        'delta_m_D_j': delta_m_D_j,
        'q_K_j': q_K_j,
        'l_c_j': l_c_j,
        'alpha_ioK_j': alpha_ioK_j,
        'alpha_iotot_j': alpha_iotot_j,
        'k_obtot_j': k_obtot_j,
        'm_dot_o_j': m_dot_o_j,
        'k_btot_j': k_btot_j,
        'c_po_j': c_po_j,
        'T_ob_j': T_ob_j,
        'T_iob_j': T_iob,
        'q_A': q_A,
        'q_C': q_C,
        'T_u': T_u,
        'contr_condens': contr_condens,
    })

    return res



def solve_conduit_minimize_method(T_in, T_u, T_uo, L, m_dot, D_h, Res_therm, alpha_a, avec_condensation=False, T_e=None, data_j_prev=None, L_tot=None, T_L=None) -> dict:

    #################################Ss
    # VALEUR PROVISOIRS CALCULS
    ##################################

    # Température moyenne des fumées
    T_m = Q_(223, 'degC')
    T_m_prev = Q_(99999999, 'degC')

    # Si il y avait de la condensation dans le segment précedent, on passe directement au calcul avec condensation
    if not data_j_prev or data_j_prev['type_calcul'] != 'segment_avec_condensation':
        # Commence le calcul sans condensation
        while (abs(T_m.to('kelvin').magnitude - T_m_prev.to('kelvin').magnitude) > 1):

            res = solve_conduit_sans_condensation(T_m=T_m, T_in=T_in, T_u=T_u, T_uo=T_uo, L=L, m_dot=m_dot, D_h=D_h, Res_therm=Res_therm, alpha_a=alpha_a, T_e=T_e, data_j_prev=data_j_prev, L_tot=L_tot)
            T_m_prev = T_m
            T_m = res['T_m']

        print("Calculs sans condensation terminés.")
        do_calcule_condensation = res['condensation']
        if do_calcule_condensation:
            print("Le segment précédent est sans condensation, mais le critère de condensation est maintenant respecté, on passe au calcul avec condensation.")
    else:
        do_calcule_condensation = True
        print("Le segment précédent est avec condensation, on passe directement au calcul avec condensation.")
        # On prend les données du segment précédent pour initialiser le calcul avec condensation
        res = data_j_prev

    if do_calcule_condensation:
        # check data_j_prev is not null
        assert data_j_prev is not None, "data_j_prev ne peut pas être None si on calcule avec condensation."
        print("Calculs avec condensation en cours...")
        # Recommence le calcul avec condensation
        T_ob_min = T_uo.to('kelvin').magnitude
        T_ob_max = data_j_prev['T_ob'] + Q_(0.1, 'delta_degC')
        T_ob_max = T_ob_max.to('kelvin').magnitude

        ratio_T_iob_min = 0.001
        ratio_T_iob_max = 0.999

        prev_T_ob = data_j_prev['T_ob'].to('kelvin').magnitude
        prev_ratio_T_iob = (prev_T_ob - data_j_prev['T_iob'].to('kelvin').magnitude) / (prev_T_ob - T_uo.to('kelvin').magnitude)

        class EarlyStop(Exception):
            pass
        eps = 1e-4

        def f_to_solve(vals):
            T_ob, ratio_T_iob = vals
            T_iob = T_ob * ratio_T_iob + T_uo.to('kelvin').magnitude * (1 - ratio_T_iob)
            T_iob_K = Q_(T_iob, 'kelvin')
            T_ob_K = Q_(T_ob, 'kelvin')
            res = solve_conduit_avec_condensation(T_iob=T_iob_K, T_ob=T_ob_K, T_in=T_in, T_u=T_u, T_uo=T_uo, L=L, m_dot=m_dot, D_h=D_h, Res_therm=Res_therm, alpha_a=alpha_a, T_e=T_e, data_j_prev=data_j_prev, L_tot=L_tot)
            T_ob_new = res['T_ob'].to('kelvin').magnitude
            ratio_T_iob_new = (res['T_iob'] - T_uo).to('kelvin').magnitude / (res['T_ob'] - T_uo).to('kelvin').magnitude
            rel_err = math.fabs((T_ob_new - T_ob) / (T_ob - T_uo.to('kelvin').magnitude)) + math.fabs(ratio_T_iob_new - ratio_T_iob)
            if rel_err < eps:
                raise EarlyStop((rel_err, (T_ob, ratio_T_iob)))
            return rel_err

        from scipy.optimize import root

        # sol = root(lambda x: f_to_solve(x[0], x[1]), x0=[(T_ob_min + T_ob_max) / 2, (ratio_T_iob_min + ratio_T_iob_max) / 2], method='hybr')
        from scipy.optimize import minimize
        bounds = [(T_ob_min, T_ob_max), (ratio_T_iob_min, ratio_T_iob_max)]

        try:
            sol = minimize(f_to_solve, x0=[prev_T_ob, prev_ratio_T_iob], bounds=bounds, method='Powell', options={'ftol': 1e-6, 'disp': False, 'maxiter': 400})
            # If we get here, early-stop was not triggered
            f_min, x_min = sol.fun, sol.x
            success = (f_min < 2e-2)
        except EarlyStop as e:
            f_min, x_min = e.args[0]
            success = True  # because we stopped exactly when f<eps

        if not success:
            raise Exception("La méthode de résolution n'a pas convergé.")

        if success:
            T_ob_sol, ratio_T_iob_sol = x_min
            T_iob_sol = T_ob_sol * ratio_T_iob_sol + T_uo.to('kelvin').magnitude * (1 - ratio_T_iob_sol)
            T_iob_K = Q_(T_iob_sol, 'kelvin')
            T_ob_K = Q_(T_ob_sol, 'kelvin')
            res = solve_conduit_avec_condensation(T_iob=T_iob_K, T_ob=T_ob_K, T_in=T_in, T_u=T_u, T_uo=T_uo, L=L, m_dot=m_dot, D_h=D_h, Res_therm=Res_therm, alpha_a=alpha_a, T_e=T_e, data_j_prev=data_j_prev)
            q_A = res['q_A']
            q_C = res['q_C']

            delta_iter_cond = math.fabs(q_A.to('watt').magnitude - q_C.to('watt').magnitude) / q_A.to('watt').magnitude
            print(f'{delta_iter_cond=}')
            if delta_iter_cond > 0.2:
                print("La méthode de résolution n'a pas convergé.")
                # raise Exception("Critère (169) pg 91 non respecté.")
        else:
            raise Exception("La méthode de résolution n'a pas convergé.")


    else:
        print("Le conduit n'est pas en condensation.")
        return res
    print("Calculs avec condensation terminés.")
    return res


def solve_conduit_root_method(*args, **kwargs):


    #################################Ss
    # VALEUR PROVISOIRS CALCULS
    ##################################

    # Température moyenne des fumées
    T_m = Q_(223, 'degC')
    T_m_prev = Q_(99999999, 'degC')


    # Commence le calcul sans condensation

    while (abs(T_m.to('kelvin').magnitude - T_m_prev.to('kelvin').magnitude) > 1):

        res = solve_conduit_sans_condensation(T_m=T_m, *args, **kwargs)
        T_m_prev = T_m
        T_m = res['T_m']

    print("Calculs sans condensation terminés.")
    if res['condensation']:
        print("Le conduit est en condensation, recalcul avec condensation...")
        # Recommence le calcul avec condensation
        T_uo = kwargs.get('T_uo')
        T_ob_min = T_uo.to('kelvin').magnitude
        T_ob_max = res['T_ob'] + Q_(0.1, 'delta_degC')
        T_ob_max = T_ob_max.to('kelvin').magnitude

        delta_T_iob_min = 1
        delta_T_iob_max = T_ob_max - T_uo.to('kelvin').magnitude

        def f_to_solve(T_ob, delta_T_iob):
            T_iob_K = Q_(T_ob - delta_T_iob, 'kelvin')
            T_ob_K = Q_(T_ob, 'kelvin')
            res = solve_conduit_avec_condensation(T_ob=T_ob_K, T_iob=T_iob_K, *args, **kwargs)
            T_ob_new = res['T_ob'].to('kelvin').magnitude
            delta_T_iob_new = T_ob_new - res['T_iob'].to('kelvin').magnitude
            r_ob = (T_ob_new - T_ob)
            r_delta_T_iob = (delta_T_iob_new - delta_T_iob)
            return r_ob, r_delta_T_iob

        from scipy.optimize import root

        sol = root(lambda x: f_to_solve(x[0], x[1]), [(T_ob_min + T_ob_max) / 2.0, (delta_T_iob_min + delta_T_iob_max) / 2.0], method='hybr')
        if sol.success:
            T_ob_sol, delta_T_iob_sol = sol.x
            T_iob_sol = T_ob_sol - delta_T_iob_sol
            T_iob_K = Q_(T_iob_sol, 'kelvin')
            T_ob_K = Q_(T_ob_sol, 'kelvin')
            res = solve_conduit_avec_condensation(T_ob=T_ob_K, T_iob=T_iob_K, *args, **kwargs)

            q_A = res['q_A']
            q_C = res['q_C']

            delta_iter_cond = math.fabs(q_A.to('watt').magnitude - q_C.to('watt').magnitude) / q_A.to('watt').magnitude
            print(f'{delta_iter_cond=}')
            if delta_iter_cond > 0.2:
                print("La méthode de résolution n'a pas convergé.")
                # raise Exception("Critère (169) pg 91 non respecté.")
        else:
            raise Exception("La méthode de résolution n'a pas convergé.")


    else:
        print("Le conduit n'est pas en condensation.")
        return res
    print("Calculs avec condensation terminés.")
    return res


def main():

    # Coefficient externe de transfert thermique exterieur
    alpha_a_ext = 25 * ureg.watt / (ureg.meter ** 2 * ureg.kelvin)
    alpha_a_int = 8 * ureg.watt / (ureg.meter ** 2 * ureg.kelvin)


    # Diamètre cheminée/appareil à combustion
    D_h = 200 * ureg.mm
    # Débit massique (fourni)
    # m_dot = 0.068 /3 * ureg.kg / ureg.s  # kg/s
    m_dot = 0.068 * ureg.kg / ureg.s  # kg/s
    # Pg 20
    # la valeur doit être connue, ou en alternative il faut prendre les resultats des formules de l'annexe B
    T_e = Q_(207, 'degC')
    # T_e = Q_(310, 'degC')
    T_L = Q_(-15, 'degC')
    all_data = []
    # Calculs circuit de raccordement

    # Nombre de segments dans le circuit de raccordement
    NsegV = 1

    print("Circuit de raccordement")
    L_v = 0.2 * ureg.m
    T_in = T_e
    T_uo = T_u = Q_(15, 'degC') # raccordement dans le batiment
    Res_therm = 0 / (ureg.watt / (ureg.meter ** 2 * ureg.kelvin))
    res_racc = solve_conduit_root_method(T_in=T_in, T_u=T_u, T_uo=T_uo, L=L_v, m_dot=m_dot, D_h=D_h, T_L=T_L, Res_therm=Res_therm, alpha_a=alpha_a_int, L_tot=L_v, data_j_prev=None, T_e=T_e)
    print(json.dumps(res_racc, indent=2, default=str))

    all_data.append({**res_racc, 'segment': 0, 'L_cum': L_v })
    # ==============================
    # Calculs conduit de fumée
    # ==============================

    # Nombre de segments dans le conduit de fumée
    NSegK = math.inf

    L_tot = 6.3 * ureg.m + 1.2 * ureg.m
    L_tot_batiment = 6.3 * ureg.m

    L_cum = L_v

    Nseg = 25
    # Dans le batiment
    T_in = res_racc['T_ob']
    res = None

    for i in range(Nseg):
        L = L_tot / Nseg
        L_cum += L
        print(f"Segment chéminée {i+1}/{Nseg}, longueur cumulée {L_cum.to('m')}")


        if L_cum <= L_tot_batiment:
            print("Le segment est entièrement à l'intérieur du bâtiment")
            T_uo = T_u = Q_(15, 'degC')  # conduit dans le batiment
            Res_therm = 0 / (ureg.watt / (ureg.meter ** 2 * ureg.kelvin))
            alpha_a = alpha_a_int
        else:
            print("Le segment est au moins en partie à l'extérieur du bâtiment")
            T_uo = T_u = Q_(-15, 'degC')  # conduit dans le batiment
            # 45 mm de laine de roche → Res_therm ~ 1 (m²K)/W
            Res_therm = 0 / (ureg.watt / (ureg.meter ** 2 * ureg.kelvin))
            alpha_a = alpha_a_ext


        res = solve_conduit_root_method(T_in=T_in, T_u=T_u, T_uo=T_uo, L=L, m_dot=m_dot, D_h=D_h, Res_therm=Res_therm, alpha_a=alpha_a, data_j_prev=res if i > 0 else None, L_tot=L_tot, T_e=T_e)
        all_data.append({**res, 'segment': i+1, 'L_cum': L_cum})
        T_in = res['T_ob']
        print(json.dumps(res, indent=2, default=str))
        if res['condensation']:
            print(f"**Attention**: Condensation détectée dans le segment {i+1}/{Nseg} du conduit")
            NSegK = min(i+1, NSegK)

    import matplotlib.pyplot as plt

    # Préparation des données
    x_array = [data['L_cum'].to('m').magnitude for data in all_data]
    T_ob_array = [data['T_ob'].to('degC').magnitude for data in all_data]
    T_iob_array = [data['T_iob'].to('degC').magnitude for data in all_data]
    contr_condens_array = [data['contr_condens'].to('delta_degC').magnitude for data in all_data]
    w_m_array = [data['w_m'].to('m/s').magnitude for data in all_data]
    Re_array = [data['Re'] for data in all_data]
    m_dot_array = [data['m_dot_o_j'].to("kg/s").magnitude for data in all_data]

    # Création d'une figure avec 4 sous-graphiques alignés verticalement
    fig, axes = plt.subplots(4, 1, figsize=(10, 18), sharex=True)

    # --- Plot 1: Températures ---
    axes[0].plot(x_array, T_ob_array, label='T_ob (°C)')
    axes[0].plot(x_array, T_iob_array, label='T_iob (°C)')
    axes[0].plot(x_array, contr_condens_array, label='contr condens', color='black', linestyle='--')
    axes[0].axhline(y=res['t_p'].to('degC').magnitude, color='r', linestyle='--', label='T condensation (°C)')
    axes[0].axhline(y=0, color='purple', linestyle='--', label='0 °C')
    axes[0].axvline(x=L_tot_batiment.to('m').magnitude, color='g', linestyle='--', label='Extérieur')
    axes[0].set_ylabel('Température (°C)')
    axes[0].set_title('Profil de température')
    axes[0].legend()
    axes[0].grid()
    axes[0].set_yticks(range(-10, 301, 10))

    # --- Plot 2: Vitesse fumées ---
    axes[1].plot(x_array, w_m_array, label='w_m (m/s)', color='orange')
    axes[1].set_ylabel('Vitesse (m/s)')
    axes[1].set_title('Vitesse des fumées')
    axes[1].legend()
    axes[1].grid()

    # --- Plot 3: Reynolds ---
    axes[2].plot(x_array, Re_array, label='Re', color='brown')
    axes[2].set_ylabel('Re')
    axes[2].set_title('Nombre de Reynolds')
    axes[2].legend()
    axes[2].grid()

    # --- Plot 4: Débit massique ---
    axes[3].plot(x_array, m_dot_array, label='m_dot (kg/s)', color='blue')
    axes[3].set_xlabel('Longueur cumulée (m)')
    axes[3].set_ylabel('Débit massique (kg/s)')
    axes[3].set_title('Débit massique des fumées')
    axes[3].legend()
    axes[3].grid()

    # Ajuster l’espacement
    plt.tight_layout()
    plt.show()


    # convert the dict to a panda dataframe
    import pandas as pd
    df = pd.DataFrame(all_data)
    # adjust df e façon à utiliser des unités SI using to._base_units()
    df_si = df.applymap(lambda x: x.to_base_units().magnitude if hasattr(x, 'to_base_units') else x)

    df_si.to_csv('calculs_condensation_maille.csv', index=False, sep=';')


    print("Done!")

if __name__ == "__main__":
    main()