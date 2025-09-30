import math

from report import *
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

    contr_condensation = (q_K_j / (m_dot_o_j * c_po_j) / denominator).to('kelvin')

    val = terme_refroidissement # + contr_condensation

    return val.to(ureg.kelvin)




def solve_conduit(T_in, T_u, T_uo, L, m_dot, D_h, Res_therm, alpha_a, avec_condensation=False, T_e=None, data_j_prev=None):
    """
    Résout les équations du conduit de fumée

    Paramètres
    ----------
    T_in : pint.Quantity
        Température fumée en entrée (K ou °C)
    T_u : pint.Quantity
        Température air autour du conduit (K ou °C)
    T_u0 : pint.Quantity
        Température air à la sortie du conduit (K ou °C)
    L : pint.Quantity
        Longueur du conduit de fumée (m)
    m_dot: pint.Quantity
        Débit massique (kg/s)
    D_h : pint.Quantity
        Diamètre cheminée/appareil à combustion (m)
    Res_therm : pint.Quantity
        Resistance thermique de la paroi (m²K/W)
    alpha_a : pint.Quantity
        Coefficient externe de transfert thermique (W/m²K) (8 pour intérieur, 25 pour extérieur)

    calculate_condensation : boolean
        Si True, calcule la présence de condensation

    T_e : pint.Quantity
        Température fumée sortie foyer (K ou °C). Si None, T_e = T_in

    data_j_prev : dict
        Résultats du segment précédent (pour calculs de condensation)

    """

    if T_e is None:
        T_e = T_in


    #################################Ss
    # VALEUR PROVISOIRS CALCULS
    ##################################

    # Température moyenne des fumées
    T_m = Q_(223, 'degC')

    if avec_condensation:

        # Température fumée sortie conduit équilibre
        T_ob = Q_(100, 'degC')
        T_ob_prev = None

        # Température de la paroi intérieure à la sortie, équilibre
        T_iob = Q_(50, 'degC')
        T_iob_prev = None

        if data_j_prev:
            T_ob_prev = data_j_prev['T_ob']
            T_iob_prev = data_j_prev['T_iob']

    ###################################

    run_while_cycle = True
    must_repeat = True
    is_condensing = data_j_prev['condensation'] if data_j_prev else False

    while run_while_cycle:

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
        # show_latex(sigma_CO2_f.__doc__)
        # show_value(sigma_CO2)

        # Coefficient de calcul de la teneur en vapeur d'eau des fumées
        # Annexe B3, valeurs pour gaz et tirage forcé
        f_w = 56
        # [56, 57]  # %

        f_m1 = [3.72, 3.75]
        f_m2 = [0.053, 0.054]

        # Teneur en vapeur d'eau des fumées (%)
        sigma_H2O = sigma_H2O_f(f_w=f_w, sigma_CO2=sigma_CO2)
        # show_latex(sigma_H2O_f.__doc__)

        # Constant de gaz de l'air
        R_L = 288 * ureg.joule / (ureg.kg * ureg.kelvin)
        R = R_f(R_L=R_L, sigma_CO2=sigma_CO2, avec_condensation=avec_condensation)

        # Hauteur au dessus du niveau de la mer
        z = 10 * ureg.m

        # Température de l'air extérieur
        T_L = Q_(-15, 'degC')
        # Prendre +15 °C pour tirage minimal ou pressio positive maximale ou niveau de l'admisson des fumées dans le conduits

        # Pression de l'air extérieur
        p_L = P_L_f(T_L=T_L, z=z, R_L=R_L)
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

        # Surface totale raccordement
        A_v = U * L_v  # m^2

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
        # show_latex(lambda_A_f.__doc__)
        # show_value(lambda_A)




        # show_value(R.to('joule/(kg*kelvin)'))

        # Densité moyenne des fumées
        rho_m = p_L / (R * T_m.to('kelvin'))
        # show_value(rho_m.to('kg/m^3'))

        # Vitesse moyenne des fumées
        w_m = m_dot / (rho_m * A)
        # show_value(w_m.to('m/s'))

        # Nombre de Reynolds
        # Cette formule est équivalente à Re = (rho_m * w_m * D_h) / eta_A
        Re = (m_dot * D_h) / (A * eta_A)
        # si < 2300 alors Re = 2300
        Re = max(Re, 2300)
        # show_value(Re, 'Re')

        f_c0 = 23  # 23.5
        f_c1 = 0.015
        f_c2 = -0.000007
        f_c3 = 0.0142  # 0.0144

        # Capacité calorifique specifique des fumées (pg 99)
        c_p = c_p_f(t_m=T_m, sigma_CO2=sigma_CO2, f_c0=f_c0, f_c1=f_c1, f_c2=f_c2, f_c3=f_c3)
        # show_latex(c_p_f.__doc__)
        # show_value(c_p, 'c_p')

        # Nombre de Prandtl
        Pr = (c_p * eta_A) / lambda_A
        # show_value(Pr, 'Pr')
        if Pr < 0.6 or Pr > 1.5:
            # TODO: Add an exception
            # show_latex("**Attention**: Pr n'est pas entre 0.6 et 1.5, les calculs peuvent être incorrects.")
            print("**Attention**: Pr n'est pas entre 0.6 et 1.5, les calculs peuvent être incorrects.")

        # valeur moyenne de rugosité tab. B.4 pg 102
        r = 0.001 * ureg.m

        # Coefficient de frottement
        psi = colebrook_psi(Re=Re, r=r, D_h=D_h)
        psi_smooth = colebrook_psi(Re=Re, r=0 * ureg.m, D_h=D_h)

        # show_value(psi, 'psi')
        # show_value(psi_smooth, 'psi_smooth')

        # Nombre de Nusselt
        Nu = Nu_f(psi=psi, psi_smooth=psi_smooth, Re=Re, Pr=Pr, D_h=D_h, L_tot=L + L_v)
        # show_value(Nu, 'Nu')

        # Coefficient intener de transfert de chaleur par convection
        alpha_i = (Nu * lambda_A) / D_h
        # convert to SI
        alpha_i = alpha_i.to_base_units()
        # show_value(alpha_i, 'alpha_i')


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


        if not avec_condensation or not is_condensing:

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
            T_iob = T_ob + (T_uo - T_ob) * (k_b / alpha_i)
            # show_value(T_iob, 'T_iob')

        p_Do_j = p_Do_j_f(T_iob_j=T_iob)

        # Calcul point de rosée
        condensation = T_iob.to('degC') < t_p.to('degC')
        if not is_condensing and condensation:
            is_condensing = True
            must_repeat = True

        if avec_condensation and is_condensing:
            print("Calculs condensation en cours...")
            f_K = f_K_f(T_e=T_e)
            # show_value(f_K, 'f_K')



            delta_m_D_j = delta_m_D_j_f(m_dot=m_dot, R=R, p_D=p_D, p_L=p_L, p_Do_prev=data_j_prev['p_Do_j'], p_Do_j=p_Do_j, f_K=f_K)

            # ENthalpie de l'eau évaporée pg 88
            r_D = 2_400_000 * ureg.joule / ureg.kg

            # Chaleur de condensation
            q_K_j = delta_m_D_j * r_D

            # l_c_j est lavProportion de surface de condensation dans le segment j
            # TODO: fix below: something is wrong with the logic
            is_first_condensing_segment = not (data_j_prev and data_j_prev.get('condensation', False))
            if not is_first_condensing_segment:
                # Ce n'est pas le premier segment avec condensation
                l_c_j = 1
            else:
                # Premier segment avec condensation
                T_iob_j_prev = data_j_prev['T_iob']
                # La température de rosée à l'entrée du premier segment (sauf que normalement elle ne dépende pas de j, donc je prends T_p)
                T_pe_1 = t_p
                # l_c_j doit être entre 0 et 1. Un crop assure que la condition soit respectée pendant les itérations
                l_c_j = max(0.1,min(1, (T_pe_1.to('degC').magnitude - T_iob.to('degC').magnitude) / (T_iob_j_prev.to('degC').magnitude - T_iob.to('degC').magnitude)))

            # Coefficient de transfert thermique par condensation
            alpha_ioK_j = q_K_j / ( l_c_j * (T_ob - T_iob) * U * L)

            # Coefficient tolla de transfert thermique, qui est donc l'effet combiné de la convection et de la condensation (pg. 86)
            alpha_iotot_j = alpha_i + alpha_ioK_j

            k_obtot_j = 1 / (1 / alpha_iotot_j + Res_therm + 1 / alpha_a)

            # Debit massique des fumées à la sortie du segment

            m_dot_o_j = data_j_prev['m_dot_o_j'] - delta_m_D_j

            ## AJOUTER ICI !!!!

            if is_first_condensing_segment:
                k_b_prev = data_j_prev['k_b'] if data_j_prev else k_b
                k_btot_j = (1 - l_c_j) * k_b_prev + l_c_j * k_obtot_j
            else:
                k_btot_prev = data_j_prev['k_btot_j']
                k_btot_j = (k_btot_prev + k_obtot_j) / 2

            # Calcul déjà fait, qui va converger
            c_po_j = c_p
            K_b = (U * k_btot_j * L) / (m_dot_o_j * c_po_j)
            K_b = K_b.to('dimensionless')

            T_ob_j = T_ob_j_f(T_ob_prev=data_j_prev['T_ob'],
                                K_b_j=K_b,
                                q_K_j=q_K_j,
                                m_dot_o_j=m_dot_o_j,
                                m_dot_o_prev=data_j_prev['m_dot_o_j'],
                                c_po_j=c_po_j,
                                c_po_j_prev=data_j_prev['c_po_j'],
                                T_u_j=T_u
                    )
            T_ob = T_ob_j

            T_m_prev = T_m
            T_m = (T_ob.to('kelvin') + T_in.to('kelvin')) / 2


            T_iob = T_ob + (T_uo - T_ob) * (k_obtot_j / alpha_iotot_j)

            # Nous sommes dans un segment avec condensation, donc T_iob < t_p
            T_iob = min(T_iob, t_p - Q_(0.1, 'delta_degC')).to('kelvin')

            q_A = data_j_prev['m_dot_o_j'] * data_j_prev['c_po_j'] * data_j_prev['T_ob'] - m_dot_o_j * c_po_j * T_ob_j + q_K_j

            q_C = alpha_iotot_j * U * L * ( data_j_prev['T_ob'] - data_j_prev['T_iob'] + T_ob_j - T_iob ) / 2

            delta_iter_cond = math.fabs(q_A.to('watt').magnitude - q_C.to('watt').magnitude) / q_A.to('watt').magnitude
            print(f'{delta_iter_cond=}')

        # Critère de convergences
        delta_iters = [T_m - T_m_prev]
        if not must_repeat and avec_condensation:
            delta_iters += [T_ob - T_ob_prev, T_iob - T_iob_prev]
        # Coeff de contraction pour s'assurer la convergence des itérations
        a = 1


        # T_ob = (1-a) * T_ob_prev + a * T_ob
        # T_iob = (1-a) * T_iob_prev + a * T_iob
        run_while_cycle = any(abs(d.magnitude) > 1 for d in delta_iters)

        T_ob_prev = T_ob
        T_iob_prev = T_iob

        # run_while_cycle = any(abs(d.magnitude) > 1 for d in delta_iters)
        run_while_cycle = run_while_cycle or must_repeat
        must_repeat = False

    res = {
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
        }

    if avec_condensation and is_condensing:
        res.update({
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
            })

    # convert everything to base units
    for key in res:
        if isinstance(res[key], pint.Quantity):
            res[key] = res[key].to_base_units()

    return res



# Température fumée sortie foyer)



# Diamètre cheminée/appareil à combustion
D_h = 200 * ureg.mm



# Débit massique (fourni)
dm = 0.008 * ureg.kg / ureg.s  # kg/s

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



# Longueur du conduit=

# Circuit de raccordement
# L = 0.2 * ureg.m
# Conduit de fumée interieur bâtiment
# L = 6.23 * ureg.m
# Conduit de fumée exterieur bâtiment
# L = 1.2 * ureg.m

# les segments doivent avoir la même longueur <= 0.5 m (pg 82 8.1)
# on va faire 21 × 0.3 m + 4 × 0.3 m


# Pg 20
# la valeur doit être connue, ou en alternative il faut prendre les resultats des formules de l'annexe B
T_e = Q_(300, 'degC')

# Coefficient externe de transfert thermique exterieur
alpha_a_ext = 25 * ureg.watt / (ureg.meter ** 2 * ureg.kelvin)
alpha_a_int = 8 * ureg.watt / (ureg.meter ** 2 * ureg.kelvin)

# Calculs circuit de raccordement

# Nombre de segments dans le circuit de raccordement
NsegV = 1

print("Circuit de raccordement")
L_v = 0.2 * ureg.m
T_in = T_e
T_uo = T_u = Q_(15, 'degC') # raccordement dans le batiment
Res_therm = 0 / (ureg.watt / (ureg.meter ** 2 * ureg.kelvin))
res_racc = solve_conduit(T_in=T_in, T_u=T_uv, T_uo=T_uo, L=L_v, m_dot=dm, D_h=D_h, Res_therm=Res_therm, alpha_a=alpha_a_int)
print(json.dumps(res_racc, indent=2, default=str))

# ==============================
# Calculs conduit de fumée
# ==============================

# Nombre de segments dans le conduit de fumée
Nseg = 21 + 4
NSegK = math.inf

L_cum = 0 * ureg.m
L_tot = 6.3 * ureg.m
Nseg_batiment = 21
# Dans le batiment
T_in = res_racc['T_ob']
for i in range(Nseg_batiment):
    L = L_tot / Nseg_batiment
    L_cum += L
    print(f"Tronçon intérieur du bâtiment {i+1}/21, longueur cumulée {L_cum.to('m')}")
    T_uo = T_u = Q_(15, 'degC') # conduit dans le batiment
    Res_therm = 0 / (ureg.watt / (ureg.meter ** 2 * ureg.kelvin))
    res = solve_conduit(T_in=T_in, T_u=T_u, T_uo=T_uo, L=L, m_dot=dm, D_h=D_h, Res_therm=Res_therm, alpha_a=alpha_a_int, avec_condensation=True, data_j_prev=res if i > 0 else None)
    T_in = res['T_ob']
    print(json.dumps(res, indent=2, default=str))
    if res['condensation']:
        print(f"**Attention**: Condensation détectée dans le segment {i+1}/21 du conduit intérieur du bâtiment")
        NSegK = min(i+1, NSegK)

print("=============================")
print("Sortie du batiment")
print("=============================")
# Exterieur batiment
for i in range(4):
    print(f"Tronçon extérieur du bâtiment {i+1}/4")
    L = 0.3 * ureg.m
    T_in = res['T_ob']
    T_uo = T_u = Q_(-15, 'degC') # conduit à l'exterieur
    Res_therm = 0 / (ureg.watt / (ureg.meter ** 2 * ureg.kelvin))
    res = solve_conduit(T_in=T_in, T_u=T_u, T_uo=T_uo, L=L, m_dot=dm, D_h=D_h, Res_therm=Res_therm, alpha_a=alpha_a_ext, avec_condensation=True, data_j_prev=res)
    print(json.dumps(res, indent=2, default=str))

print("Done!")