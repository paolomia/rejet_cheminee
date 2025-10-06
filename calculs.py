from report import *
import math
import json

def p_H_f(H, rho_L, rho_m):
    r"""
    p_H = H \cdot g \cdot (\rho_L - \rho_m)
    """
    # (31) pg 34
    # Tirage théorique du au tirage naturel
    g = 9.81 * ureg.m / ureg.s ** 2
    return H * g * (rho_L - rho_m)

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



def do_report(res: dict, filename: str = "rapport_nf_en_13384.md") -> str:
    """
    Build a Markdown+LaTeX report of the conduit calculation and save it to disk.
    The report includes normative references (EN 13384-1:2015+A1:2019).
    Returns the path to the written file.
    """
    import datetime
    import math

    # ---------- small helpers ----------
    def _is_qty(x):
        return hasattr(x, "to") and hasattr(x, "magnitude")

    def fmt(x, unit=None, digits=3):
        """
        Format pint quantities or plain numbers using SI units and up to `digits`
        significant figures. Uses scientific notation if |value| > 1e4 or 0 < |value| < 1e-3.
        """
        import math

        def _is_qty(v):
            return hasattr(v, "to") and hasattr(v, "magnitude") and hasattr(v, "units")

        def _format_sig(v: float, sig: int = 3) -> str:
            if not math.isfinite(v):
                return str(v)
            if v == 0:
                return "0"
            av = abs(v)
            use_sci = (av > 1e4) or (av < 1e-3)
            if use_sci:
                # sci notation: sig digits total -> (sig-1) after decimal
                return f"{v:.{max(sig - 1, 0)}e}"
            # fixed notation with sig digits
            exp10 = math.floor(math.log10(av))
            decimals = sig - 1 - exp10  # can be negative (round to tens, hundreds, ...)
            rounded = round(v, decimals)
            s = f"{rounded:.{max(decimals, 0)}f}"
            # normalize and strip trailing zeros
            if s.startswith("-0"):
                s = s.replace("-0", "0", 1)
            if "." in s:
                s = s.rstrip("0").rstrip(".")
            return s

        def _to_si(q):
            # If an explicit target unit was requested, try it first.
            if unit is not None:
                try:
                    return q.to(unit)
                except Exception:
                    pass  # fall back to SI suggestions

            # Try common SI targets first (then fall back to base SI).
            si_candidates = [# scalars & basics
                "1", "K", "m", "m**2", "m**3", "kg", "s", # flows, densities, velocities
                "kg/s", "kg/m**3", "m/s", "m/s**2", # thermo/fluids
                "Pa", "W", "J", "N", "W/(m*K)", "W/(m**2*K)", "Pa*s", "J/(kg*K)", "m*K/W", "m**2*K/W", ]
            for ustr in si_candidates:
                try:
                    return q.to(ustr)
                except Exception:
                    continue
            # Last resort: base SI units
            try:
                return q.to_base_units()
            except Exception:
                return q

        try:
            if _is_qty(x):
                q = _to_si(x)
                val = float(q.magnitude)
                s = _format_sig(val, digits)
                # prefer compact SI unit symbols
                try:
                    u_str = "{:~}".format(q.units)  # e.g., 'Pa', 'W/(m*K)', 'm'
                except Exception:
                    try:
                        u_str = format(q.units, "~")
                    except Exception:
                        u_str = str(q.units)
                if u_str in ("dimensionless", ""):
                    return s
                return f"{s} {u_str}"
            elif isinstance(x, (int, float)):
                return _format_sig(float(x), digits)
            else:
                return str(x)
        except Exception:
            return str(x)


    # bring selected inputs if they exist in globals()
    # (kept light: report is robust even if some are missing)
    # --- Inputs (names and symbols per EN 13384-1) ---
    ALIASES = {
        # Températures principales
        "$T_L$ (température de l’air extérieur)": "T_L",
        "$T_e$ (température des fumées à l’orifice d’admission du conduit)": "T_e",
        "$T_uo$ (température de l’air ambiant à la sortie du conduit de fumée)": "T_uo",
        "$T_u$ (température de l’air ambiant autour du conduit)": "T_u",

        # Géométrie et débits
        "$Dh$ (diamètre hydraulique intérieur du conduit)": "D_h",
        "$L$ (longueur du conduit de fumée)": "L",

        # Propriétés / paramètres
        "$r$ (rugosité moyenne intérieure)": "r",
        "$1/Λ$ (résistance thermique du conduit)": "Res_therm",

        # Puissances et coefficients (hors tableau normatif mais nécessaires)
        "$Q_F$ (puissance thermique nominale)": "Q_F",
        "$Q_N$ (puissance utile nominale)": "Q_N",
        "$f_w$ (teneur en vapeur d’eau des fumées)": "f_w",
        "$f_x1$": "f_x1",
        "$f_x2$": "f_x2",
        "$f_x3$": "f_x3",
        "$\eta$ (rendement)": "eta",
        }
    inputs = []
    for label, var in ALIASES.items():
        if var in res:
            inputs.append((label, res[var]))

    # convenience locals from results (if present)
    Re = res.get("Re")
    Pr = res.get("Pr")
    Nu = res.get("Nu")
    Tiob = res.get("T_iob")
    Tob = res.get("T_ob")
    Tm = res.get("T_m")

    # ---------- domain checks per standard ----------
    checks = []
    # Nu validity window: 2300 < Re < 1e7 and 0.6 < Pr < 1.5 (EN 13384-1 §5.8.3.2, eq. (24))
    if Re is not None:
        try:
            Re_v = float(Re.to('dimensionless').magnitude) if _is_qty(Re) else float(Re)
            ok_re = (Re_v >= 2300.0) and (Re_v <= 1.0e7)
            checks.append(f"- Domaine de validité pour Nu (eq. (24)) — 2.3e3 ≤ Re ≤ 1.0e7 : "
                          f"{'OK' if ok_re else 'HORS DOMAINE'} (Re = {Re_v:.3g})")
        except Exception:
            checks.append("- Domaine Re : impossible d’évaluer.")
    if Pr is not None:
        try:
            Pr_v = float(Pr.to('dimensionless').magnitude) if _is_qty(Pr) else float(Pr)
            ok_pr = (Pr_v > 0.6) and (Pr_v < 1.5)
            checks.append(f"- Domaine de validité pour Nu (eq. (24)) — 0.6 < Pr < 1.5 : "
                          f"{'OK' if ok_pr else 'HORS DOMAINE'} (Pr = {Pr_v:.3g})")
        except Exception:
            checks.append("- Domaine Pr : impossible d’évaluer.")

    # ---------- build Markdown ----------
    md = []
    md.append("# Rapport de calcul — Conduit de fumée")
    md.append("")
    md.append("_Selon NF EN 13384‑1:2015+A1:2019 (méthodes de calcul thermo‑aéraulique, partie 1)_")
    md.append("")
    md.append(f"**Date** : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    md.append("")
    md.append("## 1. Données d’entrée")
    if inputs:
        for k, v in inputs:
            md.append(f"- **{k}** = {fmt(v)}")
    else:
        md.append("_(Aucune donnée d’entrée supplémentaire détectée dans le module.)_")

    md.append("")
    md.append("## 2. Formules utilisées et références normatives")
    md.append("")
    md.append("**Thermo‑aéraulique principale**")
    md.append("")
    md.append("- Pression atmosphérique extérieure :")
    md.append("  $$ p_L = 97000\\,\\exp\\!\\left(\\tfrac{-g\\,z}{R_L\\,T_L}\\right) \\quad\\text{(eq. (12), §5.7.2)} $$")
    md.append("- Températures le long du conduit (sortie et moyenne) :")
    md.append("  $$ T_m = T_u + \\frac{T_e - T_u}{K}\\,\\bigl(1-e^{-K}\\bigr) \\quad\\text{(eq. (16), §5.8.1)} $$")
    md.append("  $$ T_o = T_u + (T_e - T_u)\\,e^{-K} \\quad\\text{(eq. (17), §5.8.1)} $$")
    md.append("- Coefficient de refroidissement :")
    md.append("  $$ K = \\frac{k\\,U\\,L}{\\dot{m}\\,c_p} \\quad\\text{(eq. (20), §5.8.2)} $$")
    md.append("- Coefficients de transfert thermique :")
    md.append("  $$ k_b = \\left(\\tfrac{1}{\\alpha_i} + \\tfrac{1}{\\Lambda} + \\tfrac{1}{\\alpha_a}\\right)^{-1} "
              "\\quad\\text{(eq. (21), §5.8.3.1)} $$")
    md.append("  $$ k = \\left(\\tfrac{1}{\\alpha_i} + \\tfrac{1}{\\Lambda} + \\tfrac{1}{\\alpha_a}\\right)^{-1} "
              "\\text{ avec } S_H\\ \\text{selon §5.7.7 (eq. (22), §5.8.3.1)} $$")
    md.append("- Transfert interne par convection :")
    md.append("  $$ \\alpha_i = \\frac{\\mathrm{Nu}\\,\\lambda_A}{D_h} \\quad\\text{(eq. (23), §5.8.3.2)} $$")
    md.append("- Nombre de Nusselt :")
    md.append("  $$ \\mathrm{Nu} = \\left(\\tfrac{\\psi}{\\psi_{\\text{smooth}}}\\right)^{0.67}\\,0.0214\\,"
              "\\bigl(\\mathrm{Re}^{0.8}-100\\bigr)\\,\\mathrm{Pr}^{0.4}\\,"
              "\\left[1+\\left(\\tfrac{D_h}{L_{tot}}\\right)^{0.67}\\right] \\quad\\text{(eq. (24), §5.8.3.2)} $$")
    md.append("- Nombres sans dimension :")
    md.append("  $$ \\mathrm{Pr} = \\tfrac{c_p\\,\\eta_A}{\\lambda_A} \\quad\\text{(eq. (25), §5.8.3.2)} "
              "\\qquad \\mathrm{Re} = \\tfrac{\\rho_m\\,w_m\\,D_h}{\\eta_A} \\quad\\text{(eq. (26), §5.8.3.2)} $$")
    md.append("  _(Forme utilisée dans le code pour Re : "
              "$\\mathrm{Re} = \\tfrac{\\dot{m}\\,D_h}{A\\,\\eta_A}$, équivalente à l’eq. (26).)_")
    md.append("- Grandeurs de base :")
    md.append("  $$ \\rho_m = \\tfrac{p_L}{R\\,T_m} \\quad\\text{(eq. (27), §5.9.1)} \\qquad "
              "w_m = \\tfrac{\\dot{m}}{\\rho_m\\,A} \\quad\\text{(eq. (28), §5.9.2)} $$")
    md.append("- Frottement (équation de Colebrook sous la forme de la norme) :")
    md.append("  $$ \\frac{1}{\\sqrt{\\psi}} = -2\\,\\log_{10}\\!\\left(\\frac{2.51}{\\mathrm{Re}\\,\\sqrt{\\psi}} "
              "+ \\frac{r}{3.71\\,D_h}\\right) \\quad\\text{(eq. (35), §5.10.3.3)} $$")
    md.append("")
    md.append("**Propriétés des fumées (Annexe B)**")
    md.append("")
    md.append("- $R$ (constante des gaz fumées), $c_p(T_m,\\sigma\\_{CO_2})$, $\\lambda_A(T_m)$ : "
              "voir §5.7.3.2, §5.7.5 et §5.8.3.2 (formules données en Annexe B, tab. B.1/B.3/B.8).")
    md.append("- $\\eta_A(T_m)$ (viscosité dynamique) : formule (B.10), Annexe B.")
    md.append("- $\\sigma\\_{CO_2}$, $\\sigma\\_{H_2O}$ : valeurs de l’Annexe B selon le combustible/appareil.")
    md.append("- Température de paroi en sortie $T_{iob}$ : voir §5.12 (méthode de calcul).")

    md.append("")
    md.append("## 3. Résultats principaux")
    # show a curated subset first
    main_keys = [("A", "$A$ (Section interne)"), ("A_t", "$A_t$ (Surface totale du conduit de fumée)"), ("A_v", "$A_v$ (Surface totale du conduit de raccordement)"), ("p_L", "$p_L$ (Pression extérieure)"), ("R", "$R$ (Constante des gaz des fumées)"), ("c_p", "$c_p$ (Capacité calorifique spécifique des fumées)"), ("lambda_A", "$\\lambda_A$ (Conductivité thermique des fumées)"), ("eta_A", "$\\eta_A$ (Viscosité dynamique des fumées)"), ("rho_m", "$\\rho_m$ (Masse volumique moyenne des fumées)"), ("w_m", "$w_m$ (Vitesse moyenne des fumées)"), ("Re", "$\\mathrm{Re}$ (Nombre de Reynolds)"), ("Pr", "$\\mathrm{Pr}$ (Nombre de Prandtl)"), ("psi", "$\\psi$ (Coefficient de frottement)"), ("psi_smooth", "$\\psi_{smooth}$ (Coefficient de frottement pour conduit lisse)"), ("Nu", "$\\mathrm{Nu}$ (Nombre de Nusselt)"), ("alpha_i", "$\\alpha_i$ (Coefficient convectif interne)"),
        ("Res_therm", "$1/\\Lambda$ (Résistance thermique du conduit)"), ("k", "$k$ (Coefficient de transfert thermique, non équilibré)"), ("k_b", "$k_b$ (Coefficient de transfert thermique, équilibré)"), ("K", "$K$ (Coefficient de refroidissement, non équilibré)"), ("K_b", "$K_b$ (Coefficient de refroidissement, équilibré)"), ("K_ob", "$K_{ob}$ (Coefficient de refroidissement à la sortie)"), ("T_m", "$T_m$ (Température moyenne des fumées)"), ("T_ob", "$T_{ob}$ (Température des fumées à la sortie, équilibrée)"), ("T_iob", "$T_{iob}$ (Température de la paroi interne à la sortie)"), ("m_dot", "$\\dot{m}$ (Débit massique des fumées)"), ]
    for key, label in main_keys:
        if key in res:
            md.append(f"- **{label}** = {fmt(res[key])}")

    md.append("")
    md.append("## 4. Vérifications de domaine (conformité des équations)")
    if checks:
        md.extend(checks)
    else:
        md.append("_Aucune vérification disponible._")

    md.append("")
    md.append("## 5. Remarques et cas particuliers")
    md.append("- **αa (extérieur/intérieur)** : la norme fixe 23 W/(m²·K) à l’extérieur et 8 W/(m²·K) à "
              "l’intérieur (§5.8.3.3). Nous utilisons une interpolation linéaire entre ces deux valeurs, en fonction de la proportion de cheminée"
                " en intérieur/extérieur (cf. §5.8.3.3).")
    md.append("- **T_u = 15 °C** car environnement chauffé avec surface de cheminée située à l'exterieur < 1/4 de la surface totale (§ 5.7.1.3 NOTE 1)")
    md.append("- Les résultats obtenus sont calculés par approximation itérative, avec une valeur de T_m initiale de 223 °C à la première itération.")
    md.append('')
    md.append("---")


    # ---------- write file ----------
    content = "\n".join(md) + "\n"
    with open(filename, "w", encoding="utf-8") as f:
        f.write(content)

    # Provide a console hint for the caller (kept minimal)
    print(f"Markdown report written to: {filename}")
    return filename

def m_dot_f(Q_F, f_m1, f_m2, sigma_CO2):
    r"""
    \dot{m} = \left(\frac{f_{m1}}{\sigma(\text{CO}_2)} + f_{m2}\right) \cdot Q_F
    """
    res = (f_m1 / sigma_CO2 + f_m2) * Q_F.to("kW").magnitude
    return Q_(res/1000, "kg/s")



def solve_conduit(T_L, T_e, T_u, T_uo, L, D_h, Res_therm, L_ext, Q, eta, z, H, type_calcul):
    """
    T_L: Température de l'air extérieur (K)
    T_e: Température des fumées à l’orifice d’admission du conduit (K)
    T_u: Température de l'air ambiant autour du conduit (K)
    T_uo: Température de l'air ambiant à la sortie du conduit de fumée (K)
    L: Longueur du conduit de fumée (m)
    D_h: Diamètre hydraulique intérieur du conduit (m)
    Res_therm: Résistance thermique du conduit (m²·K/W)
    L_ext: Longueur extérieure du conduit (m)
    Q: Puissance utile nominale (kW)
    eta: Rendement de l'appareil à combustion (%)
    z: Altitude (m)
    H: Hauteur totale du conduit (m)
    type_calcul: Type de calcul (str) - "pression_maxi" pour le calcul sur le tirage minimal,
        ou "pression_mini" pour le calcul sur le tirage maximal et les exigences de température de paroi.
    """


    #################################Ss
    # VALEUR PROVISOIRS CALCULS
    ##################################

    # Température moyenne des fumées
    T_m = Q_(223, 'degC')
    T_m_prev = None

    ###################################

    # Choix des coefficients selon le type de calcul
    if type_calcul == "pression_maxi":
        S_H = 0.5    # 5.7.7 pg 28
        S_E = 1.2
    elif type_calcul == "pression_mini":
        S_H = 1
        S_E = 1

    while (T_m_prev is None or abs((T_m - T_m_prev)).magnitude > 1e-1):
        T_m_prev = T_m
        # Circonférence intérieure du conduit de fumée
        U = 3.1416 * D_h  # m

        # Section cheminée/appareil à combustion
        A = 3.1416 * (D_h / 2) ** 2  # m^2
        # show_value(A)

        # Surface totale conduit de fumée
        A_t = U * L  # m^2
        # show_value(A_t)

        # Coefficient de calcul de la teneur en vapeur d'eau des fumées
        # Annexe B3, valeurs pour gaz et tirage forcé
        f_w = 56
        # [56, 57]  # %


        # Puissance thermique nominale
        Q_F = Q / eta

        # 48 - 195

        # Annexe B3
        f_x1 = 8.6
        f_x2 = 0.078
        f_x3 = 10.2

        # Puissance utile nominale
        Q_N = Q


        # Coefficient de calcul de la teneur en CO2 des fumées
        sigma_CO2 = sigma_CO2_f(Q_N=Q_N, f_x1=f_x1, f_x2=f_x2, f_x3=f_x3)
        # show_latex(sigma_CO2_f.__doc__)
        # show_value(sigma_CO2)

        # rendement de l'appareil à combustion (%)
        # eta = eta_w(Q_N=Q_N)
        # show_latex(eta_w.__doc__)
        # show_value(eta)

        # On calcul le debit massique des fumées à avec B.1
        f_m1 = 3.735
        f_m2 = 0.0535
        m_dot = m_dot_f(Q_F, f_m1, f_m2, sigma_CO2)


        # Teneur en vapeur d'eau des fumées (%)
        sigma_H2O = sigma_H2O_f(f_w=f_w, sigma_CO2=sigma_CO2)
        # show_latex(sigma_H2O_f.__doc__)

        # Constant de gaz de l'air
        R_L = 288 * ureg.joule / (ureg.kg * ureg.kelvin)

        # Pression de l'air extérieur
        p_L = p_L_f(T_L=T_L, z=z, R_L=R_L)
        # show_latex(P_L_f.__doc__)
        # show_value(p_L)



        # Viscosité dynamique des fumées
        eta_A = eta_A_f(t_m=T_m)
        # show_latex(eta_A_f.__doc__)
        # show_value(eta_A, 'eta_A')

        # Coefficient de conductivité thermique des fumées
        lambda_A = lambda_A_f(t_m=T_m)
        # show_latex(lambda_A_f.__doc__)
        # show_value(lambda_A)

        # coefficient de calcul de la constante des gaz des fumÈes, en 1/%
        # pg 96
        f_R = 0.0033
        # ou 0.0032

        R = R_f(f_R=f_R, R_L=R_L, sigma_CO2=sigma_CO2)
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
        show_value(Re, 'Re')

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
        # if Pr < 0.6 or Pr > 1.5:
        #     show_latex("**Attention**: Pr n'est pas entre 0.6 et 1.5, les calculs peuvent être incorrects.")

        # valeur moyenne de rugosité tab. B.4 pg 102
        r = 0.001 * ureg.m

        # Coefficient de frottement
        psi = colebrook_psi(Re=Re, r=r, D_h=D_h)
        psi_smooth = colebrook_psi(Re=Re, r=0 * ureg.m, D_h=D_h)

        # show_value(psi, 'psi')
        # show_value(psi_smooth, 'psi_smooth')

        # Nombre de Nusselt
        Nu = Nu_f(psi=psi, psi_smooth=psi_smooth, Re=Re, Pr=Pr, D_h=D_h, L_tot=L)
        # show_value(Nu, 'Nu')

        # Coefficient intener de transfert de chaleur par convection
        alpha_i = (Nu * lambda_A) / D_h
        # show_value(alpha_i, 'alpha_i')



        # Coefficient externe de transfert thermique exterieur
        alpha_a_ext = 25 * ureg.watt / (ureg.meter ** 2 * ureg.kelvin)
        alpha_a_int = 8 * ureg.watt / (ureg.meter ** 2 * ureg.kelvin)
        # on prend une interpolation
        alpha_e = ((L - L_ext) * alpha_a_int + L_ext * alpha_a_ext) / L

        # Coefficient de transfert thermique température NON équilibrée
        k = 1 / (1 / alpha_i + S_H * (Res_therm + 1 / alpha_e))
        # show_value(k, 'k')

        # Coefficient de transfert thermique température équilibrée
        k_b = 1 / (1 / alpha_i + Res_therm + 1 / alpha_e)
        # show_value(k_b, 'k_b')

        # Coefficient de refroidissement température NON équilibrée
        K = (k * U * L) / (m_dot * c_p)
        # show_value(K, 'K')
        # Coefficient de refroidissement température équilibrée
        K_b = (k_b * U * L) / (m_dot * c_p)

        # Coeffiecient de refroidissement température à la sortie
        k_ob = 1 / (1 / alpha_i + Res_therm + 1 / alpha_a_ext)

        # Température des fumées à la sortie  équilibrée
        T_ob = T_u + (T_e - T_u) * math.exp(-K_b)
        # show_value(T_ob, 'T_ob')

        # Température moyenne des fumées
        T_m = T_u + (T_e - T_u) * (1 - math.exp(-K)) / K
        # show_value(T_m, 'T_m')

        # Température de la paroi intérieure à la sortie
        T_iob = T_ob + (T_uo - T_ob) * (k_ob / alpha_i)
        # show_value(T_iob, 'T_iob')

        rho_1 = p_L / (R * T_e.to('kelvin'))
        w_1 = m_dot / (rho_1 * A)
        rho_2 = p_L / (R * T_u.to('kelvin'))
        w_2 = m_dot / (rho_2 * A)

        # Température de rosée
        p_D = sigma_H2O / 100 * p_L
        t_p = 4077.9 / (23.6448 - math.log(p_D.magnitude)) - 236.67
        t_p = Q_(t_p, 'degC')

    # Vérification exigeance de température minimale de la paroi interne à la sortie (T_iob) - (6)
    condensation = T_iob < t_p

    # CALCULS DE PRESSIONS
    rho_L = p_L / (R_L * T_L.to('kelvin'))
    # Tirage naturel théorique
    p_H = p_H_f(H, rho_L, rho_m).to('pascal')

    # Difference de pression du au changement de vitesse
    P_G = 0.5 * (rho_2 * w_2 ** 2 - rho_1 * w_1 ** 2).to('pascal')

    # Perte de charge par frottement
    P_R_frottement = (psi * (L / D_h) * 0.5 * rho_m * w_m ** 2).to('pascal')

    S_EG = S_E if P_G.magnitude >= 0 else 1

    # Perte de charge du au raccordement
    # pg 106 avec m_dot_1 = 0
    Zeta_racc = 1.2

    # Nous avons une seule changement de direction/raccord, entre le circuit de raccordement et le circuit de fumée
    sum_Zeta_i = Zeta_racc

    P_R_raccord = (sum_Zeta_i * 0.5 * rho_m * w_m ** 2).to('pascal')
    P_E = P_R_frottement + P_R_raccord
    P_R = S_E * P_E + S_EG * P_G



    return {
        'A': A,
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
        'c_p': c_p.to('joule/(kg*kelvin)'),
        'psi': psi,
        'psi_smooth': psi_smooth,
        'Nu': Nu,
        'alpha_i': alpha_i.to('watt/(m^2*kelvin)'),
        'Res_therm': Res_therm.to('m^2*K/watt'),
        'k': k,
        'k_b': k_b,
        'K': K,
        'K_b': K_b,
        'T_ob': T_ob.to('kelvin'),
        'T_iob': T_iob.to('kelvin'),
        'z': z.to('m'),
        'rho_1': rho_1.to('kg/m^3'),
        'w_1': w_1.to('m/s'),
        'rho_2': rho_2,
        'w_2': w_2.to('m/s'),
        'm_dot': m_dot,
        't_p': t_p,
        'H': L,  # hauteur totale
        'D_h': D_h,
        'T_e': T_e,
        'T_u': T_u,
        'T_uo': T_uo,
        'T_L': T_uo,
        'Q_N': Q_N,
        'condensation': condensation,
        'P_H': p_H,
        'P_G': P_G,
        'P_E': P_E,
        'P_R': P_R,
        }


def main():
    from json_pint import json_pint_load, json_pint_dump

    # ╔╦╗╔═╗╔╗╔╔╗╔╔═╗╔═╗╔═╗  ╔═╗╔╗╔╔╦╗╦═╗╔═╗╔═╗
    #  ║║║ ║║║║║║║║╣ ║╣ ╚═╗  ║╣ ║║║ ║ ╠╦╝║╣ ║╣
    # ═╩╝╚═╝╝╚╝╝╚╝╚═╝╚═╝╚═╝  ╚═╝╝╚╝ ╩ ╩╚═╚═╝╚═╝
    # ==============================================================
    # LIRE ATTENTIVEMENT LES DONNÉES CI-DESSOUS ET LES MODIFIER SELON VOTRE CAS

    # --- Données thermiques de l’appareil ---

    # Température des fumées à la sortie du générateur (Te)
    T_e = Q_(310, "degC")

    # Puissance utile nominale (Q̇N) — selon la fiche technique du générateur
    Q_N = Q_(140, "kW")

    # Rendement de l’appareil à combustion (η)
    # Valeur comprise en général entre 0.8 et 0.9 pour les chaudières gaz.
    eta = 0.86

    # --- Géométrie du conduit de fumée ---

    # Diamètre hydraulique intérieur du conduit (Dh)
    D_h = 200 * ureg.mm

    # Longueur totale du conduit de fumée (L)
    # Inclut la section de raccordement et la cheminée
    L = 0.2 * ureg.m + 7.5 * ureg.m

    # Difference d'hauteur entre l'orifice d'admission du conduit et la sortie du conduit (H)
    H = L - 0.2 * ureg.m  # j'enleve le raccord

    # Longueur de la partie extérieure du conduit (L_ext)
    L_ext = 1.2 * ureg.m

    # Résistance thermique du conduit (1/Λ) - 1 mm de tôle simple c'est négligeable
    Res_therm = 0 / (ureg.watt / (ureg.meter ** 2 * ureg.kelvin))

    # --- Conditions environnementales ---

    # Température de l’air ambiant autour de la cheminée (Tu)
    # On prend 15 °C si la surface de cheminée située à l'extérieur est < 1/4 de la surface totale (§ 5.7.1.3 NOTE 1)
    T_u = Q_(15, "degC")

    # Hauteur au-dessus du niveau de la mer (z)
    z = 41 * ureg.m

    # Pression de la vitesse du vent - voir pg 36, 5.10.4
    # P_L = 25 * ureg.pascal
    P_L = 0 * ureg.pascal

    # ==============================================================
    # ╔═╗╔═╗╦  ╔═╗╦ ╦╦  ╔═╗
    # ║  ╠═╣║  ║  ║ ║║  ╚═╗
    # ╚═╝╩ ╩╩═╝╚═╝╚═╝╩═╝╚═╝
    # ==============================================================
    # Vérification des critère de température, tirage maximal / pression positive minimal
    # Nous prenons le cas le plus défavorable dans ce contexte, il fait froid (T_L = -15 °C) donc plus de risque de condensation
    # et le tirage est au maximum, sans le vent, donc la pression à l'entrée du conduit est la plus basse possible (P_ZOmin)
    # ==============================================================
    #
    # Température de l’air extérieur (TL)
    T_L = Q_(-15, "degC")

    # Température de l’air ambiant à la sortie du conduit (Tuo)
    # 0 °C pour ambiance sèche, sans condensation - c'est donc différent de T_L, ce n'est pas une erreur
    T_uo = Q_(0, "degC")
    res = solve_conduit(T_L=T_L, T_e=T_e, T_u=T_u, T_uo=T_uo, L=L, D_h=D_h, Res_therm=Res_therm, L_ext=L_ext, Q=Q_N, eta=eta, z=z, H=H, type_calcul="pression_mini")

    if res['condensation']:
        print("**ATTENTION**: Condensation probable dans le conduit de fumée (T_iob < t_p). Repeter les calculs pour le cas avec condensation (calculs_condensation_temperature.py). "
              "Bon courage !!")
    else:
        print("Pas de condensation dans le conduit de fumée, critère (6) vérifié.")

    P_ZOmin = res['P_R'] - res['P_H']
    print(f"P_ZOmin (pression minimale) = {P_ZOmin:.1f}")

    with open("results_tirage_maxi.json", "w") as f:
        json_pint_dump(res, f, ureg, indent=2)

    # ==============================================================
    # Vérification des critère de tirage minimal / pression positive maximale
    # Il fait chaud et il y a du vent, donc la pression à l'entrée du conduit est la plus haute possible (P_ZO)
    # Il faut que cette valeur soit inférieure à la valeur max admissible, sinon nous avons un retour des fumées
    # ==============================================================
    # Température de l’air extérieur (TL)
    T_L = Q_(15, "degC") # Pas très élevé, mais c'est le valeur donnée dans la norme (5.7.1.3)

    T_uo = T_L
    res = solve_conduit(T_L=T_L, T_e=T_e, T_u=T_u, T_uo=T_uo, L=L, D_h=D_h, Res_therm=Res_therm, L_ext=L_ext, Q=Q_N, eta=eta, z=z, H=H, type_calcul="pression_maxi")

    P_ZO = res['P_R'] - res['P_H'] + P_L

    print(f"P_ZO (pression maximale) = {P_ZO:.1f}")

    with open("results_tirage_minim.json", "w") as f:
        json_pint_dump(res, f, ureg, indent=2)


if __name__ == "__main__":
    main()
