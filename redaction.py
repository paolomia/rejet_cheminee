def do_report(res_min: dict, res_max: dict, filename: str = "/note_calculs/note_de_calcul.md") -> str:
    """
    Génère un rapport Markdown+LaTeX comparatif du calcul de conduit et le sauvegarde.
    Le rapport inclut les références normatives (EN 13384-1:2015+A1:2019) pour
    les hypothèses de tirage minimal et maximal.
    Retourne le chemin du fichier écrit.
    """
    import datetime
    import os
    import math

    # ---------- Petits assistants ----------
    def _is_qty(x):
        return hasattr(x, "to") and hasattr(x, "magnitude")

    def fmt(x, unit=None, digits=3):
        """
        Formate les quantités pint ou les nombres simples en utilisant les unités SI
        et jusqu'à `digits` chiffres significatifs. Utilise la notation scientifique
        si |valeur| > 1e4 ou 0 < |valeur| < 1e-3.
        """

        def _is_qty(v):
            return hasattr(v, "to") and hasattr(v, "magnitude") and hasattr(v, "units")

        def _format_sig(v: float, sig: int = 3) -> str:
            if not math.isfinite(v): return str(v)
            if v == 0: return "0"
            av = abs(v)
            use_sci = (av > 1e4) or (av > 0 and av < 1e-3)
            if use_sci:
                return f"{v:.{max(sig - 1, 0)}e}"
            exp10 = math.floor(math.log10(av)) if av > 0 else 0
            decimals = sig - 1 - exp10
            rounded = round(v, int(decimals))
            s = f"{rounded:.{max(int(decimals), 0)}f}"
            if "." in s: s = s.rstrip("0").rstrip(".")
            return s

        def _to_si(q):
            if unit is not None:
                try:
                    return q.to(unit)
                except Exception:
                    pass
            si_candidates = ["1", "K", "m", "m**2", "m**3", "kg", "s", "kg/s", "kg/m**3", "m/s", "m/s**2", "Pa", "W", "J", "N", "W/(m*K)", "W/(m**2*K)", "Pa*s", "J/(kg*K)", "m*K/W", "m**2*K/W", "degC"]
            for ustr in si_candidates:
                try:
                    return q.to(ustr)
                except Exception:
                    continue
            try:
                return q.to_base_units()
            except Exception:
                return q

        try:
            if _is_qty(x):
                q = _to_si(x)
                val = float(q.magnitude)
                s = _format_sig(val, digits)
                try:
                    u_str = "{:~}".format(q.units)
                except Exception:
                    u_str = str(q.units)
                if u_str in ("dimensionless", ""): return s
                return f"{s} {u_str}"

            elif isinstance(x, bool):
                return "Oui" if x else "Non"
            elif isinstance(x, (int, float)):
                return _format_sig(float(x), digits)
            else:
                return str(x)
        except Exception:
            return str(x)

    # --- Alias pour les tables (libellé, symbole, clé_dict) ---
    ALIASES_INPUTS = [("Température de l’air extérieur", "$T_L$", "T_L"), ("Température des fumées à l’admission", "$T_e$", "T_e"), ("Température de l’air ambiant (conduit)", "$T_u$", "T_u"), ("Température de l’air ambiant (sortie)", "$T_{uo}$", "T_uo"), ("Diamètre hydraulique intérieur", "$D_h$", "D_h"), ("Longueur du conduit de fumée", "$L$", "L"), ("Résistance thermique du conduit", "$1/\\Lambda$", "Res_therm"), ("Puissance utile nominale", "$Q_N$", "Q_N"), ("Rendement", "$\\eta$", "eta"), ("Pression du vent", "$P_L$ (vent)", "P_L_vent"), ]
    ALIASES_RESULTS = [("Tirage naturel théorique", "$P_H$", "P_H"), ("Pertes de charge (frottement + raccords)", "$P_E$", "P_E"), ("Variation de pression (vitesse)", "$P_G$", "P_G"), ("Pertes de charge totales corrigées", "$P_R$", "P_R"), ("Condensation", "", "condensation"), ("Température de rosée", "$t_p$", "t_p"), ("Température moyenne des fumées", "$T_m$", "T_m"), ("Température des fumées à la sortie", "$T_{ob}$", "T_ob"), ("Température de paroi interne à la sortie", "$T_{iob}$", "T_iob"), ("Débit massique des fumées", "$\\dot{m}$", "m_dot"), ("Masse volumique moyenne des fumées", "$\\rho_m$", "rho_m"), ("Vitesse moyenne des fumées", "$w_m$", "w_m"), ("Constante des gaz des fumées", "$R$", "R"), ("Capacité calorifique spécifique", "$c_p$", "c_p"), ("Conductivité thermique", "$\\lambda_A$", "lambda_A"), ("Viscosité dynamique", "$\\eta_A$", "eta_A"),
        ("Teneur en CO2 des fumées", "$\\sigma_{CO_2}$", "sigma_CO2"), ("Teneur en H2O des fumées", "$\\sigma_{H_2O}$", "sigma_H2O"), ("Nombre de Reynolds", "$\\mathrm{Re}$", "Re"), ("Nombre de Prandtl", "$\\mathrm{Pr}$", "Pr"), ("Nombre de Nusselt", "$\\mathrm{Nu}$", "Nu"), ("Coefficient de frottement", "$\\psi$", "psi"), ("Coeff. frottement (lisse)", "$\\psi_{smooth}$", "psi_smooth"), ("Coefficient convectif interne", "$\\alpha_i$", "alpha_i"), ("Coeff. de transfert thermique (non éq.)", "$k$", "k"), ("Coeff. de transfert thermique (éq.)", "$k_b$", "k_b"), ("Coefficient de refroidissement (non éq.)", "$K$", "K"), ("Coefficient de refroidissement (éq.)", "$K_b$", "K_b"), ("Section interne", "$A$", "A"), ("Surface totale du conduit", "$A_t$", "A_t"), ("Pression atmosphérique extérieure", "$p_L$", "p_L"), ]

    # ---------- Construction du Markdown ----------
    md = []
    md.append("# Rapport de calcul — Conduit de fumée")
    md.append(f"_Selon NF EN 13384‑1:2015+A1:2019 (méthodes de calcul thermo‑aéraulique, partie 1)_")
    md.append(f"**Date** : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # --- Section 1: Données d'entrée ---
    md.append("\n## 1. Données d’entrée")
    md.append("| Paramètre | Symbole | Tirage minimale | Tirage maximale |")
    md.append("|:---|:---:|:---:|:---:|")
    for label, symbol, var in ALIASES_INPUTS:
        val_min = res_min.get(var, "N/A")
        val_max = res_max.get(var, "N/A")
        md.append(f"| {label} | {symbol} | {fmt(val_min)} | {fmt(val_max)} |")

    # --- Section 2: Formules ---
    md.append("\n## 2. Formules utilisées et références normatives")
    md.append("\n**Aéraulique**")
    md.append(r"- Tirage naturel théorique : $$ p_H = H \cdot g \cdot (\rho_L - \rho_m) \quad\text{(eq. (31), §5.10.2)} $$")
    md.append(r"- Variation de pression (vitesse) : $$ P_G = \tfrac{1}{2} (\rho_2 w_2^2 - \rho_1 w_1^2) \quad\text{(eq. (32), §5.10.3.1)} $$")
    md.append(r"- Pertes de charge totales : $$ P_R = S_E (P_{R,\text{frottement}} + P_{R,\text{raccords}}) + S_{EG} P_G \quad\text{(§5.10.1)} $$")
    md.append(r"- Frottement (Colebrook) : $$ \frac{1}{\sqrt{\psi}} = -2\,\log_{10}\!\left(\frac{2.51}{\mathrm{Re}\,\sqrt{\psi}} + \frac{r}{3.71\,D_h}\right) \quad\text{(eq. (35), §5.10.3.3)} $$")

    md.append("\n**Thermique**")
    md.append(r"- Températures (moyenne et sortie) : $$ T_m = T_u + \frac{T_e - T_u}{K}\,\bigl(1-e^{-K}\bigr) \quad\text{(eq. (16), §5.8.1)} $$ $$ T_o = T_u + (T_e - T_u)\,e^{-K} \quad\text{(eq. (17), §5.8.1)} $$")
    md.append(r"- Coefficient de refroidissement : $$ K = \frac{k\,U\,L}{\dot{m}\,c_p} \quad\text{(eq. (20), §5.8.2)} $$")
    md.append(r"- Coefficient de transfert thermique : $$ k = \left(\tfrac{1}{\alpha_i} + \tfrac{1}{\Lambda} + \tfrac{1}{\alpha_a}\right)^{-1} \text{ avec } S_H \text{ selon §5.7.7 (eq. (22), §5.8.3.1)} $$")
    md.append(r"- Transfert interne par convection : $$ \alpha_i = \frac{\mathrm{Nu}\,\lambda_A}{D_h} \quad\text{(eq. (23), §5.8.3.2)} $$")
    md.append(r"- Nombre de Nusselt : $$ \mathrm{Nu} = \left(\tfrac{\psi}{\psi_{\text{smooth}}}\right)^{0.67}\,0.0214\,\bigl(\mathrm{Re}^{0.8}-100\bigr)\,\mathrm{Pr}^{0.4}\,\left[1+\left(\tfrac{D_h}{L_{tot}}\right)^{0.67}\right] \quad\text{(eq. (24), §5.8.3.2)} $$")

    md.append("\n**Grandeurs physiques et sans dimension**")
    md.append(r"- Nombres sans dimension : $$ \mathrm{Pr} = \tfrac{c_p\,\eta_A}{\lambda_A} \quad\text{(eq. (25))} \qquad \mathrm{Re} = \tfrac{\rho_m\,w_m\,D_h}{\eta_A} \quad\text{(eq. (26), §5.8.3.2)} $$")
    md.append(r"- Grandeurs de base : $$ \rho_m = \tfrac{p_L}{R\,T_m} \quad\text{(eq. (27), §5.9.1)} \qquad w_m = \tfrac{\dot{m}}{\rho_m\,A} \quad\text{(eq. (28), §5.9.2)} $$")
    md.append(r"- Pression atmosphérique : $$ p_L = 97000\,\exp\!\left(\tfrac{-g\,z}{R_L\,T_L}\right) \quad\text{(eq. (12), §5.7.2)} $$")
    md.append("- Les propriétés des fumées ($R, c_p, \lambda_A, \eta_A, \sigma_{CO_2}, \dots$) sont issues de l'**Annexe B**.")

    # --- Section 3: Résultats ---
    md.append("\n## 3. Résultats principaux")
    md.append("| Grandeur | Symbole | Tirage minimale | Tirage maximale |")
    md.append("|:---|:---:|:---:|:---:|")
    for label, symbol, key in ALIASES_RESULTS:
        val_min = res_min.get(key, "N/A")
        val_max = res_max.get(key, "N/A")
        md.append(f"| {label} | {symbol} | {fmt(val_min)} | {fmt(val_max)} |")

    # --- Section 4: Vérifications ---
    md.append("\n## 4. Vérifications de domaine (conformité des équations)")
    re_min, re_max = res_min.get('Re', 0), res_max.get('Re', 0)
    ok_re_min, ok_re_max = (2300 <= re_min <= 1e7), (2300 <= re_max <= 1e7)
    md.append(f"- **Domaine de validité pour Nu (Re)** : `2300 ≤ Re ≤ 1e7`")
    md.append(f"  - Tirage min: {'OK' if ok_re_min else 'HORS DOMAINE'} (Re = {fmt(re_min)})")
    md.append(f"  - Tirage max: {'OK' if ok_re_max else 'HORS DOMAINE'} (Re = {fmt(re_max)})")
    pr_min, pr_max = res_min.get('Pr', 0), res_max.get('Pr', 0)
    ok_pr_min, ok_pr_max = (0.6 < pr_min < 1.5), (0.6 < pr_max < 1.5)
    md.append(f"- **Domaine de validité pour Nu (Pr)** : `0.6 < Pr < 1.5`")
    md.append(f"  - Tirage min: {'OK' if ok_pr_min else 'HORS DOMAINE'} (Pr = {fmt(pr_min)})")
    md.append(f"  - Tirage max: {'OK' if ok_pr_max else 'HORS DOMAINE'} (Pr = {fmt(pr_max)})")

    # --- Section 5: Remarques ---
    md.append("\n## 5. Remarques et cas particuliers")
    md.append("- **αa (extérieur/intérieur)** : la norme fixe 23 W/(m²·K) à l’extérieur et 8 W/(m²·K) à l’intérieur (§5.8.3.3). Une interpolation linéaire est utilisée.")
    md.append("- **T_u = 15 °C** car environnement chauffé avec surface de cheminée extérieure < 1/4 de la surface totale (§ 5.7.1.3 NOTE 1).")
    md.append("- Les résultats sont obtenus par approximation itérative (valeur initiale de T_m = 223 °C).")
    md.append("---")

    # ---------- Écriture du fichier ----------
    output_dir = os.path.dirname(filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    content = "\n".join(md) + "\n"
    with open(filename, "w", encoding="utf-8") as f:
        f.write(content)

    print(f"Rapport Markdown généré : {filename}")
    return filename


def do_report_english(res_min: dict, res_max: dict, filename: str = "/note_calculs/note_de_calcul_english.md") -> str:
    """
    Generates a comparative Markdown+LaTeX report for the duct calculation and saves it.
    The report includes normative references (EN 13384-1:2015+A1:2019) for
    minimum and maximum draft conditions.
    Returns the path to the written file.
    """
    import datetime
    import os
    import math

    # ---------- Helper functions ----------
    def _is_qty(x):
        return hasattr(x, "to") and hasattr(x, "magnitude")

    def fmt(x, unit=None, digits=3):
        """
        Formats pint quantities or plain numbers using SI units
        and up to `digits` significant figures. Uses scientific notation
        if |value| > 1e4 or 0 < |value| < 1e-3.
        """

        def _is_qty(v):
            return hasattr(v, "to") and hasattr(v, "magnitude") and hasattr(v, "units")

        def _format_sig(v: float, sig: int = 3) -> str:
            if not math.isfinite(v): return str(v)
            if v == 0: return "0"
            av = abs(v)
            use_sci = (av > 1e4) or (av > 0 and av < 1e-3)
            if use_sci:
                return f"{v:.{max(sig - 1, 0)}e}"
            exp10 = math.floor(math.log10(av)) if av > 0 else 0
            decimals = sig - 1 - exp10
            rounded = round(v, int(decimals))
            s = f"{rounded:.{max(int(decimals), 0)}f}"
            if "." in s: s = s.rstrip("0").rstrip(".")
            return s

        def _to_si(q):
            if unit is not None:
                try:
                    return q.to(unit)
                except Exception:
                    pass
            si_candidates = ["1", "K", "m", "m**2", "m**3", "kg", "s", "kg/s", "kg/m**3", "m/s", "m/s**2", "Pa", "W", "J", "N", "W/(m*K)", "W/(m**2*K)", "Pa*s", "J/(kg*K)", "m*K/W", "m**2*K/W", "degC"]
            for ustr in si_candidates:
                try:
                    return q.to(ustr)
                except Exception:
                    continue
            try:
                return q.to_base_units()
            except Exception:
                return q

        try:
            if _is_qty(x):
                q = _to_si(x)
                val = float(q.magnitude)
                s = _format_sig(val, digits)
                try:
                    u_str = "{:~}".format(q.units)
                except Exception:
                    u_str = str(q.units)
                if u_str in ("dimensionless", ""): return s
                return f"{s} {u_str}"

            elif isinstance(x, bool):
                return "Yes" if x else "No"
            elif isinstance(x, (int, float)):
                return _format_sig(float(x), digits)
            else:
                return str(x)
        except Exception:
            return str(x)

    # --- Aliases for tables (label, symbol, dict_key) ---
    ALIASES_INPUTS = [("Outside Air Temperature", "$T_L$", "T_L"), ("Flue Gas Inlet Temperature", "$T_e$", "T_e"), ("Ambient Air Temperature (Duct)", "$T_u$", "T_u"), ("Ambient Air Temperature (Outlet)", "$T_{uo}$", "T_uo"), ("Internal Hydraulic Diameter", "$D_h$", "D_h"), ("Flue Duct Length", "$L$", "L"), ("Thermal Resistance of the Duct", "$1/\\Lambda$", "Res_therm"), ("Nominal Heat Output", "$Q_N$", "Q_N"), ("Efficiency", "$\\eta$", "eta"), ("Wind Pressure", "$P_L$ (wind)", "P_L_vent"), ]
    ALIASES_RESULTS = [("Theoretical Natural Draft", "$P_H$", "P_H"), ("Pressure Loss (Friction + Fittings)", "$P_E$", "P_E"), ("Pressure Variation (Velocity)", "$P_G$", "P_G"), ("Total Corrected Pressure Loss", "$P_R$", "P_R"), ("Condensation", "", "condensation"), ("Dew Point Temperature", "$t_p$", "t_p"), ("Mean Flue Gas Temperature", "$T_m$", "T_m"), ("Flue Gas Outlet Temperature", "$T_{ob}$", "T_ob"), ("Inner Wall Temperature at Outlet", "$T_{iob}$", "T_iob"), ("Flue Gas Mass Flow Rate", "$\\dot{m}$", "m_dot"), ("Mean Flue Gas Density", "$\\rho_m$", "rho_m"), ("Mean Flue Gas Velocity", "$w_m$", "w_m"), ("Flue Gas Constant", "$R$", "R"), ("Specific Heat Capacity", "$c_p$", "c_p"), ("Thermal Conductivity", "$\\lambda_A$", "lambda_A"), ("Dynamic Viscosity", "$\\eta_A$", "eta_A"),
        ("Flue Gas CO2 Content", "$\\sigma_{CO_2}$", "sigma_CO2"), ("Flue Gas H2O Content", "$\\sigma_{H_2O}$", "sigma_H2O"), ("Reynolds Number", "$\\mathrm{Re}$", "Re"), ("Prandtl Number", "$\\mathrm{Pr}$", "Pr"), ("Nusselt Number", "$\\mathrm{Nu}$", "Nu"), ("Friction Factor", "$\\psi$", "psi"), ("Friction Factor (Smooth)", "$\\psi_{smooth}$", "psi_smooth"), ("Internal Convective Coefficient", "$\\alpha_i$", "alpha_i"), ("Heat Transfer Coeff. (Unbalanced)", "$k$", "k"), ("Heat Transfer Coeff. (Balanced)", "$k_b$", "k_b"), ("Cooling Coefficient (Unbalanced)", "$K$", "K"), ("Cooling Coefficient (Balanced)", "$K_b$", "K_b"), ("Internal Cross-Section", "$A$", "A"), ("Total Duct Surface Area", "$A_t$", "A_t"), ("Outside Atmospheric Pressure", "$p_L$", "p_L"), ]

    # ---------- Building the Markdown ----------
    md = []
    md.append("# Calculation Report — Flue Duct")
    md.append(f"_According to NF EN 13384‑1:2015+A1:2019 (Thermal and fluid dynamic calculation methods, Part 1)_")
    md.append(f"**Date**: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # --- Section 1: Input Data ---
    md.append("\n## 1. Input Data")
    md.append("| Parameter | Symbol | Minimum Draft | Maximum Draft |")
    md.append("|:---|:---:|:---:|:---:|")
    for label, symbol, var in ALIASES_INPUTS:
        val_min = res_min.get(var, "N/A")
        val_max = res_max.get(var, "N/A")
        md.append(f"| {label} | {symbol} | {fmt(val_min)} | {fmt(val_max)} |")

    # --- Section 2: Formulas Used and Standard References ---
    md.append("\n## 2. Formulas Used and Standard References")
    md.append("\n**Aerodynamics**")
    md.append(r"- Theoretical natural draft: $$ p_H = H \cdot g \cdot (\rho_L - \rho_m) \quad\text{(eq. (31), §5.10.2)} $$")
    md.append(r"- Pressure variation (velocity): $$ P_G = \tfrac{1}{2} (\rho_2 w_2^2 - \rho_1 w_1^2) \quad\text{(eq. (32), §5.10.3.1)} $$")
    md.append(r"- Total pressure loss: $$ P_R = S_E (P_{R,\text{friction}} + P_{R,\text{fittings}}) + S_{EG} P_G \quad\text{(§5.10.1)} $$")
    md.append(r"- Friction (Colebrook): $$ \frac{1}{\sqrt{\psi}} = -2\,\log_{10}\!\left(\frac{2.51}{\mathrm{Re}\,\sqrt{\psi}} + \frac{r}{3.71\,D_h}\right) \quad\text{(eq. (35), §5.10.3.3)} $$")

    md.append("\n**Thermodynamics**")
    md.append(r"- Temperatures (mean and outlet): $$ T_m = T_u + \frac{T_e - T_u}{K}\,\bigl(1-e^{-K}\bigr) \quad\text{(eq. (16), §5.8.1)} $$ $$ T_o = T_u + (T_e - T_u)\,e^{-K} \quad\text{(eq. (17), §5.8.1)} $$")
    md.append(r"- Cooling coefficient: $$ K = \frac{k\,U\,L}{\dot{m}\,c_p} \quad\text{(eq. (20), §5.8.2)} $$")
    md.append(r"- Heat transfer coefficient: $$ k = \left(\tfrac{1}{\alpha_i} + \tfrac{1}{\Lambda} + \tfrac{1}{\alpha_a}\right)^{-1} \text{ with } S_H \text{ from §5.7.7 (eq. (22), §5.8.3.1)} $$")
    md.append(r"- Internal convection heat transfer: $$ \alpha_i = \frac{\mathrm{Nu}\,\lambda_A}{D_h} \quad\text{(eq. (23), §5.8.3.2)} $$")
    md.append(r"- Nusselt number: $$ \mathrm{Nu} = \left(\tfrac{\psi}{\psi_{\text{smooth}}}\right)^{0.67}\,0.0214\,\bigl(\mathrm{Re}^{0.8}-100\bigr)\,\mathrm{Pr}^{0.4}\,\left[1+\left(\tfrac{D_h}{L_{tot}}\right)^{0.67}\right] \quad\text{(eq. (24), §5.8.3.2)} $$")

    md.append("\n**Physical Quantities and Dimensionless Numbers**")
    md.append(r"- Dimensionless numbers: $$ \mathrm{Pr} = \tfrac{c_p\,\eta_A}{\lambda_A} \quad\text{(eq. (25))} \qquad \mathrm{Re} = \tfrac{\rho_m\,w_m\,D_h}{\eta_A} \quad\text{(eq. (26), §5.8.3.2)} $$")
    md.append(r"- Basic quantities: $$ \rho_m = \tfrac{p_L}{R\,T_m} \quad\text{(eq. (27), §5.9.1)} \qquad w_m = \tfrac{\dot{m}}{\rho_m\,A} \quad\text{(eq. (28), §5.9.2)} $$")
    md.append(r"- Atmospheric pressure: $$ p_L = 97000\,\exp\!\left(\tfrac{-g\,z}{R_L\,T_L}\right) \quad\text{(eq. (12), §5.7.2)} $$")
    md.append("- Flue gas properties ($R, c_p, \lambda_A, \eta_A, \sigma_{CO_2}, \dots$) are taken from **Annex B**.")

    # --- Section 3: Main Results ---
    md.append("\n## 3. Main Results")
    md.append("| Quantity | Symbol | Minimum Draft | Maximum Draft |")
    md.append("|:---|:---:|:---:|:---:|")
    for label, symbol, key in ALIASES_RESULTS:
        val_min = res_min.get(key, "N/A")
        val_max = res_max.get(key, "N/A")
        md.append(f"| {label} | {symbol} | {fmt(val_min)} | {fmt(val_max)} |")

    # --- Section 4: Domain Checks (Equation Conformity) ---
    md.append("\n## 4. Domain Checks (Equation Conformity)")
    re_min, re_max = res_min.get('Re', 0), res_max.get('Re', 0)
    ok_re_min, ok_re_max = (2300 <= re_min <= 1e7), (2300 <= re_max <= 1e7)
    md.append(f"- **Validity Domain for Nu (Re)**: `2300 ≤ Re ≤ 1e7`")
    md.append(f"  - Min Draft: {'OK' if ok_re_min else 'OUT OF DOMAIN'} (Re = {fmt(re_min)})")
    md.append(f"  - Max Draft: {'OK' if ok_re_max else 'OUT OF DOMAIN'} (Re = {fmt(re_max)})")
    pr_min, pr_max = res_min.get('Pr', 0), res_max.get('Pr', 0)
    ok_pr_min, ok_pr_max = (0.6 < pr_min < 1.5), (0.6 < pr_max < 1.5)
    md.append(f"- **Validity Domain for Nu (Pr)**: `0.6 < Pr < 1.5`")
    md.append(f"  - Min Draft: {'OK' if ok_pr_min else 'OUT OF DOMAIN'} (Pr = {fmt(pr_min)})")
    md.append(f"  - Max Draft: {'OK' if ok_pr_max else 'OUT OF DOMAIN'} (Pr = {fmt(pr_max)})")

    # --- Section 5: Notes and Special Cases ---
    md.append("\n## 5. Notes and Special Cases")
    md.append("- **αa (external/internal)**: the standard specifies 23 W/(m²·K) for the exterior and 8 W/(m²·K) for the interior (§5.8.3.3). A linear interpolation is used.")
    md.append("- **T_u = 15 °C** due to a heated environment with an exterior chimney surface area < 1/4 of the total surface area (§ 5.7.1.3 NOTE 1).")
    md.append("- The results are obtained by iterative approximation (initial T_m value of 223 °C).")
    md.append("---")

    # ---------- Writing the file ----------
    output_dir = os.path.dirname(filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    content = "\n".join(md) + "\n"
    with open(filename, "w", encoding="utf-8") as f:
        f.write(content)

    print(f"Markdown report generated: {filename}")
    return filename