# Rapport de calcul — Conduit de fumée
_Selon NF EN 13384‑1:2015+A1:2019 (méthodes de calcul thermo‑aéraulique, partie 1)_
**Date** : 2025-10-06 16:47:32

## 1. Données d’entrée
| Paramètre | Symbole | Tirage minimale | Tirage maximale |
|:---|:---:|:---:|:---:|
| Température de l’air extérieur | $T_L$ | 273 K | 288 K |
| Température des fumées à l’admission | $T_e$ | 583 K | 583 K |
| Température de l’air ambiant (conduit) | $T_u$ | 288 K | 288 K |
| Température de l’air ambiant (sortie) | $T_{uo}$ | 273 K | 288 K |
| Diamètre hydraulique intérieur | $D_h$ | 0.2 m | 0.2 m |
| Longueur du conduit de fumée | $L$ | 7.7 m | 7.7 m |
| Résistance thermique du conduit | $1/\Lambda$ | 0 K * m ** 2 / W | 0 K * m ** 2 / W |
| Puissance utile nominale | $Q_N$ | 1.40e+05 W | 1.40e+05 W |
| Rendement | $\eta$ | 0.86 | 0.86 |
| Pression du vent | $P_L$ (vent) | N/A | N/A |

## 2. Formules utilisées et références normatives

**Aéraulique**
- Tirage naturel théorique : $$ p_H = H \cdot g \cdot (\rho_L - \rho_m) \quad\text{(eq. (31), §5.10.2)} $$
- Variation de pression (vitesse) : $$ P_G = \tfrac{1}{2} (\rho_2 w_2^2 - \rho_1 w_1^2) \quad\text{(eq. (32), §5.10.3.1)} $$
- Pertes de charge totales : $$ P_R = S_E (P_{R,\text{frottement}} + P_{R,\text{raccords}}) + S_{EG} P_G \quad\text{(§5.10.1)} $$
- Frottement (Colebrook) : $$ \frac{1}{\sqrt{\psi}} = -2\,\log_{10}\!\left(\frac{2.51}{\mathrm{Re}\,\sqrt{\psi}} + \frac{r}{3.71\,D_h}\right) \quad\text{(eq. (35), §5.10.3.3)} $$

**Thermique**
- Températures (moyenne et sortie) : $$ T_m = T_u + \frac{T_e - T_u}{K}\,\bigl(1-e^{-K}\bigr) \quad\text{(eq. (16), §5.8.1)} $$ $$ T_o = T_u + (T_e - T_u)\,e^{-K} \quad\text{(eq. (17), §5.8.1)} $$
- Coefficient de refroidissement : $$ K = \frac{k\,U\,L}{\dot{m}\,c_p} \quad\text{(eq. (20), §5.8.2)} $$
- Coefficient de transfert thermique : $$ k = \left(\tfrac{1}{\alpha_i} + \tfrac{1}{\Lambda} + \tfrac{1}{\alpha_a}\right)^{-1} \text{ avec } S_H \text{ selon §5.7.7 (eq. (22), §5.8.3.1)} $$
- Transfert interne par convection : $$ \alpha_i = \frac{\mathrm{Nu}\,\lambda_A}{D_h} \quad\text{(eq. (23), §5.8.3.2)} $$
- Nombre de Nusselt : $$ \mathrm{Nu} = \left(\tfrac{\psi}{\psi_{\text{smooth}}}\right)^{0.67}\,0.0214\,\bigl(\mathrm{Re}^{0.8}-100\bigr)\,\mathrm{Pr}^{0.4}\,\left[1+\left(\tfrac{D_h}{L_{tot}}\right)^{0.67}\right] \quad\text{(eq. (24), §5.8.3.2)} $$

**Grandeurs physiques et sans dimension**
- Nombres sans dimension : $$ \mathrm{Pr} = \tfrac{c_p\,\eta_A}{\lambda_A} \quad\text{(eq. (25))} \qquad \mathrm{Re} = \tfrac{\rho_m\,w_m\,D_h}{\eta_A} \quad\text{(eq. (26), §5.8.3.2)} $$
- Grandeurs de base : $$ \rho_m = \tfrac{p_L}{R\,T_m} \quad\text{(eq. (27), §5.9.1)} \qquad w_m = \tfrac{\dot{m}}{\rho_m\,A} \quad\text{(eq. (28), §5.9.2)} $$
- Pression atmosphérique : $$ p_L = 97000\,\exp\!\left(\tfrac{-g\,z}{R_L\,T_L}\right) \quad\text{(eq. (12), §5.7.2)} $$
- Les propriétés des fumées ($R, c_p, \lambda_A, \eta_A, \sigma_{CO_2}, \dots$) sont issues de l'**Annexe B**.

## 3. Résultats principaux
| Grandeur | Symbole | Tirage minimale | Tirage maximale |
|:---|:---:|:---:|:---:|
| Tirage naturel théorique | $P_H$ | 51.2 Pa | 40.2 Pa |
| Pertes de charge (frottement + raccords) | $P_E$ | 10 Pa | 9.77 Pa |
| Variation de pression (vitesse) | $P_G$ | -2.15 Pa | -2.15 Pa |
| Pertes de charge totales corrigées | $P_R$ | 7.88 Pa | 9.57 Pa |
| Condensation |  | Non | Non |
| Température de rosée | $t_p$ | 328 K | 328 K |
| Température moyenne des fumées | $T_m$ | 539 K | 525 K |
| Température des fumées à la sortie | $T_{ob}$ | 499 K | 499 K |
| Température de paroi interne à la sortie | $T_{iob}$ | 343 K | 353 K |
| Débit massique des fumées | $\dot{m}$ | 0.0683 kg / s | 0.0683 kg / s |
| Masse volumique moyenne des fumées | $\rho_m$ | 0.602 kg / m ** 3 | 0.617 kg / m ** 3 |
| Vitesse moyenne des fumées | $w_m$ | 3.61 m / s | 3.52 m / s |
| Constante des gaz des fumées | $R$ | 298 J / K / kg | 298 J / K / kg |
| Capacité calorifique spécifique | $c_p$ | 1150 J / K / kg | 1150 J / K / kg |
| Conductivité thermique | $\lambda_A$ | 0.0396 W / K / m | 0.0387 W / K / m |
| Viscosité dynamique | $\eta_A$ | 2.61e-05 Pa * s | 2.56e-05 Pa * s |
| Teneur en CO2 des fumées | $\sigma_{CO_2}$ | 10.2 | 10.2 |
| Teneur en H2O des fumées | $\sigma_{H_2O}$ | 16.5 | 16.5 |
| Nombre de Reynolds | $\mathrm{Re}$ | 1.67e+04 | 1.70e+04 |
| Nombre de Prandtl | $\mathrm{Pr}$ | 0.757 | 0.757 |
| Nombre de Nusselt | $\mathrm{Nu}$ | 56.7 | 57.7 |
| Coefficient de frottement | $\psi$ | 0.0351 | 0.0351 |
| Coeff. frottement (lisse) | $\psi_{smooth}$ | 0.0271 | 0.0269 |
| Coefficient convectif interne | $\alpha_i$ | 11.2 W / K / m ** 2 | 11.2 W / K / m ** 2 |
| Coeff. de transfert thermique (non éq.) | $k$ | 5.46 W / K / m ** 2 | 7.32 W / K / m ** 2 |
| Coeff. de transfert thermique (éq.) | $k_b$ | 5.46 W / K / m ** 2 | 5.45 W / K / m ** 2 |
| Coefficient de refroidissement (non éq.) | $K$ | 0.337 | 0.453 |
| Coefficient de refroidissement (éq.) | $K_b$ | 0.337 | 0.337 |
| Section interne | $A$ | 0.0314 m ** 2 | 0.0314 m ** 2 |
| Surface totale du conduit | $A_t$ | 4.84 m ** 2 | 4.84 m ** 2 |
| Pression atmosphérique extérieure | $p_L$ | 9.65e+04 Pa | 9.65e+04 Pa |

## 4. Vérifications de domaine (conformité des équations)
- **Domaine de validité pour Nu (Re)** : `2300 ≤ Re ≤ 1e7`
  - Tirage min: OK (Re = 1.67e+04)
  - Tirage max: OK (Re = 1.70e+04)
- **Domaine de validité pour Nu (Pr)** : `0.6 < Pr < 1.5`
  - Tirage min: OK (Pr = 0.757)
  - Tirage max: OK (Pr = 0.757)

## 5. Remarques et cas particuliers
- **αa (extérieur/intérieur)** : la norme fixe 23 W/(m²·K) à l’extérieur et 8 W/(m²·K) à l’intérieur (§5.8.3.3). Une interpolation linéaire est utilisée.
- **T_u = 15 °C** car environnement chauffé avec surface de cheminée extérieure < 1/4 de la surface totale (§ 5.7.1.3 NOTE 1).
- Les résultats sont obtenus par approximation itérative (valeur initiale de T_m = 223 °C).
---
