# Calculation Report — Flue Duct
_According to NF EN 13384‑1:2015+A1:2019 (Thermal and fluid dynamic calculation methods, Part 1)_
**Date**: 2025-10-06 16:47:32

## 1. Input Data
| Parameter | Symbol | Minimum Draft | Maximum Draft |
|:---|:---:|:---:|:---:|
| Outside Air Temperature | $T_L$ | 273 K | 288 K |
| Flue Gas Inlet Temperature | $T_e$ | 583 K | 583 K |
| Ambient Air Temperature (Duct) | $T_u$ | 288 K | 288 K |
| Ambient Air Temperature (Outlet) | $T_{uo}$ | 273 K | 288 K |
| Internal Hydraulic Diameter | $D_h$ | 0.2 m | 0.2 m |
| Flue Duct Length | $L$ | 7.7 m | 7.7 m |
| Thermal Resistance of the Duct | $1/\Lambda$ | 0 K * m ** 2 / W | 0 K * m ** 2 / W |
| Nominal Heat Output | $Q_N$ | 1.40e+05 W | 1.40e+05 W |
| Efficiency | $\eta$ | 0.86 | 0.86 |
| Wind Pressure | $P_L$ (wind) | N/A | N/A |

## 2. Formulas Used and Standard References

**Aerodynamics**
- Theoretical natural draft: $$ p_H = H \cdot g \cdot (\rho_L - \rho_m) \quad\text{(eq. (31), §5.10.2)} $$
- Pressure variation (velocity): $$ P_G = \tfrac{1}{2} (\rho_2 w_2^2 - \rho_1 w_1^2) \quad\text{(eq. (32), §5.10.3.1)} $$
- Total pressure loss: $$ P_R = S_E (P_{R,\text{friction}} + P_{R,\text{fittings}}) + S_{EG} P_G \quad\text{(§5.10.1)} $$
- Friction (Colebrook): $$ \frac{1}{\sqrt{\psi}} = -2\,\log_{10}\!\left(\frac{2.51}{\mathrm{Re}\,\sqrt{\psi}} + \frac{r}{3.71\,D_h}\right) \quad\text{(eq. (35), §5.10.3.3)} $$

**Thermodynamics**
- Temperatures (mean and outlet): $$ T_m = T_u + \frac{T_e - T_u}{K}\,\bigl(1-e^{-K}\bigr) \quad\text{(eq. (16), §5.8.1)} $$ $$ T_o = T_u + (T_e - T_u)\,e^{-K} \quad\text{(eq. (17), §5.8.1)} $$
- Cooling coefficient: $$ K = \frac{k\,U\,L}{\dot{m}\,c_p} \quad\text{(eq. (20), §5.8.2)} $$
- Heat transfer coefficient: $$ k = \left(\tfrac{1}{\alpha_i} + \tfrac{1}{\Lambda} + \tfrac{1}{\alpha_a}\right)^{-1} \text{ with } S_H \text{ from §5.7.7 (eq. (22), §5.8.3.1)} $$
- Internal convection heat transfer: $$ \alpha_i = \frac{\mathrm{Nu}\,\lambda_A}{D_h} \quad\text{(eq. (23), §5.8.3.2)} $$
- Nusselt number: $$ \mathrm{Nu} = \left(\tfrac{\psi}{\psi_{\text{smooth}}}\right)^{0.67}\,0.0214\,\bigl(\mathrm{Re}^{0.8}-100\bigr)\,\mathrm{Pr}^{0.4}\,\left[1+\left(\tfrac{D_h}{L_{tot}}\right)^{0.67}\right] \quad\text{(eq. (24), §5.8.3.2)} $$

**Physical Quantities and Dimensionless Numbers**
- Dimensionless numbers: $$ \mathrm{Pr} = \tfrac{c_p\,\eta_A}{\lambda_A} \quad\text{(eq. (25))} \qquad \mathrm{Re} = \tfrac{\rho_m\,w_m\,D_h}{\eta_A} \quad\text{(eq. (26), §5.8.3.2)} $$
- Basic quantities: $$ \rho_m = \tfrac{p_L}{R\,T_m} \quad\text{(eq. (27), §5.9.1)} \qquad w_m = \tfrac{\dot{m}}{\rho_m\,A} \quad\text{(eq. (28), §5.9.2)} $$
- Atmospheric pressure: $$ p_L = 97000\,\exp\!\left(\tfrac{-g\,z}{R_L\,T_L}\right) \quad\text{(eq. (12), §5.7.2)} $$
- Flue gas properties ($R, c_p, \lambda_A, \eta_A, \sigma_{CO_2}, \dots$) are taken from **Annex B**.

## 3. Main Results
| Quantity | Symbol | Minimum Draft | Maximum Draft |
|:---|:---:|:---:|:---:|
| Theoretical Natural Draft | $P_H$ | 51.2 Pa | 40.2 Pa |
| Pressure Loss (Friction + Fittings) | $P_E$ | 10 Pa | 9.77 Pa |
| Pressure Variation (Velocity) | $P_G$ | -2.15 Pa | -2.15 Pa |
| Total Corrected Pressure Loss | $P_R$ | 7.88 Pa | 9.57 Pa |
| Condensation |  | No | No |
| Dew Point Temperature | $t_p$ | 328 K | 328 K |
| Mean Flue Gas Temperature | $T_m$ | 539 K | 525 K |
| Flue Gas Outlet Temperature | $T_{ob}$ | 499 K | 499 K |
| Inner Wall Temperature at Outlet | $T_{iob}$ | 343 K | 353 K |
| Flue Gas Mass Flow Rate | $\dot{m}$ | 0.0683 kg / s | 0.0683 kg / s |
| Mean Flue Gas Density | $\rho_m$ | 0.602 kg / m ** 3 | 0.617 kg / m ** 3 |
| Mean Flue Gas Velocity | $w_m$ | 3.61 m / s | 3.52 m / s |
| Flue Gas Constant | $R$ | 298 J / K / kg | 298 J / K / kg |
| Specific Heat Capacity | $c_p$ | 1150 J / K / kg | 1150 J / K / kg |
| Thermal Conductivity | $\lambda_A$ | 0.0396 W / K / m | 0.0387 W / K / m |
| Dynamic Viscosity | $\eta_A$ | 2.61e-05 Pa * s | 2.56e-05 Pa * s |
| Flue Gas CO2 Content | $\sigma_{CO_2}$ | 10.2 | 10.2 |
| Flue Gas H2O Content | $\sigma_{H_2O}$ | 16.5 | 16.5 |
| Reynolds Number | $\mathrm{Re}$ | 1.67e+04 | 1.70e+04 |
| Prandtl Number | $\mathrm{Pr}$ | 0.757 | 0.757 |
| Nusselt Number | $\mathrm{Nu}$ | 56.7 | 57.7 |
| Friction Factor | $\psi$ | 0.0351 | 0.0351 |
| Friction Factor (Smooth) | $\psi_{smooth}$ | 0.0271 | 0.0269 |
| Internal Convective Coefficient | $\alpha_i$ | 11.2 W / K / m ** 2 | 11.2 W / K / m ** 2 |
| Heat Transfer Coeff. (Unbalanced) | $k$ | 5.46 W / K / m ** 2 | 7.32 W / K / m ** 2 |
| Heat Transfer Coeff. (Balanced) | $k_b$ | 5.46 W / K / m ** 2 | 5.45 W / K / m ** 2 |
| Cooling Coefficient (Unbalanced) | $K$ | 0.337 | 0.453 |
| Cooling Coefficient (Balanced) | $K_b$ | 0.337 | 0.337 |
| Internal Cross-Section | $A$ | 0.0314 m ** 2 | 0.0314 m ** 2 |
| Total Duct Surface Area | $A_t$ | 4.84 m ** 2 | 4.84 m ** 2 |
| Outside Atmospheric Pressure | $p_L$ | 9.65e+04 Pa | 9.65e+04 Pa |

## 4. Domain Checks (Equation Conformity)
- **Validity Domain for Nu (Re)**: `2300 ≤ Re ≤ 1e7`
  - Min Draft: OK (Re = 1.67e+04)
  - Max Draft: OK (Re = 1.70e+04)
- **Validity Domain for Nu (Pr)**: `0.6 < Pr < 1.5`
  - Min Draft: OK (Pr = 0.757)
  - Max Draft: OK (Pr = 0.757)

## 5. Notes and Special Cases
- **αa (external/internal)**: the standard specifies 23 W/(m²·K) for the exterior and 8 W/(m²·K) for the interior (§5.8.3.3). A linear interpolation is used.
- **T_u = 15 °C** due to a heated environment with an exterior chimney surface area < 1/4 of the total surface area (§ 5.7.1.3 NOTE 1).
- The results are obtained by iterative approximation (initial T_m value of 223 °C).
---
