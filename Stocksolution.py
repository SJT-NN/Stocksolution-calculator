#!/usr/bin/env python3
import math
import streamlit as st
from molmass import Formula

# --- Molar mass from multiple formulas ---
def molar_mass_multi(formulas_str):
    total_mass = 0.0
    for f_str in formulas_str.split(','):
        f_str = f_str.strip()
        if f_str:
            f = Formula(f_str)
            total_mass += f.mass
    return total_mass

# --- Mass required given purity ---
def mass_required_molar(M, conc_mol_L, vol_L, purity):
    pure_mass = M * conc_mol_L * vol_L
    return pure_mass / purity

# --- Convert mg/L to mol/L ---
def mgL_to_molL(conc_mg_L, M):
    return conc_mg_L / 1000.0 / M

# --- Convert mol/L to mg/L ---
def molL_to_mgL(conc_mol_L, M):
    return conc_mol_L * M * 1000.0

# --- Uncertainty of concentration ---
def conc_uncertainty(M, m_weighed, u_mass, vol_L, u_vol_L, purity, u_purity):
    dc_dm = purity / (M * vol_L)
    dc_dV = -(m_weighed * purity) / (M * vol_L**2)
    dc_dp = m_weighed / (M * vol_L)
    return math.sqrt(
        (dc_dm * u_mass) ** 2 +
        (dc_dV * u_vol_L) ** 2 +
        (dc_dp * u_purity) ** 2
    )

# --- Safe float parser ---
def parse_float(text, default=0.0):
    try:
        return float(text)
    except ValueError:
        return default

# --- UI setup ---
st.set_page_config(page_title="Solution Prep Calculator", page_icon="🧪")
st.title("🧪 Multi‑Component Solution Preparation Calculator")

# Chemical formulas
formulas_str = st.text_input(
    "Chemical formula(s) (comma‑separated)",
    "NaCl"
)

# Target concentration: choose unit
conc_col, unit_col = st.columns([2, 1])
with conc_col:
    target_conc = parse_float(st.text_input("Target concentration", "0.1"))
with unit_col:
    conc_unit = st.selectbox("Concentration unit", ["mol/L", "mg/L"])

# Volume input (unless using mass/density)
use_mass_density = st.checkbox("Enter solution mass & density instead of volume", value=False)

if not use_mass_density:
    vol_col, unit_vol_col = st.columns([2, 1])
    with vol_col:
        vol_value = parse_float(st.text_input("Solution volume", "100.0"))
    with unit_vol_col:
        vol_unit = st.selectbox("Volume unit", ["mL", "L"], index=0)
    vol_L = vol_value / 1000.0 if vol_unit == "mL" else vol_value
    u_vol_L = parse_float(st.text_input("Volume uncertainty", "0.1")) / (1000.0 if vol_unit == "mL" else 1)
else:
    sol_mass_g = parse_float(st.text_input("Total solution mass [g]", "100.0"))
    sol_density_gmL = parse_float(st.text_input("Solution density [g/mL]", "1.00"))
    vol_L = sol_mass_g / (sol_density_gmL * 1000.0)
    u_vol_L = 0.0  # could add density & mass uncertainties later

# Purity & uncertainties
col3, col4, col5 = st.columns(3)
with col3:
    purity_percent = parse_float(st.text_input("Purity [%]", "99.5"))
with col4:
    u_purity_percent = parse_float(st.text_input("Purity uncertainty [%]", "0.05"))
with col5:
    u_scale = parse_float(st.text_input("Scale uncertainty [g]", "0.001"))

# --- Convert purity to fraction ---
purity = purity_percent / 100.0
u_purity = u_purity_percent / 100.0

# --- Computations ---
if formulas_str and target_conc > 0 and vol_L > 0 and purity > 0:
    M_total = molar_mass_multi(formulas_str)
    # Convert all to mol/L internally
    conc_mol_L = target_conc if conc_unit == "mol/L" else mgL_to_molL(target_conc, M_total)
    m_req = mass_required_molar(M_total, conc_mol_L, vol_L, purity)
    u_c_solution = conc_uncertainty(M_total, m_req, u_scale, vol_L, u_vol_L, purity, u_purity)
    rel_u_c_solution = (u_c_solution / conc_mol_L) * 100

    # --- Output ---
    st.subheader("Results")
    st.write(f"**Total molar mass:** {M_total:.5f} g/mol")
    st.write(f"**Required total mass to weigh:** {m_req:.5f} g")
    st.write(f"**Target concentration (mol/L):** {conc_mol_L:.6f}")
    st.write(f"**Target concentration (mg/L):** {molL_to_mgL(conc_mol_L, M_total):.3f}")
    st.write(f"**Uncertainty of resulting concentration:** ± {u_c_solution:.6f} mol/L")
    st.write(f"**Relative concentration uncertainty:** {rel_u_c_solution:.3f} %")
