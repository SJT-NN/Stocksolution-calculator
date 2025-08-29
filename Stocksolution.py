#!/usr/bin/env python3
import math
import streamlit as st
from molmass import Formula

# --- Molar mass from multiple formulas ---
def molar_mass_multi(formulas_str):
    """Takes comma-separated formulas and returns total molar mass."""
    total_mass = 0.0
    for f_str in formulas_str.split(','):
        f_str = f_str.strip()
        if f_str:
            f = Formula(f_str)
            total_mass += f.mass
    return total_mass

# --- Required mass, purity‑corrected ---
def mass_required(M, conc, vol_L, purity):
    pure_mass = M * conc * vol_L
    return pure_mass / purity  # weigh more if purity < 100%

# --- Uncertainty of concentration ---
def conc_uncertainty(M, m_weighed, u_mass,
                     vol_L, u_vol_L,
                     purity, u_purity):
    """
    Standard uncertainty in concentration [mol/L] from:
    weighing (u_mass), volume delivery (u_vol_L),
    and purity uncertainty (fractional).
    """
    # c = (m * p) / (M * V)
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

# --- UI ---
st.set_page_config(page_title="Solution Prep Calculator", page_icon="🧪")
st.title("🧪 Multi‑Component Solution Preparation Calculator")

# Multiple formulas
formulas_str = st.text_input(
    "Chemical formula(s) (comma‑separated, e.g. 'NaCl, H2O')",
    "NaCl"
)

# Target concentration
conc = parse_float(st.text_input("Target concentration [mol/L]", "0.1"))

# Volume and units - same line
vol_col, unit_col = st.columns([2, 1])
with vol_col:
    vol_value = parse_float(st.text_input("Solution volume", "100.0"))
with unit_col:
    vol_unit = st.selectbox("Unit", ["mL", "L"], index=0)

# Volume uncertainty and its units - same line
uvol_col, uunit_col = st.columns([2, 1])
with uvol_col:
    u_vol_value = parse_float(st.text_input("Volume uncertainty", "0.1"))
with uunit_col:
    u_vol_unit = st.selectbox("Unit ", ["mL", "L"], index=0)

# Purity percentage & uncertainty, plus scale uncertainty
col3, col4, col5 = st.columns(3)
with col3:
    purity_percent = parse_float(st.text_input("Purity of chemical [%]", "99.5"))
with col4:
    u_purity_percent = parse_float(st.text_input("Purity uncertainty [%]", "0.05"))
with col5:
    u_scale = parse_float(st.text_input("Scale uncertainty [g]", "0.001"))

# --- Convert to base units ---
vol_L = vol_value / 1000.0 if vol_unit == "mL" else vol_value
u_vol_L = u_vol_value / 1000.0 if u_vol_unit == "mL" else u_vol_value
purity = purity_percent / 100.0
u_purity = u_purity_percent / 100.0

# --- Compute & display ---
if formulas_str and conc > 0 and vol_L > 0 and purity > 0:
    M_total = molar_mass_multi(formulas_str)
    m_req = mass_required(M_total, conc, vol_L, purity)
    u_c_solution = conc_uncertainty(M_total, m_req, u_scale, vol_L, u_vol_L, purity, u_purity)
    rel_u_c_solution = (u_c_solution / conc) * 100

    st.subheader("Results")
    st.write(f"**Total molar mass (all components):** {M_total:.5f} g/mol")
    st.write(f"**Required total mass to weigh:** {m_req:.5f} g")
    st.write(f"**Uncertainty of resulting concentration:** ± {u_c_solution:.6f} mol/L")
    st.write(f"**Relative concentration uncertainty:** {rel_u_c_solution:.3f} %")
    st.caption(
        "Purity affects both the mass calculation and the propagated concentration uncertainty."
    )
