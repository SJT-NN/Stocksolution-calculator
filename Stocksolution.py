#!/usr/bin/env python3
import math
import streamlit as st
from molmass import Formula

# --- Molar mass (no atomic weight uncertainty) ---
def molar_mass(formula_str):
    f = Formula(formula_str)
    return f.mass

# --- Required mass, purity‑corrected ---
def mass_required(M, conc, vol_L, purity):
    pure_mass = M * conc * vol_L
    return pure_mass / purity

# --- Uncertainty of concentration from uncertainties in mass, vol, purity ---
def conc_uncertainty(M, m_weighed, u_mass,
                     vol_L, u_vol_L,
                     purity, u_purity):
    """
    Standard uncertainty in concentration [mol/L] from
    weighing (u_mass), volume delivery (u_vol_L), and purity (u_purity)
    """
    # c = (m * p) / (M * V)
    dc_dm = purity / (M * vol_L)                   # mass term
    dc_dV = -(m_weighed * purity) / (M * vol_L**2) # volume term
    dc_dp = m_weighed / (M * vol_L)                # purity term

    return math.sqrt(
        (dc_dm * u_mass) ** 2 +     # from scale
        (dc_dV * u_vol_L) ** 2 +    # from volume
        (dc_dp * u_purity) ** 2     # from purity
    )

# --- Safe float parser ---
def parse_float(text, default=0.0):
    try:
        return float(text)
    except ValueError:
        return default

# --- UI ---
st.set_page_config(page_title="Solution Prep Calculator", page_icon="🧪")
st.title("🧪 Solution Preparation Calculator")
st.markdown(
    "Enter your target solution parameters — get the required mass, "
    "the uncertainty of the resulting concentration, and its relative %."
)

# Formula
formula = st.text_input("Chemical formula", "NaCl")

# Target concentration (absolute target uncertainty kept for display only)
col1, col2 = st.columns(2)
with col1:
    conc = parse_float(st.text_input("Target concentration [mol/L]", "0.1"))
with col2:
    u_conc_input = parse_float(st.text_input("Target conc uncertainty [mol/L]", "0.0005"))

# Volume and units
col3, col4 = st.columns(2)
with col3:
    vol_value = parse_float(st.text_input("Solution volume", "100.0"))
    vol_unit = st.selectbox("Volume unit", ["mL", "L"], index=0)
with col4:
    u_vol_value = parse_float(st.text_input("Volume uncertainty", "0.1"))
    u_vol_unit = st.selectbox("Uncertainty volume unit", ["mL", "L"], index=0)

# Purity & scale uncertainties
col5, col6 = st.columns(2)
with col5:
    purity_percent = parse_float(st.text_input("Purity of chemical [%]", "99.5"))
with col6:
    u_purity_percent = parse_float(st.text_input("Purity uncertainty [%]", "0.05"))

u_scale = parse_float(st.text_input("Scale uncertainty [g]", "0.001"))

# --- Convert to base units ---
vol_L = vol_value / 1000.0 if vol_unit == "mL" else vol_value
u_vol_L = u_vol_value / 1000.0 if u_vol_unit == "mL" else u_vol_value
purity = purity_percent / 100.0
u_purity = u_purity_percent / 100.0

# --- Compute & display ---
if formula and conc > 0 and vol_L > 0 and purity > 0:
    M = molar_mass(formula)
    m_req = mass_required(M, conc, vol_L, purity)
    u_c_solution = conc_uncertainty(M, m_req, u_scale, vol_L, u_vol_L, purity, u_purity)
    rel_u_c_solution = (u_c_solution / conc) * 100

    st.subheader("Results")
    st.write(f"**Molar mass:** {M:.5f} g/mol")
    st.write(f"**Required mass to weigh:** {m_req:.5f} g")
    st.write(f"**Uncertainty of resulting concentration:** ± {u_c_solution:.6f} mol/L")
    st.write(f"**Relative concentration uncertainty:** {rel_u_c_solution:.3f} %")
    st.caption(
        "Concentration uncertainty includes contributions from weighing, "
        "volume, and purity uncertainties")
    
