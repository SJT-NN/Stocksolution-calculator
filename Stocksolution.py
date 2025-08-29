import streamlit as st
from molmass import Formula
import math

# --- Molar mass & its uncertainty from molmass ---
def molar_mass_and_uncertainty(formula_str):
    f = Formula(formula_str)
    molar_mass = f.mass
    u_sq = 0.0
    for elem, count in f.composition():
        atomic_unc = elem.mass_uncertainty or 0.0
        u_sq += (count * atomic_unc) ** 2
    return molar_mass, math.sqrt(u_sq)

# --- Required mass, purity‑corrected ---
def mass_required(M, conc, vol_L, purity):
    """
    M in g/mol, conc in mol/L, vol_L in L, purity as fraction (0–1)
    """
    pure_mass = M * conc * vol_L
    return pure_mass / purity  # weigh more if purity < 100%

# --- Full uncertainty propagation ---
def propagated_uncertainty(M, u_M, conc, u_conc, vol_L, u_vol_L, purity, u_purity, u_scale):
    """
    m = (M * c * V) / p
    Inputs:
        M: molar mass
        u_M: molar mass uncertainty
        c: concentration
        u_c: concentration uncertainty
        V: volume in L
        u_V: volume uncertainty in L
        p: purity (fraction, e.g. 0.995)
        u_purity: purity uncertainty (fraction)
        u_scale: weighing device uncertainty in g (independent term)
    """
    dm_dM = (conc * vol_L) / purity
    dm_dc = (M * vol_L) / purity
    dm_dV = (M * conc) / purity
    dm_dp = -(M * conc * vol_L) / (purity ** 2)

    # Combine in quadrature: from formula parameters + independent scale noise
    return math.sqrt(
        (dm_dM * u_M) ** 2 +
        (dm_dc * u_conc) ** 2 +
        (dm_dV * u_vol_L) ** 2 +
        (dm_dp * u_purity) ** 2 +
        (u_scale) ** 2
    )

# --- Streamlit UI ---
st.set_page_config(page_title="Solution Prep Calculator", page_icon="🧪")
st.title("🧪 Solution Preparation Calculator (with full uncertainty)")

formula = st.text_input("Chemical formula", "NaCl")

col1, col2 = st.columns(2)
with col1:
    conc = st.number_input("Target concentration", value=0.1, min_value=0.0, step=0.01, format="%.4f")
    u_conc = st.number_input("Concentration uncertainty", value=0.0005, min_value=0.0, step=0.0001, format="%.5f")
with col2:
    conc_unit = st.selectbox("Concentration unit", ["mol/L"], index=0)  # easy to expand later

# Volume entry with units
vol_value = st.number_input("Solution volume", value=100.0, min_value=0.0, step=1.0, format="%.3f")
vol_unit = st.selectbox("Volume unit", ["mL", "L"], index=0)
u_vol_value = st.number_input("Volume uncertainty", value=0.1, min_value=0.0, step=0.01, format="%.4f")
u_vol_unit = st.selectbox("Uncertainty volume unit", ["mL", "L"], index=0)

# Purity & scale
purity_percent = st.number_input("Purity of chemical [%]", value=99.5, min_value=0.0, max_value=100.0, step=0.1, format="%.3f")
u_purity_percent = st.number_input("Purity uncertainty [%]", value=0.05, min_value=0.0, step=0.01, format="%.4f")
u_scale = st.number_input("Scale uncertainty [g]", value=0.001, min_value=0.0, step=0.0001, format="%.5f")

# Convert all units to base SI for calc
vol_L = vol_value / 1000.0 if vol_unit == "mL" else vol_value
u_vol_L = u_vol_value / 1000.0 if u_vol_unit == "mL" else u_vol_value
purity = purity_percent / 100.0
u_purity = u_purity_percent / 100.0

# Run calculations
if formula and conc > 0 and vol_L > 0 and purity > 0:
    M, u_M = molar_mass_and_uncertainty(formula)
    m_req = mass_required(M, conc, vol_L, purity)
    u_m_req = propagated_uncertainty(M, u_M, conc, u_conc, vol_L, u_vol_L, purity, u_purity, u_scale)

    st.subheader("Results")
    st.write(f"**Molar mass:** {M:.5f} ± {u_M:.5f} g/mol")
    st.write(f"**Required mass to weigh:** {m_req:.5f} ± {u_m_req:.5f} g")
    st.caption("Uncertainty combines molar mass, concentration, volume, purity, and scale tolerance.")
