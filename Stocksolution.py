import streamlit as st
from molmass import Formula
import math

# --- Core molar mass calculation ---
def molar_mass_and_uncertainty(formula_str):
    f = Formula(formula_str)
    molar_mass = f.mass
    u_sq = 0.0
    for elem, count in f.composition():
        atomic_unc = elem.mass_uncertainty or 0.0
        u_sq += (count * atomic_unc) ** 2
    return molar_mass, math.sqrt(u_sq)

# --- Mass required ---
def mass_required(M, conc, vol):
    return M * conc * vol  # g

# --- Error propagation ---
def propagated_uncertainty(M, u_M, conc, u_conc, vol, u_vol):
    """
    Mass = M * conc * vol
    Total uncertainty combines effects from M, conc, and vol:

    u_m = sqrt( (∂m/∂M * u_M)^2 +
                (∂m/∂c * u_c)^2 +
                (∂m/∂V * u_V)^2 )
    """
    dm_dM = conc * vol
    dm_dc = M * vol
    dm_dV = M * conc
    return math.sqrt(
        (dm_dM * u_M) ** 2 +
        (dm_dc * u_conc) ** 2 +
        (dm_dV * u_vol) ** 2
    )

# --- Streamlit UI ---
#st.set_page_config(page_title="Solution Prep Calculator", page_icon="🧪")
st.title("🧪 Solution Preparation Calculator")
st.markdown("Enter formula, target concentration, and volume to compute required mass — including uncertainty from all three inputs.")

formula = st.text_input("Chemical formula", "NaCl")
conc = st.number_input("Target concentration [mol/L]", value=0.1, min_value=0.0, step=0.01, format="%.4f")
u_conc = st.number_input("Concentration uncertainty [mol/L]", value=0.0005, min_value=0.0, step=0.0001, format="%.5f")
vol = st.number_input("Solution volume [L]", value=1.0, min_value=0.0, step=0.1, format="%.4f")
u_vol = st.number_input("Volume uncertainty [L]", value=0.001, min_value=0.0, step=0.0001, format="%.5f")

if formula and conc > 0 and vol > 0:
    M, u_M = molar_mass_and_uncertainty(formula)
    m = mass_required(M, conc, vol)
    u_m = propagated_uncertainty(M, u_M, conc, u_conc, vol, u_vol)

    st.subheader("Results")
    st.write(f"**Molar mass:** {M:.5f} ± {u_M:.5f} g/mol")
    st.write(f"**Required mass:** {m:.5f} ± {u_m:.5f} g")
    st.caption("Uncertainty is combined from molar mass, concentration, and volume measurements.")
