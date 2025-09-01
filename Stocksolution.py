#!/usr/bin/env python3
import math
import streamlit as st
from molmass import Formula
from collections import defaultdict
import pandas as pd

# ---------- Helpers ----------
def parse_float(x, default=0.0):
    try:
        return float(x)
    except ValueError:
        return default

def mass_required_molar(M, conc_mol_L, vol_L, purity):
    pure_mass = M * conc_mol_L * vol_L
    return pure_mass / purity

def mgL_to_molL(conc_mg_L, M):
    return conc_mg_L / 1000.0 / M

def molL_to_mgL(conc_mol_L, M):
    return conc_mol_L * M * 1000.0

def conc_uncertainty_component(M, m_weighed, u_mass, vol_L, u_vol_L, purity, u_purity):
    dc_dm = purity / (M * vol_L)
    dc_dV = -(m_weighed * purity) / (M * vol_L**2)
    dc_dp = m_weighed / (M * vol_L)
    return math.sqrt((dc_dm * u_mass)**2 + (dc_dV * u_vol_L)**2 + (dc_dp * u_purity)**2)

# ---------- UI ----------
st.set_page_config(page_title="Multi‑Solute Solution Prep", page_icon="🧪")
st.title("🧪 Multi‑Component Solution Preparation")

# Debug mode toggle
debug_mode = st.checkbox("Enable debug mode")

# Volume input
vol_value = parse_float(st.text_input("Solution volume", "1.0"))
vol_unit = st.selectbox("Volume unit", ["L", "mL"], index=0)
vol_L = vol_value / 1000.0 if vol_unit == "mL" else vol_value

u_vol_value = parse_float(st.text_input("Volume uncertainty", "0.001"))
u_vol_unit = st.selectbox("Uncertainty volume unit", ["L", "mL"], index=0)
u_vol_L = u_vol_value / 1000.0 if u_vol_unit == "mL" else u_vol_value

# ---------- Solute Inputs ----------
n_solutes = st.number_input("Number of solutes", min_value=1, value=2, step=1)
results = []
element_mgL_total = defaultdict(float)

for i in range(int(n_solutes)):
    st.markdown(f"### Solute {i+1}")
    formula_raw = st.text_input(f"Formula {i+1}", "NaCl", key=f"f_{i}")
    formula = formula_raw.strip()  # ✅ Strip whitespace
    conc_val = parse_float(st.text_input(f"Target conc {i+1}", "0.1"), 0.0)
    conc_unit = st.selectbox(f"Concentration unit {i+1}", ["mol/L", "mg/L"], key=f"cu_{i}")
    purity_pct = parse_float(st.text_input(f"Purity [%] {i+1}", "99.5"))
    u_purity_pct = parse_float(st.text_input(f"Purity uncertainty [%] {i+1}", "0.05"))
    u_mass_g = parse_float(st.text_input(f"Scale uncertainty [g] for solute {i+1}", "0.001"))

    if formula and conc_val > 0:
        try:
            M = Formula(formula).mass
            conc_molL = conc_val if conc_unit == "mol/L" else mgL_to_molL(conc_val, M)
            purity = purity_pct / 100.0
            u_purity = u_purity_pct / 100.0
            m_req = mass_required_molar(M, conc_molL, vol_L, purity)
            u_c = conc_uncertainty_component(M, m_req, u_mass_g, vol_L, u_vol_L, purity, u_purity)
            conc_mgL_val = molL_to_mgL(conc_molL, M)

            if debug_mode:
                st.write(f"🔍 Debug: Formula={formula}, M={M:.4f}, conc_molL={conc_molL:.6f}, conc_mgL={conc_mgL_val:.3f}")

            # --- Per-element breakdown ---
            element_conc_mgL = {}
            atoms = Formula(formula).atoms
            for sym, count in atoms.items():
                M_elem = Formula(sym).mass
                frac = (count * M_elem) / M
                elem_mgL = conc_mgL_val * frac
                element_conc_mgL[sym] = elem_mgL
                element_mgL_total[sym] += elem_mgL

            results.append({
                "formula": formula,
                "M": M,
                "conc_molL": conc_molL,
                "conc_mgL": conc_mgL_val,
                "m_req": m_req,
                "u_c": u_c,
                "element_conc_mgL": element_conc_mgL
            })

        except Exception as e:
            st.warning(f"⚠️ Could not parse formula '{formula}': {e}")
            if debug_mode:
                st.exception(e)

# ---------- Output ----------
if results:
    if element_mgL_total:
        st.markdown("### 💡 Total Element Concentrations in Solution (mg/L)")
        df_elements = pd.DataFrame(
            sorted(element_mgL_total.items(), key=lambda kv: kv[1], reverse=True),
            columns=["Element", "Concentration (mg/L)"]
        )
        df_elements
