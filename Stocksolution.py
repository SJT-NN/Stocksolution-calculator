#!/usr/bin/env python3
import math
import streamlit as st
from molmass import Formula
from collections import defaultdict
import pandas as pd
from openpyxl import load_workbook
from io import BytesIO

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

# Volume input
vol_value = parse_float(st.text_input("Solution volume", "1.0"))
vol_unit = st.selectbox("Volume unit", ["L", "mL"], index=0)
vol_L = vol_value / 1000.0 if vol_unit == "mL" else vol_value

u_vol_value = parse_float(st.text_input("Volume uncertainty", "0.001"))
u_vol_unit = st.selectbox("Uncertainty volume unit", ["L", "mL"], index=0)
u_vol_L = u_vol_value / 1000.0 if u_vol_unit == "mL" else u_vol_value

# ---------- NEW: Mode selector ----------
mode = st.radio(
    "Select input mode",
    ["Manual entry", "Premade standard solutions"],
    index=0
)

# ---------- NEW: Define premade solutions ----------
standard_solutions = [
    {
        "name": "0.1 M NaCl",
        "formula": "NaCl",
        "conc_val": 0.1,
        "conc_unit": "mol/L",
        "purity_pct": 100,
        "u_purity_pct": 0.05,
        "u_mass_g": 0.001
    },
    {
        "name": "1000 mg/L KNO3",
        "formula": "KNO3",
        "conc_val": 1000,
        "conc_unit": "mg/L",
        "purity_pct": 99.5,
        "u_purity_pct": 0.1,
        "u_mass_g": 0.001
    }
]

# ---------- Processing ----------
results = []
element_mgL = defaultdict(float)

if mode == "Manual entry":
    n_solutes = st.number_input("Number of solutes", min_value=1, value=2, step=1)
    for i in range(int(n_solutes)):
        st.markdown(f"### Solute {i+1}")
        formula = st.text_input(f"Formula {i+1}", "NaCl", key=f"f_{i}")
        conc_val = parse_float(st.text_input(f"Target conc {i+1}", "0.1"), 0.0)
        conc_unit = st.selectbox(f"Concentration unit {i+1}", ["mol/L", "mg/L"], key=f"cu_{i}")
        purity_pct = parse_float(st.text_input(f"Purity [%] {i+1}", "100"))
        u_purity_pct = parse_float(st.text_input(f"Purity uncertainty [%] {i+1}", "0.05"))
        u_mass_g = parse_float(st.text_input(f"Scale uncertainty [g] for solute {i+1}", "0.001"))

        if formula and conc_val > 0:
            f = Formula(formula)
            M = f.mass
            conc_molL = conc_val if conc_unit == "mol/L" else mgL_to_molL(conc_val, M)
            conc_mgL_val = molL_to_mgL(conc_molL, M)
            purity = purity_pct / 100.0
            u_purity = u_purity_pct / 100.0
            m_req = mass_required_molar(M, conc_molL, vol_L, purity)
            u_c = conc_uncertainty_component(M, m_req, u_mass_g, vol_L, u_vol_L, purity, u_purity)

            # Elemental breakdown
            atoms = getattr(f, "atoms", None)
            if isinstance(atoms, int) or atoms is None:
                if hasattr(f, "composition"):
                    atoms = f.composition()
                else:
                    raise TypeError(f"Cannot extract element breakdown from {type(f)}")

            element_breakdown = {}
            for el_key, comp_item in atoms.items():
                sym = getattr(el_key, "symbol", el_key if isinstance(el_key, str) else str(el_key))
                el_mass = getattr(el_key, "mass", None)
                if el_mass is None:
                    el_mass = Formula(sym).mass
                count_val = getattr(comp_item, "count", comp_item)
                frac_mass = (el_mass * count_val) / M
                elem_conc_mgL = conc_mgL_val * frac_mass

                element_mgL[sym] += elem_conc_mgL
                element_breakdown[sym] = elem_conc_mgL

            results.append({
                "formula": formula,
                "M": M,
                "conc_molL": conc_molL,
                "conc_mgL": conc_mgL_val,
                "m_req": m_req,
                "u_c": u_c,
                "elements": element_breakdown
            })

else:  # Premade standard solutions
    choice = st.selectbox(
        "Choose a standard solution",
        [s["name"] for s in standard_solutions]
    )
    selected = next(s for s in standard_solutions if s["name"] == choice)

    formula = selected["formula"]
    conc_val = selected["conc_val"]
    conc_unit = selected["conc_unit"]
    purity_pct = selected["purity_pct"]
    u_purity_pct = selected["u_purity_pct"]
    u_mass_g = selected["u_mass_g"]

    f = Formula(formula)
    M = f.mass
    conc_molL = conc_val if conc_unit == "mol/L" else mgL_to_molL(conc_val, M)
    conc_mgL_val = molL_to_mgL(conc_molL, M)
    purity = purity_pct / 100.0
    u_purity = u_purity_pct / 100.0
    m_req = mass_required_molar(M, conc_molL, vol_L, purity)
    u_c = conc_uncertainty_component(M, m_req, u_mass_g, vol_L, u_vol_L, purity, u_purity)

    atoms = getattr(f, "atoms", None)
    if isinstance(atoms, int) or atoms is None:
        if hasattr(f, "composition"):
            atoms = f.composition()
        else:
            raise TypeError(f"Cannot extract element breakdown from {type(f)}")

    element_breakdown = {}
    for el_key, comp_item in atoms.items():
        sym = getattr(el_key, "symbol", el_key if isinstance(el_key, str) else str(el_key))
        el_mass = getattr(el_key, "mass", None)
        if el_mass is None:
            el_mass = Formula(sym).mass
        count_val = getattr(comp_item, "count", comp_item)
        frac_mass = (el_mass * count_val) / M
        elem_conc_mgL = conc_mgL_val * frac_mass

        element_mgL[sym] += elem_conc_mgL
        element_breakdown[sym] = elem_conc_mgL

    results.append({
        "formula": formula,
        "M": M,
        "conc_molL": conc_molL,
        "conc_mgL": conc_mgL_val,
        "m_req": m_req,
        "u_c": u_c,
        "elements": element_breakdown
    })

# ---------- Output ----------
if results:
    if element_mgL:
        st.markdown("### 💡 Element concentrations in solution (mg/L)")
        df_elements = pd.DataFrame(
            sorted(element_mgL.items(), key=lambda kv: kv[1], reverse=True),
            columns=["Element", "Concentration (mg/L)"]
        )
        df_elements["Concentration (mg/L)"] = df_elements["Concentration (mg/L)"].map(lambda x: f"{x:.3f}")
        st.dataframe(df_elements, use_container_width=True)

    st.subheader("Component‑wise Results
