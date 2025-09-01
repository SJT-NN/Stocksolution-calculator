#!/usr/bin/env python3
import math
import streamlit as st
from molmass import Formula
from collections import defaultdict
import pandas as pd
from io import BytesIO
from openpyxl import load_workbook

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

# ---------- Mode selector ----------
mode = st.radio(
    "Select input mode",
    ["Manual entry", "Premade standard solutions"],
    index=0
)

# ---------- Premade multi‑component solutions ----------
standard_solutions = [
    {
        "name": "V3000",
        "components": [
            {
                "formula": "NaCl",
                "conc_val": 0.1000,          # M
                "conc_unit": "mol/L",
                "purity_pct": 100,           # corrected from 1000
                "u_purity_pct": 0.0,
                "u_mass_g": 0.0005
            }
        ]
    },
    {
        "name": "Q26xx",
        "components": [
            {
                "formula": "NaCl",
                "conc_val": 16.48982*1000,
                "conc_unit": "mg/L",
                "purity_pct": 100,
                "u_purity_pct": 0.0,
                "u_mass_g": 0.0005
            },
            {
                "formula": "NaBr",
                "conc_val": 5.796375*1000,
                "conc_unit": "mg/L",
                "purity_pct": 100,
                "u_purity_pct": 0.0,
                "u_mass_g": 0.0005
            },
            {
                "formula": "LiOH*H2O",
                "conc_val": 9.06946*1000,
                "conc_unit": "mg/L",
                "purity_pct": 100,
                "u_purity_pct": 0.0,
                "u_mass_g": 0.0005
            },
            {
                "formula": "H3BO3",
                "conc_val": 5.72034*1000,
                "conc_unit": "mg/L",
                "purity_pct": 100,
                "u_purity_pct": 0.0,
                "u_mass_g": 0.0005
            }
        ]
    },
    {
        "name": "Q31xx",
        "components": [
            {
                "formula": "NaH2PO4",
                "conc_val": 25270,
                "conc_unit": "mg/L",
                "purity_pct": 100,
                "u_purity_pct": 0.0,
                "u_mass_g": 0.0005
            }
        ]
    },
    {
        "name": "Q35xx",
        "components": [
            {
                "formula": "NH4Cl",
                "conc_val": 114570,
                "conc_unit": "mg/L",
                "purity_pct": 100,
                "u_purity_pct": 0.0,
                "u_mass_g": 0.0005
            }
        ]
    }
]

# If premade mode, pick one and set defaults list
if mode == "Premade standard solutions":
    choice = st.selectbox(
        "Choose a standard solution",
        [s["name"] for s in standard_solutions]
    )
    defaults_list = next(s for s in standard_solutions if s["name"] == choice)["components"]
else:
    defaults_list = []

# ---------- Number of solutes ----------
n_solutes = st.number_input(
    "Number of solutes",
    min_value=1,
    value=len(defaults_list) if defaults_list else 2,
    step=1
)

# ---------- Data processing ----------
results = []
element_mgL = defaultdict(float)

for i in range(int(n_solutes)):
    st.markdown(f"### Solute {i+1}")
    defaults = defaults_list[i] if i < len(defaults_list) else {}

    # Allow '*' in formula but ignore for calculations
    raw_formula = st.text_input(f"Formula {i+1}", defaults.get("formula", "NaCl"), key=f"f_{i}")
    formula = raw_formula.replace("*", "")

    conc_val = parse_float(st.text_input(f"Target conc {i+1}", str(defaults.get("conc_val", 0.1)), key=f"tc_{i}"))
    conc_unit = st.selectbox(
        f"Concentration unit {i+1}",
        ["mol/L", "mg/L"],
        index=(["mol/L", "mg/L"].index(defaults.get("conc_unit", "mol/L"))),
        key=f"cu_{i}"
    )
    purity_pct = parse_float(st.text_input(f"Purity [%] {i+1}", str(defaults.get("purity_pct", 100)), key=f"p_{i}"))
    u_purity_pct = parse_float(st.text_input(f"Purity uncertainty [%] {i+1}", str(defaults.get("u_purity_pct", 0.05)), key=f"up_{i}"))
    u_mass_g = parse_float(st.text_input(f"Scale uncertainty [g] for solute {i+1}", str(defaults.get("u_mass_g", 0.001)), key=f"um_{i}"))

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
            "formula": raw_formula,  # keep original with '*' if present
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

    st.subheader("Component‑wise Results")
    for r in results:
        st.markdown(f"**{r['formula']}**")
        st.write(f"- Molar mass: {r['M']:.5f} g/mol")
        st.write(f"- Required mass: {r['m_req']:.5f} g")
        st.write(f"- Target concentration: {r['conc_molL']:.6f} mol/L  ({r['conc_mgL']:.3f} mg/L)")
        st.write(f"- Uncertainty in concentration: ± {r['u_c']:.6f} mol/L")

        if r["elements"]:
            st.write("  **Elemental breakdown (mg/L):**")
            for elem, val in r["elements"].items():
                st.write(f"    - {elem}: {val:.3f}")

if results:
    # Flatten results into rows
    table_rows = []
    for r in results:
        for elem, val in r["elements"].items():
            table_rows.append({
                "Formula": r["formula"],
                "Molar mass (g/mol)": r["M"],
                "Required mass (g)": r["m_req"],
                "Target conc (mol/L)": r["conc_molL"],
                "Target conc (mg/L)": r["conc_mgL"],
                "Uncertainty (mol/L)": r["u_c"],
                "Element": elem,
                "Element conc (mg/L)": val
            })

    df_all = pd.DataFrame(table_rows)

    # Preview table
    st.markdown("### 📊 Component‑wise Results Preview")
    st.dataframe(df_all.style.format({
        "Molar mass (g/mol)": "{:.5f}",
        "Required mass (g)": "{:.5f}",
        "Target conc (mol/L)": "{:.6f}",
        "Target conc (mg/L)": "{:.3f}",
        "Uncertainty (mol/L)": "± {:.6f}",
        "Element conc (mg/L)": "{:.3f}"
    }), use_container_width=True)

    # Export to Excel
    from io import BytesIO
    import openpyxl  # make sure this is installed

    output = BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        df_all.to_excel(writer, index=False, sheet_name="Results")

    st.download_button(
        label="💾 Download results as Excel",
        data=output.getvalue(),
        file_name="solution_results.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )
