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

def safe_sheet_name(name, max_len=31):
    # Replace * with .
    name = name.replace("*", ".")
    # Trim to Excel's max length
    return name[:max_len]

# ---------- UI ----------
st.set_page_config(page_title="Stock solution prep helper", page_icon="🧪")
st.title("🧪 Stock solution helper")

# Volume input
# ---------- Solution size input ----------
prep_mode = st.radio("Volumetric or gravimetric:",["By volume", "By mass"],index=0)

if prep_mode == "By volume":
    vol_value = parse_float(st.text_input("Solution volume", "1.0"))
    vol_unit = st.selectbox("Volume unit", ["L", "mL"], index=0)
    vol_L = vol_value / 1000.0 if vol_unit == "mL" else vol_value

    u_vol_value = parse_float(st.text_input("Volume uncertainty", "0.001"))
    u_vol_unit = st.selectbox("Uncertainty volume unit", ["L", "mL"], index=0)
    u_vol_L = u_vol_value / 1000.0 if u_vol_unit == "mL" else u_vol_value

elif prep_mode == "By mass":
    mass_value = parse_float(st.text_input("Solution mass", "1000.0"))  # g
    mass_unit = st.selectbox("Mass unit", ["g", "kg"], index=0)
    mass_g = mass_value * (1000.0 if mass_unit == "kg" else 1.0)

    u_mass_value = parse_float(st.text_input("Scale uncertainty", "10"))  # g
    u_mass_unit = st.selectbox("Scale uncertainty unit", ["g", "kg"], index=0)
    u_mass_g = u_mass_value * (1000.0 if mass_unit == "kg" else 1.0)

    density_value = parse_float(st.text_input("Solution density", "1.000"))  # g/mL
    u_density_value = parse_float(st.text_input("Density uncertainty", "0.001"))  # g/mL

    # Convert to volume in litres
    vol_L = (mass_g / density_value) / 1000.0

    # Propagate uncertainty in volume from mass and density uncertainties
    u_vol_L = math.sqrt(
        ((1.0 / density_value) * u_mass_g)**2 +
        ((-mass_g / (density_value**2)) * u_density_value)**2
    ) / 1000.0

# ---------- Number of solutes ----------
n_solutes = st.number_input(
    "Number of solutes",
    min_value=1,
    value= 1,
    step=1
)

# ---------- Data processing ----------
results = []
element_mgL_target = defaultdict(float)
element_mgL_realised = defaultdict(float)

for i in range(int(n_solutes)):
    st.markdown(f"### Solute {i+1}")
    defaults = {}

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

        # NEW: Actual weighed mass input
        actual_mass_g = parse_float(
            st.text_input(f"Actual weighed mass [g] for solute {i+1}", f"{m_req:.5f}", key=f"am_{i}")
        )

        # Realised concentrations
        realised_conc_molL = (actual_mass_g * purity) / (M * vol_L)
        realised_conc_mgL = molL_to_mgL(realised_conc_molL, M)

        u_c = conc_uncertainty_component(M, m_req, u_mass_g, vol_L, u_vol_L, purity, u_purity)
        u_c_realised = conc_uncertainty_component(M, actual_mass_g, u_mass_g, vol_L, u_vol_L, purity, u_purity)

        # Elemental breakdowns
        atoms = getattr(f, "atoms", None)
        if isinstance(atoms, int) or atoms is None:
            if hasattr(f, "composition"):
                atoms = f.composition()
            else:
                raise TypeError(f"Cannot extract element breakdown from {type(f)}")

        element_breakdown_target = {}
        element_breakdown_realised = {}
        for el_key, comp_item in atoms.items():
            sym = getattr(el_key, "symbol", el_key if isinstance(el_key, str) else str(el_key))
            el_mass = getattr(el_key, "mass", None) or Formula(sym).mass
            count_val = getattr(comp_item, "count", comp_item)
            frac_mass = (el_mass * count_val) / M

            elem_conc_target = conc_mgL_val * frac_mass
            elem_conc_realised = realised_conc_mgL * frac_mass

            element_mgL_target[sym] += elem_conc_target
            element_mgL_realised[sym] += elem_conc_realised

            element_breakdown_target[sym] = elem_conc_target
            element_breakdown_realised[sym] = elem_conc_realised

        results.append({
            "formula": raw_formula,  # keep original with '*' if present
            "M": M,
            "conc_molL": conc_molL,
            "conc_mgL": conc_mgL_val,
            "m_req": m_req,
            "u_c": u_c,
            "actual_mass_g": actual_mass_g,
            "realised_conc_molL": realised_conc_molL,
            "realised_conc_mgL": realised_conc_mgL,
            "u_c_realised": u_c_realised,
            "elements_target": element_breakdown_target,
            "elements_realised": element_breakdown_realised
        })

# ---------- Output ----------
if results:
    # Element totals across all solutes
    if element_mgL_target:
        st.markdown("### 💡 Element concentrations in solution (mg/L)")

        # Build both DataFrames
        df_elements_target = pd.DataFrame(
            sorted(element_mgL_target.items(), key=lambda kv: kv[0]),  # sort by element name
            columns=["Element", "Target conc (mg/L)"]
        )
        df_elements_realised = pd.DataFrame(
            sorted(element_mgL_realised.items(), key=lambda kv: kv[0]),
            columns=["Element", "Realised conc (mg/L)"]
        )

        # Merge into one DataFrame on "Element"
        df_elements = pd.merge(df_elements_target, df_elements_realised, on="Element", how="outer")

        # Sort by Target conc (numeric) descending
        df_elements = df_elements.sort_values(by="Target conc (mg/L)", ascending=False)

        # Format numeric columns to 5 decimal places
        for col in ["Target conc (mg/L)", "Realised conc (mg/L)"]:
            df_elements[col] = df_elements[col].map(lambda x: f"{x:.5f}")

        # Show single table
        st.dataframe(df_elements, use_container_width=True)

       

    # Component-wise details
#st.subheader("Component‑wise Results")

for r in results:
    # Table 1: Parameters for this solute
    df_solute = pd.DataFrame([
        ["Molar mass (g/mol)", f"{r['M']:.5f}"],
        ["Target mass (g)", f"{r['m_req']:.5f}"],
        ["Actual mass (g)", f"{r['actual_mass_g']:.5f}"],
        ["Target conc (mol/L)", f"{r['conc_molL']:.6f}"],
        ["Target conc (mg/L)", f"{r['conc_mgL']:.3f}"],
        ["Realised conc (mol/L)", f"{r['realised_conc_molL']:.6f}"],
        ["Realised conc (mg/L)", f"{r['realised_conc_mgL']:.3f}"],
        ["Uncertainty target (mol/L)", f"± {r['u_c']:.6f}"],
        ["Uncertainty realised (mol/L)", f"± {r['u_c_realised']:.6f}"]
    ], columns=["Parameter", "Value"])

 #   st.markdown(f"**{r['formula']}**")
  #  st.table(df_solute)

    # Table 2: Elemental breakdown
    df_elements = pd.DataFrame([
        [elem, f"{r['elements_target'][elem]:.3f}", f"{r['elements_realised'][elem]:.3f}"]
        for elem in r["elements_target"]
    ], columns=["Element", "Target (mg/L)", "Realised (mg/L)"])
   # st.caption("Elemental breakdown")
    #st.table(df_elements)

# --- Combined summary table ---
st.subheader("Summary Table")

df_summary = pd.DataFrame([
    {
        "Formula": r["formula"],
        "Molar mass (g/mol)": f"{r['M']:.5f}",
        "Target mass (g)": f"{r['m_req']:.5f}",
        "Actual mass (g)": f"{r['actual_mass_g']:.5f}",
        "Target conc (mol/L)": f"{r['conc_molL']:.6f}",
        "Target conc (mg/L)": f"{r['conc_mgL']:.3f}",
        "Realised conc (mol/L)": f"{r['realised_conc_molL']:.6f}",
        "Realised conc (mg/L)": f"{r['realised_conc_mgL']:.3f}",
        "Uncertainty target (mol/L)": f"± {r['u_c']:.6f}",
        "Uncertainty realised (mol/L)": f"± {r['u_c_realised']:.6f}"
    }
    for r in results
])
st.markdown("### 📊 Summary")
st.dataframe(df_summary, use_container_width=True)
# Flatten results into rows for export/preview
table_rows = []
for r in results:
   for elem in r["elements_target"]:
        table_rows.append({
            "Formula": raw_formula,
            "Molar mass (g/mol)": M,
            "Target conc (mol/L)": conc_molL,
            "Target conc (mg/L)": conc_mgL_val,
            "Mass required (g)": m_req,
            "Uncertainty target (mol/L)": u_c,
            "Actual mass (g)": actual_mass_g,
            "Realised conc (mol/L)": realised_conc_molL,
            "Realised conc (mg/L)": realised_conc_mgL,
            "Uncertainty realised (mol/L)": u_c_realised})
            
df_all = pd.DataFrame(table_rows)

# Ensure dicts exist to avoid .keys() errors
element_mgL_target = element_mgL_target if isinstance(element_mgL_target, dict) else {}
element_mgL_realised = element_mgL_realised if isinstance(element_mgL_realised, dict) else {}

# Build main DataFrame with uncertainties
df_all = pd.DataFrame(table_rows)  # table_rows should now include u_c and u_c_realised

# Export to Excel
output = BytesIO()
with pd.ExcelWriter(output, engine="openpyxl") as writer:
    # Main results sheet
    df_all.to_excel(writer, index=False, sheet_name="Results")

    # Element totals sheet
    elements = sorted(set(element_mgL_target.keys()) | set(element_mgL_realised.keys()))
    df_totals = pd.DataFrame({
        "Element": elements,
        "Target conc (mg/L)": [element_mgL_target.get(e, 0.0) for e in elements],
        "Realised conc (mg/L)": [element_mgL_realised.get(e, 0.0) for e in elements]
    }).sort_values("Realised conc (mg/L)", ascending=False)
    df_totals.to_excel(writer, index=False, sheet_name="Element Totals")

    # One sheet per solute
    for idx, r in enumerate(results, start=1):
        # Parameters for this solute
        df_solute = pd.DataFrame([
            ["Molar mass (g/mol)", f"{r['M']:.5f}"],
            ["Target mass (g)", f"{r['m_req']:.5f}"],
            ["Actual mass (g)", f"{r['actual_mass_g']:.5f}"],
            ["Target conc (mol/L)", f"{r['conc_molL']:.5f}"],
            ["Target conc (mg/L)", f"{r['conc_mgL']:.5f}"],
            ["Realised conc (mol/L)", f"{r['realised_conc_molL']:.5f}"],
            ["Realised conc (mg/L)", f"{r['realised_conc_mgL']:.5f}"],
            ["Uncertainty target (mol/L)", f"± {r['u_c']:.5f}"],
            ["Uncertainty realised (mol/L)", f"± {r['u_c_realised']:.5f}"]
        ], columns=["Parameter", "Value"])

        # Elemental breakdown
        df_elements = pd.DataFrame([
            [elem, f"{r['elements_target'][elem]:.5f}", f"{r['elements_realised'][elem]:.5f}"]
            for elem in r["elements_target"]
        ], columns=["Element", "Target (mg/L)", "Realised (mg/L)"])

        # Write both tables to the same sheet, one below the other
        sheet_name = f"Solute {idx} - {r['formula']}"
        sheet_name = safe_sheet_name(f"Solute {idx} - {r['formula']}")
        df_solute.to_excel(writer, index=False, sheet_name=sheet_name, startrow=0)
        df_elements.to_excel(writer, index=False, sheet_name=sheet_name, startrow=len(df_solute) + 2)


# Download button
st.download_button(
    label="💾 Download results as Excel",
    data=output.getvalue(),
    file_name="solution_results.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)
    
