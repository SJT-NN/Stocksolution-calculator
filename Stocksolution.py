import csv
import io
import streamlit as st
from molmass import Formula

# --- Core calculation ---
DEFAULT_UNIT = "g/mol"  # options: g/mol, kg/mol, amu

def molar_mass_and_uncertainty(formula_str):
    f = Formula(formula_str)
    molar_mass = f.mass
    uncertainty_sq = 0.0
    for elem, count in f.composition():
        atomic_unc = elem.mass_uncertainty or 0.0
        uncertainty_sq += (count * atomic_unc) ** 2
    uncertainty = uncertainty_sq ** 0.5
    return molar_mass, uncertainty

def convert_units(value, from_unit="g/mol", to_unit="g/mol"):
    if from_unit == to_unit:
        return value
    if from_unit == "g/mol" and to_unit == "kg/mol":
        return value / 1000
    if from_unit == "g/mol" and to_unit == "amu":
        return value / 1.0
    if from_unit == "kg/mol" and to_unit == "g/mol":
        return value * 1000
    if from_unit == "amu" and to_unit == "g/mol":
        return value * 1.0
    raise ValueError(f"Conversion {from_unit} → {to_unit} not supported")

def process_formulas(formulas, unit=DEFAULT_UNIT):
    results = []
    for formula in formulas:
        mm, unc = molar_mass_and_uncertainty(formula)
        results.append({
            "Formula": formula,
            "Molar Mass": convert_units(mm, "g/mol", unit),
            "Uncertainty": convert_units(unc, "g/mol", unit),
            "Unit": unit
        })
    return results

# --- Streamlit UI ---
st.set_page_config(page_title="Batch Molar Mass Calculator", page_icon="🧪")
st.title("🧪 Batch Molar Mass Calculator")
st.markdown("Upload a CSV of chemical formulas or enter one manually to get molar mass and uncertainties.")

# Unit selection
unit = st.selectbox("Output unit", ["g/mol", "kg/mol", "amu"], index=0)

# Manual input
manual_formula = st.text_input("Enter a chemical formula (e.g., H2O, C6H12O6)")

# File upload
uploaded_file = st.file_uploader("Or upload a CSV file (first column = formulas)", type=["csv"])

results = []

if manual_formula:
    results = process_formulas([manual_formula], unit)

elif uploaded_file:
    content = uploaded_file.read().decode("utf-8")
    reader = csv.reader(io.StringIO(content))
    formulas = [row[0].strip() for row in reader if row]
    results = process_formulas(formulas, unit)

# Display results
if results:
    st.subheader("Results")
    st.dataframe(results)

    # Download as CSV
    output_csv = io.StringIO()
    writer = csv.DictWriter(output_csv, fieldnames=["Formula", "Molar Mass", "Uncertainty", "Unit"])
    writer.writeheader()
    writer.writerows(results)

    st.download_button(
        label="Download results as CSV",
        data=output_csv.getvalue(),
        file_name="molar_masses.csv",
        mime="text/csv"
    )
