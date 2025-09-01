import math
import streamlit as st
from molmass import Formula
from collections import defaultdict
import pandas as pd

st.title("Molecular Solution Calculator")

formula_str = st.text_input("Enter molecular formula (e.g., NaCl, C6H12O6)", "NaCl")
concentration = st.number_input("Target concentration (mol/L)", min_value=0.0, value=0.1, step=0.01)
volume_L = st.number_input("Volume of stock solution (L)", min_value=0.0, value=1.0, step=0.1)

if st.button("Calculate"):
    try:
        f = Formula(formula_str)
        molar_mass = f.mass  # g/mol

        mass_needed_g = concentration * volume_L * molar_mass

        element_conc = defaultdict(float)
        for elem, count, _ in f.composition():
            count = float(count)  # convert to numeric
            atomic_mass = Formula(elem).mass
            element_conc[elem] = concentration * count * atomic_mass * 1000  # mg/L

        st.subheader("Results")
        st.write(f"**Molar mass of {formula_str}:** {molar_mass:.4f} g/mol")
        st.write(f"**Mass needed:** {mass_needed_g:.4f} g for {volume_L} L at {concentration} mol/L")

        df = pd.DataFrame({
            "Element": list(element_conc.keys()),
            "Concentration (mg/L)": list(element_conc.values())
        })
        st.table(df)

    except Exception as e:
        st.error(f"Error: {e}")
