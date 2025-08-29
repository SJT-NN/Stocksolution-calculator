#!/usr/bin/env python3
import math
import streamlit as st
from molmass import Formula
from collections import defaultdict
import pandas as pd
import numpy as np

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

# Solve reverse problem: element mg/L → molecule masses
def solve_masses(element_targets, molecules):
    """
    element_targets: dict {element_symbol: target mg/L}
    molecules: list of (formula, purity_frac)
    Returns: list of masses (g/L) of each molecule to satisfy targets (least squares)
    """
    elems = list(element_targets.keys())
    A = np.zeros((len(elems), len(molecules)))
    for j, (formula, purity) in enumerate(molecules):
        try:
            atom_dict = Formula(formula).atoms
            for i, el in enumerate(elems):
                if el in atom_dict:
                    A[i, j] = atom_dict[el] * Formula(el).mass
        except Exception:
            pass
    # Convert targets to g/L
    b = np.array([element_targets[e] / 1000.0 for e in elems])
    # Adjust by purity
    for j, (_, purity) in enumerate(molecules):
        A[:, j] *= purity
    mols, *_ = np.linalg.lstsq(A, b, rcond=None)
    masses_gL = mols * [Formula(f).mass for f, _ in molecules]
    return masses_gL

# ---------- UI ----------
st.set_page_config(page_title="Multi‑Solute Solution Prep", page_icon="🧪")
st.title("🧪 Multi‑Component Solution Preparation")

mode = st.radio("Select mode:", ["Forward: By Molecule Concentration", "Reverse: By Element Concentration"])

# Volume / Mass selection
use_mass_density = st.checkbox("Enter total solution mass & density instead of volume", value=False)

if not use_mass_density:
    vol_value = parse_float(st.text_input("Solution volume", "1.0"))
    vol_unit = st.selectbox("Volume unit", ["L", "mL"], index=0)
    vol_L = vol_value / 1000.0 if vol_unit == "mL" else vol_value

    u_vol_value = parse_float(st.text_input("Volume uncertainty", "0.001"))
    u_vol_unit = st.selectbox("Uncertainty volume unit", ["L", "mL"], index=0)
    u_vol_L = u_vol_value / 1000.0 if u_vol_unit == "mL" else u_vol_value
else:
    sol_mass_val = parse_float(st.text_input("Total solution mass", "100.0"))
    mass_unit = st.selectbox("Mass unit", ["g", "kg"], index=0)
    sol_mass_g = sol_mass_val * (1000.0 if mass_unit == "kg" else 1.0)

    u_mass_sol_val = parse_float(st.text_input("Scale uncertainty for solution mass", "0.01"))
    u_mass_unit = st.selectbox("Uncertainty mass unit", ["g", "kg"], index=0)
    u_mass_sol = u_mass_sol_val * (1000.0 if u_mass_unit == "kg" else 1.0)

    sol_density_gmL = parse_float(st.text_input("Solution density [g/mL]", "1.00"))

    vol_L = sol_mass_g / (sol_density_gmL * 1000.0)
    u_vol_L = u_mass_sol / (sol_density_gmL * 1000.0)

# ---------- Mode 1: Forward ----------
if mode == "Forward: By Molecule Concentration":
    n_solutes = st.number_input("Number of solutes", min_value=1, value=2, step=1)
    results = []
    element_mgL = defaultdict(float)

    for i in range(int(n_solutes)):
        st.markdown(f"### Solute {i+1}")
        formula = st.text_input(f"Formula {i+1}", "NaCl", key=f"f_{i}")
        conc_val = parse_float(st.text_input(f"Target conc {i+1}", "0.1"), 0.0)
        conc_unit = st.selectbox(f"Concentration unit {i+1}", ["mol/L", "mg/L"], key=f"cu_{i}")
        purity_pct = parse_float(st.text_input(f"Purity [%] {i+1}", "99.5"))
        u_purity_pct = parse_float(st.text_input(f"Purity uncertainty [%] {i+1}", "0.05"))
        u_mass_g = parse_float(st.text_input(f"Scale uncertainty [g] for solute {i+1}", "0.001"))

        if formula and conc_val > 0:
            M = Formula(formula).mass
            conc_molL = conc_val if conc_unit == "mol/L" else mgL_to_molL(conc_val, M)
            purity = purity_pct / 100.0
            u_purity = u_purity_pct / 100.0
            m_req = mass_required_molar(M, conc_molL, vol_L, purity)
            u_c = conc_uncertainty_component(M, m_req, u_mass_g, vol_L, u_vol_L, purity, u_purity)

            results.append({
                "formula": formula,
                "M": M,
                "conc_molL": conc_molL,
                "conc_mgL": molL_to_mgL(conc_molL, M),
                "m_req": m_req,
                "u_c": u_c
            })

            try:
                atoms = Formula(formula).atoms
                for sym, count in atoms.items():
                    element_mgL[sym] += conc_molL * count * Formula(sym).mass * 1000.0
            except Exception:
                pass

    if results:
        st.subheader("Component‑wise Results")
        for r in results:
            st.markdown(f"**{r['formula']}**")
            st.write(f"- Molar mass: {r['M']:.5f} g/mol")
            st.write(f"- Required mass: {r['m_req']:.5f} g")
            st.write(f"- Target concentration: {r['conc_molL']:.6f} mol/L  ({r['conc_mgL']:.3f} mg/L)")
            st.write(f"- Uncertainty in concentration: ± {r['u_c']:.6f} mol/L")

        if element_mgL:
            df_elements = pd.DataFrame(sorted(element_mgL.items(), key=lambda kv: kv[1], reverse=True),
                                       columns=["Element", "Concentration (mg/L)"])
            df_elements["Concentration (mg/L)"] = df_elements["Concentration (mg/L)"].map(lambda x: f"{x:.3f}")
            st.markdown("### Element concentrations in solution (mg/L)")
            st.dataframe(df_elements, use_container_width=True)

# ---------- Mode 2: Reverse ----------
if mode == "Reverse: By Element Concentration":
    # Step 1: Get solutes
    n_solutes = st.number_input("Number of available solutes", min_value=1, value=2, step=1)
    molecules = []
    for i in range(int(n_solutes)):
        col1, col2 = st.columns([2, 1])
        with col1:
            formula = st.text_input(f"Solute {i+1} formula", key=f"sol_formula_{i}")
        with col2:
            purity_pct = parse_float(st.text_input(f"Purity [%] for {formula or f'Solute {i+1}'}", "100.0", key=f"pur_{i}"))
        if formula:
            molecules.append((formula, purity_pct/100.0))

    if molecules:
        # Step 2: Unique element list
        unique_elements = set()
        for formula, _ in molecules:
            try:
                unique_elements.update(Formula(formula).atoms.keys())
            except Exception:
                pass
        unique_elements = sorted(unique_elements)

        # Step 3: Desired concentration for each element
        st.subheader("Target element concentrations")
        element_targets = {}
        for el in unique_elements:
            col1, col2 = st.columns([2, 1])
            with col1:
                conc_value = parse_float(
                    st.text_input(f"{el} target concentration", "0.0", key=f"el_target_val_{el}")
                )
            with col2:
                conc_unit = st.selectbox(
                    "Unit",
                    ["mg/L", "mol/L"],
                    key=f"el_target_unit_{el}"
                )

            # Convert mol/L → mg/L if needed
            if conc_unit == "mol/L":
                conc_value = molL_to_mgL(conc_value, Formula(el).mass)

            element_targets[el] = conc_value  # always stored as mg/L

        # Step 4: Compute masses
        if st.button("Calculate required solute masses"):
            masses_gL = solve_masses(element_targets, molecules)
            st.subheader("Required solute masses")
            for (formula, _), mass in zip(molecules, masses_gL):
                total_mass = mass * vol_L  # scale to actual total solution volume
                st.write(f"- **{formula}**: {mass:.5f} g/L → weigh {total_mass:.5f} g for your batch")
