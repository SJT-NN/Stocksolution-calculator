#!/usr/bin/env python3
"""
Batch molar mass calculator using molmass 2025.4.14

- Reads formulas from a list or CSV file
- Calculates molar mass + combined standard uncertainty
- Supports unit conversion (g/mol, kg/mol, amu)
"""

import csv
from molmass import Formula

# --- Settings ---
DEFAULT_UNIT = "g/mol"  # options: g/mol, kg/mol, amu

# --- Core calculation ---
def molar_mass_and_uncertainty(formula_str):
    """
    Calculate molar mass and combined standard uncertainty.
    """
    f = Formula(formula_str)
    molar_mass = f.mass

    # Quadrature sum of uncertainties from each element
    uncertainty_sq = 0.0
    for elem, count in f.composition():
        atomic_unc = elem.mass_uncertainty or 0.0
        uncertainty_sq += (count * atomic_unc) ** 2

    uncertainty = uncertainty_sq ** 0.5
    return molar_mass, uncertainty

# --- Unit conversion ---
def convert_units(value, from_unit="g/mol", to_unit="g/mol"):
    if from_unit == to_unit:
        return value
    if from_unit == "g/mol" and to_unit == "kg/mol":
        return value / 1000
    if from_unit == "g/mol" and to_unit == "amu":
        return value / 1.0  # 1 g/mol = 1 amu exactly (by definition)
    if from_unit == "kg/mol" and to_unit == "g/mol":
        return value * 1000
    if from_unit == "amu" and to_unit == "g/mol":
        return value * 1.0
    raise ValueError(f"Conversion {from_unit} → {to_unit} not supported")

# --- Batch processing ---
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

# --- CSV helpers ---
def read_formulas_from_csv(file_path):
    formulas = []
    with open(file_path, newline="") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row:
                formulas.append(row[0].strip())
    return formulas

def write_results_to_csv(results, file_path):
    with open(file_path, mode="w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["Formula", "Molar Mass", "Uncertainty", "Unit"])
        writer.writeheader()
        writer.writerows(results)

# --- Main script logic ---
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Batch molar mass calculator using molmass")
    parser.add_argument("input", help="Formula or path to CSV file")
    parser.add_argument("--unit", default=DEFAULT_UNIT, help="Output unit (g/mol, kg/mol, amu)")
    parser.add_argument("--output", help="Optional CSV output path")
    args = parser.parse_args()

    # Detect if input is a CSV file or single formula
    if args.input.lower().endswith(".csv"):
        formulas = read_formulas_from_csv(args.input)
    else:
        formulas = [args.input]

    results = process_formulas(formulas, unit=args.unit)

    # Output to terminal
    for r in results:
        print(f"{r['Formula']}: {r['Molar Mass']:.5f} ± {r['Uncertainty']:.5f} {r['Unit']}")

    # Optionally save to CSV
    if args.output:
        write_results_to_csv(results, args.output)
        print(f"\nResults saved to {args.output}")
