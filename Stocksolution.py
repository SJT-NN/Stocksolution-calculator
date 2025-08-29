import streamlit as st
import pubchempy as pcp
#from rdkit import Chem
#from rdkit.Chem import Draw
import re
import math

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    rdkit_available = True
except ImportError:
    rdkit_available = False

# --- Atomic weights (g/mol) with uncertainties ---
atomic_weights = {
    "H": (1.00794, 0.00007),
    "C": (12.0107, 0.0008),
    "N": (14.0067, 0.0002),
    "O": (15.9994, 0.0003),
    "Na": (22.98976928, 0.00000002),
    "Cl": (35.453, 0.002),
    "K": (39.0983, 0.0001),
    "S": (32.065, 0.005),
    # Extend as needed...
}

# --- Helpers ---
def parse_formula(formula):
    # Simple regex parser: Element + optional number
    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
    comp = []
    for (elem, count) in tokens:
        count = int(count) if count else 1
        comp.append((elem, count))
    return comp

def molar_mass_and_uncertainty(formula):
    comp = parse_formula(formula)
    total_mass = 0
    total_var = 0
    for elem, n in comp:
        mass, umass = atomic_weights.get(elem, (0,0))
        total_mass += n * mass
        total_var += (n * umass)**2
    u_total = math.sqrt(total_var)
    return total_mass, u_total

def fetch_smiles(name_or_formula):
    try:
        compounds = pcp.get_compounds(name_or_formula, 'name')
        if compounds:
            return compounds[0].isomeric_smiles
    except:
        return None
    return None

# --- Streamlit UI ---
st.title("🧪 Stock Solution Mass Calculator")

chem_input = st.text_input("Chemical name or formula", "NaCl")
molarity = st.number_input("Target molarity (mol/L)", value=1.0, min_value=0.0)
volume = st.number_input("Target volume (L)", value=1.0, min_value=0.0)
purity = st.number_input("Purity (%)", value=100.0, min_value=0.0, max_value=100.0)
u_purity = st.number_input("Purity uncertainty (%)", value=0.0, min_value=0.0)

if chem_input:
    mmass, u_mmass = molar_mass_and_uncertainty(chem_input)
    if mmass > 0:
        # Ideal mass
        ideal_mass = molarity * volume * mmass
        # Adjusted for purity
        adjusted_mass = ideal_mass / (purity/100 if purity else 1)
        # Uncertainty propagation
        rel_u_molar = u_mmass / mmass if mmass else 0
        rel_u_purity = (u_purity/100) / (purity/100) if purity else 0
        total_u = adjusted_mass * math.sqrt(rel_u_molar**2 + rel_u_purity**2)

        st.subheader("Results")
        st.write(f"**Molar mass:** {mmass:.4f} ± {u_mmass:.4f} g/mol")
        st.write(f"Ideal mass: {ideal_mass:.4f} g")
        st.write(f"Mass to weigh (purity adjusted): {adjusted_mass:.4f} g")
        st.write(f"Combined uncertainty: ±{total_u:.4f} g")

        smiles = fetch_smiles(chem_input)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol, size=(300, 300))
            st.image(img, caption="Molecular structure")
        else:
            st.info("Could not fetch molecular structure.")
    else:
        st.error("Unknown elements in formula. Please extend the atomic weight table.")
