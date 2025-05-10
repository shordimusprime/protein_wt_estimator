import streamlit as st
from collections import Counter
from math import floor
import io
import pandas as pd 

# Amino acid and PTM weights
amino_acid_wt = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
    'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
    'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
    'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
}

ptm_wt = {
    "Phosphorylation": 79.97,
    "Acetylation": 42.01,
    "Methylation": 14.02,
    "Ubiquitination": 8564.8,
    "Glycosylation": 203.1,
    "Disulfide bond": -2.015  # per disulfide bond
}

def calculate_weight(sequence, ptms):
    sequence = sequence.upper()
    valid_seq = ''.join([aa for aa in sequence if aa in amino_acid_wt])
    invalid = [aa for aa in sequence if aa not in amino_acid_wt]

    weight = sum(amino_acid_wt[aa] for aa in valid_seq)
    composition = Counter(valid_seq)

    if len(valid_seq) > 1:
        weight -= (len(valid_seq) - 1) * 18.015  # Water loss from peptide bonds

    ptm_applied = []
    for ptm in ptms:
        if ptm == "Disulfide bond" and 'C' in composition:
            bonds = floor(composition['C'] / 2)
            weight += ptm_wt[ptm] * bonds
            ptm_applied.append((ptm, bonds))
        else:
            weight += ptm_wt[ptm]
            ptm_applied.append((ptm, 1))

    return weight, composition, invalid, ptm_applied

def parse_fasta(file):
    sequences = []
    current_seq = ""
    for line in file:
        line = line.decode("utf-8").strip()
        if line.startswith(">"):
            if current_seq:
                sequences.append(current_seq)
                current_seq = ""
        else:
            current_seq += line
    if current_seq:
        sequences.append(current_seq)
    return sequences

# Streamlit UI
st.title("🔬 Protein Molecular Weight Calculator")
st.markdown("Paste a protein sequence or upload a `.txt` or `.fasta` file. You can also select post-translational modifications (PTMs).")

# Upload option
uploaded_file = st.file_uploader("Upload a protein sequence file (txt or FASTA):", type=["txt", "fasta"])

# Manual sequence input
sequence_input = st.text_area("Or manually enter a protein sequence (1-letter code):", height=150)

ptm_selected = st.multiselect("Select Post-Translational Modifications (PTMs):", options=list(ptm_wt.keys()))

# Button
if st.button("Calculate"):
    sequences_to_process = []

    if uploaded_file:
        content = uploaded_file.readlines()
        if uploaded_file.name.endswith(".fasta"):
            sequences_to_process = parse_fasta(content)
        else:
            sequences_to_process = [line.decode("utf-8").strip() for line in content if line.strip()]
    elif sequence_input:
        sequences_to_process = [sequence_input.strip()]
    else:
        st.error("Please provide a sequence or upload a file.")
    
    results =[]

    for idx, seq in enumerate(sequences_to_process, start=1):
        weight, composition, invalid, ptms_applied = calculate_weight(seq, ptm_selected)

        st.subheader(f"Result for Sequence {idx}:")
        st.success(f"**{weight:.2f} Da**  (**{weight / 1000:.2f} kDa**)")

        st.markdown("**Amino Acid Composition**")
        st.write(dict(composition))

        if ptms_applied:
            st.markdown("**Applied PTMs:**")
            for ptm, count in ptms_applied:
                delta = ptm_wt[ptm] * count
                st.write(f"• {ptm} (×{count}) → {delta:+.2f} Da")

        if invalid:
            st.warning(f"Ignored invalid characters: {', '.join(invalid)}")

        #prepare data to be downloaded in csv format
        ptm_str = "; ".join([f"{ptm} x{count}" for ptm, count in ptms_applied]) if ptms_applied else "None"
        comp_str = "; ".join([f"{aa}:{count}" for aa, count in composition.items()])
        invalid_str = ','.join(invalid) if invalid else ""

        results.append({
            "Sequence": seq,
            "Molecular Weight (Da)": round(weight, 2),
            "Molecular Weight (kDa)": round(weight / 1000, 2),
            "Amino Acid Composition": comp_str,
            "Applied PTMs": ptm_str,
            "Invalid Characters": invalid_str
        })

    # Convert results to DataFrame and offer CSV download
    if results:
        df = pd.DataFrame(results)
        csv = df.to_csv(index=False)

        st.download_button(
            label="📥 Download Results as CSV",
            data=csv,
            file_name="protein_weight_results.csv",
            mime="text/csv"
        )
