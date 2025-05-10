# Protein Molecular Weight Estimator

This web application calculates the molecular weight of a protein based on its amino acid sequence and optional post-translational modifications (PTMs). The app is built using **Streamlit** and allows users to input protein sequences manually or upload a sequence file in `.txt` or `.fasta` format. The app also supports various PTMs, including phosphorylation, acetylation, methylation, ubiquitination, glycosylation, and disulfide bond formation.

## Features

- **Manual Input**: Paste a protein sequence using 1-letter amino acid codes.
- **File Upload**: Upload a `.txt` or `.fasta` file with protein sequences.
- **Post-Translational Modifications (PTMs)**: Select from a list of common PTMs, including disulfide bond formation.
- **Molecular Weight Calculation**: Calculate the molecular weight in Da (Daltons) and kDa (kilodaltons).
- **Amino Acid Composition**: View the composition of amino acids in the protein sequence.
- **Result Download**: Download the results as a CSV file.

## Installation

To run this app locally, follow these steps:

1. Clone the repository:
    ```bash
    git clone https://github.com/shordimusprime/protein_wt_estimator.git
    cd protein_wt_estimator
    ```

2. Create a virtual environment (optional but recommended):
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```

3. Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Running the App Locally

After setting up the environment, run the Streamlit app locally:

```bash
streamlit run app.py
