# Protein Structure Analysis Pipeline

**Description:** This pipeline analyzes protein structures by integrating UniProt and structural data. Given a list of UniProt IDs, it fetches protein metadata (gene, name, cross-references, variants, etc.), downloads experimental PDB structures and AlphaFold predictions, detects structural gaps/coils, aligns predicted models to experiments to patch missing segments, and summarizes everything in a multi-sheet Excel report.

**Project Structure:**

- `pipeline.py`: Main script orchestrating the workflow (parsing args, calling modules, writing Excel).
- `uniprot.py`: Fetches and parses UniProt entries (API queries, caching).
- `alphafold.py`: Fetches AlphaFold models and metadata via the AlphaFold API.
- `pdb_download.py`: Downloads PDB files from RCSB.
- `pdb_details.py`: Retrieves metadata (method, resolution, conformation, entities,etc) from RCSB REST API.
- `structure_analysis.py`: Parses PDB files for helix/sheet info, missing residues, and small molecules.
- `alignment.py`: Aligns AlphaFold models to experimental PDBs (using MDAnalysis) and patches gaps.
- `tests/`: Pytest unit tests for key functions.
- `requirements.txt`: Python dependencies.
- `README.md`: This documentation.

**Installation:**

1. Clone the repository.
2. (Optional) Create a virtual environment.
3. Install dependencies:
    pip install -r requirements.txt
    
**Usage:**

The main script is `pipeline.py`. Prepare an input Excel (e.g. `uniprot_template.xlsx`) with UniProt IDs in the first column (header optional). Then run:

python pipeline.py --input uniprot_template.xlsx --output uniprot_results.xlsx

Additional options:
- `--base PATH`: Base output directory (default: current directory).
- `--window N`: Window size around gaps for alignment (default: 5).
 - `--input`,Input Excel file with UniProt IDs ( default="uniprot_template.xlsx")
- `--output`, Output Excel file (default="uniprot_results.xlsx")
- `--ids`:  One or more UniProt IDs to process directly (bypass Excel template)
- `-v`: Verbose logging.

**Output:**

- `pipeline.log`: Log of all actions and errors.
- `pdb_structures_exp/`: Downloaded experimental PDBs organized by `Gene_UniProt/`.
- `pdb_structures_pred/`: Downloaded AlphaFold PDBs by `Gene_UniProt/`.
- `pdb_structures_af_aligned/`: Aligned and patched PDBs by `Gene/`.
- `uniprot_results.xlsx`: Excel report with sheets:
  - **UniProt Results:** Summary of each UniProt entry (hyperlinked ID and gene, protein name, PDB/AF entries, variants, etc.).
  - **PDB_Details:** For each structure (experimental or predicted): method, resolution,environment exp, conformation (state), polymer entities, organisms, PubMed ref  and listed small molecules (ligands).
  - **Structural_Integrity:** Detected gaps and low-confidence coils for each PDB (type, residue range, length, secondary structure).
  - **Alignment_Summary:** Results of aligning AF models to PDBs Exp (percentage of gaps patched, output file paths, success status).

**Example:**

With `uniprot_template.xlsx` containing:

| UniProt ID |
|------------|
| P12345     |
| Q8N158     |

Run the pipeline:

python pipeline.py --input uniprot_template.xlsx --output results.xlsx

This will fill `results.xlsx` with all data for P12345 and Q8N158, download relevant structures, and align models. Check `pipeline.log` for details. 

Each sheet contains clickable hyperlinks (blue, underlined) to UniProt, RCSB, or AlphaFold entries.

**Notes:**

- The pipeline caches UniProt API responses in `cache/uniprot/` to avoid repeated queries.
- MDAnalysis and BioPython are used for structure handling; ensure they are installed.
