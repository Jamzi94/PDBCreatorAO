
# Protein Structure Analysis Pipeline (PDBCreatorAO)

**Description:**  
PDBCreatorAO is an advanced pipeline for protein structure analysis, alignment, patching, and orientation. For any set of UniProt IDs, it automates retrieval of protein metadata (gene names, descriptions, interaction partners), downloads experimental structures from PDB/OPM and AlphaFold predictions, analyzes secondary structures, detects gaps/coils, aligns and patches models, orients structures using DUM references and compiles detailed Excel reports with hyperlinks and quality metrics.

---

## Project Structure

- `pipeline.py`: Main coordinator script (argument parsing, workflow orchestration, Excel/report writing).
- `uniprot.py`: Fetches and parses UniProt entries (gene names, protein descriptions, interaction partners, API queries, caching).
- `pdb_download.py`: Downloads PDB files from RCSB and OPM.
- `alphafold.py`: Fetches AlphaFold models and metadata.
- `structure_analysis.py`: Parses PDBs for helix/sheet info, gaps, low-confidence regions, and small molecules.
- `pdb_details.py`: Retrieves structure metadata (experimental method, resolution, conformational state, experimental environment, polymer entities, ligands, etc.) from RCSB REST API.
- `alignment.py`: Aligns AlphaFold models to experimental PDBs (using MDAnalysis/BioPython), patches missing segments, manages homomeric and heteromeric complexes, detects and resolves chain correspondence.
- `oriented.py`: Handles detection of DUM-containing reference structures and orients all structures for each gene or complex accordingly.
- `requirements.txt`: Python dependencies.
- `tests/`: Pytest unit tests for key functions.
- `README.md`: This documentation.

---

## Installation

1. Clone the repository.
2. *(Optional)* Create a Python virtual environment.
3. Install dependencies:  
   ```bash
   pip install -r requirements.txt
   ```

---

## Pipeline Workflow

1. **Data Retrieval from UniProt (`uniprot.py`):**
    - Fetches protein information (gene names, protein description, interaction partners, related UniProt IDs for heteromers).
    - Extracts all relevant metadata for downstream analysis.
    - Caches responses in `cache/uniprot/` for speed.

2. **Structure Download (`pdb_download.py`, `alphafold.py`):**
    - Downloads experimental structures from RCSB and OPM.
    - Retrieves AlphaFold predictions for each protein.
    - Structures are stored in a directory tree:  
      - `pdb_structures_exp/Gene_UniProt/`
      - `pdb_structures_pred/Gene_UniProt/`

3. **Structure Analysis (`structure_analysis.py`, `pdb_details.py`):**
    - Parses secondary structure annotations (helix, sheet, coil).
    - Detects gaps, missing regions, and low-confidence segments (pLDDT).
    - Lists small molecules/ligands present.
    - Retrieves detailed metadata (method, resolution, environment, polymer entities, conformation/state, source organism, PubMed references).

4. **Structure Alignment and Patching (`alignment.py`):**
    - Aligns AlphaFold models to each experimental structure.
    - Patches missing regions in experimental structures using predicted models.
    - Handles both homomeric (single gene) and heteromeric (multi-gene/partner) complexes.
    - Detects and resolves chain correspondence, selecting optimal mapping (lowest RMSD, valid alignment).

5. **Structure Orientation (`oriented.py`):**
   - Among all available structures for a gene or complex, the pipeline automatically identifies those containing DUM (by analyzing small molecule/ligand annotations).
   - For each group (gene or complex), the first DUM-containing structure is selected as the orientation reference.
   - All other experimental or predicted structures in the group are superimposed onto this reference by precise structural alignment (usually on Cα atoms), ensuring a common spatial frame.
   - The pipeline records RMSD values (before/after alignment) and logs all steps for reproducibility.

6. **Results Compilation (`pipeline.py`):**
    - Generates a comprehensive Excel report with multiple sheets, including clickable hyperlinks to all relevant database entries.
    - Compiles alignment quality metrics, RMSD values, gap-patching summaries, and other analysis results.
    - Organizes all output structures in standardized folders for easy retrieval.

---

## Usage

Prepare an input Excel (e.g. `uniprot_template.xlsx`) with UniProt IDs in the first column (header optional).  
Run the pipeline:

```bash
python pipeline.py 
```

Or with direct UniProt IDs:
```bash
python pipeline.py --ids P12345 P67890
```

Options:
- `--base PATH`: Base output directory (default: current directory)
- `--window N`: Window size for gap alignment (default: 5)
- `--input`: Input Excel file with UniProt IDs (default: "uniprot_template.xlsx")
- `--output`: Output Excel file (default: "uniprot_results.xlsx")
- `--ids`: One or more UniProt IDs to process directly (bypass Excel template)
- `-v`: Verbose logging

---

## Output

- `pipeline.log`: Log of all actions, warnings, and errors.
- `pdb_structures_exp/`: Experimental PDB structures, sorted by `Gene_UniProt/`.
- `pdb_structures_pred/`: AlphaFold models, sorted by `Gene_UniProt/`.
- `pdb_structures_af_aligned/`: Aligned and patched structures by `Gene/`.
- `pdb_structures_oriented/`: Oriented structures (after DUM alignment), by `Gene/`.
- `uniprot_results.xlsx`: Multi-sheet Excel report with:
  - **UniProt Results:** Summary of protein entries (hyperlinks, partners, cross-references, variants).
  - **PDB_Details:** Structure info (method, resolution, state, environment, polymer entities, source organism, PubMed, ligands).
  - **Structural_Integrity:** Gaps and coils detected in each structure (type, range, secondary structure context).
  - **Alignment_Summary:** Details of alignments (percentage patched, paths to outputs, patching success, errors).
  - **Alignment_RMSD:** RMSD values and orientation metrics for each aligned structure.

All sheets have clickable hyperlinks (UniProt, RCSB, AlphaFold, PubMed…).

  
---

## Warning 

All the patched structures need to be minimized in vacuum prior to their usage


---

## Notes

- UniProt API results are cached locally for efficiency.
- Structure alignment uses MDAnalysis and BioPython—**make sure they are installed**.
- The pipeline is robust for both single chains and multi-chain complexes, handling complicated chain correspondences and orientation requirements.
- Each output directory is named systematically for downstream compatibility and traceability.
- All steps are logged in `pipeline.log` for transparency and reproducibility.

---

# Version Française

**Description**  
PDBCreatorAO est un pipeline avancé pour l’analyse structurale, l’alignement, le raccommodage et l’orientation des protéines. À partir d’IDs UniProt, il automatise la récupération des métadonnées protéiques (noms de gènes, descriptions, partenaires d’interaction), le téléchargement des structures expérimentales (PDB/OPM) et des prédictions AlphaFold, l’analyse des structures secondaires, la détection des lacunes et régions de faible confiance, l’alignement et le raccommodage, l’orientation via références DUM et la génération de rapports Excel multi-feuilles détaillés.

---

## Structure du projet :

- `pipeline.py` : Coordination principale (arguments, workflow, rapports Excel)
- `uniprot.py` : Récupération/parsing UniProt (gènes, partenaires, requêtes API, cache)
- `pdb_download.py` : Téléchargement de PDB depuis RCSB/OPM
- `alphafold.py` : Récupération des modèles AlphaFold et métadonnées
- `structure_analysis.py` : Parsing PDB (hélice, feuillet, lacune, coil, petites molécules)
- `pdb_details.py` : Métadonnées structurales (méthode, résolution, état conformationnel, environnement, entités polymères, ligands…)
- `alignment.py` : Alignement, raccommodage (gestion homomère/hétéromère, correspondance de chaînes)
- `oriented.py` : Détection des structures de référence DUM et orientation systématique de toutes les structures du même gène/complexe.
- `requirements.txt` : Dépendances Python
- `tests/` : Tests unitaires Pytest
- `README.md` : Cette documentation

---

## Installation :

1. Cloner le dépôt.
2. *(Optionnel)* Créer un environnement virtuel.
3. Installer les dépendances :
   ```bash
   pip install -r requirements.txt
   ```

---

## Flux de travail :

1. **Récupération des données UniProt (`uniprot.py`)** :
   - Récupère infos protéiques (noms de gènes, description, partenaires d’interaction, complexes hétéromériques).
   - Extrait toutes les métadonnées pertinentes pour l’analyse.
   - Stocke en cache dans `cache/uniprot/`.

2. **Téléchargement de structures (`pdb_download.py`, `alphafold.py`)** :
   - Télécharge structures expérimentales (RCSB/OPM).
   - Récupère modèles AlphaFold.
   - Arborescence organisée :
     - `pdb_structures_exp/Gene_UniProt/`
     - `pdb_structures_pred/Gene_UniProt/`

3. **Analyse structurale (`structure_analysis.py`, `pdb_details.py`)** :
   - Extrait hélices, feuillets, coils, gaps/lacunes, molécules/ligands.
   - Récupère métadonnées détaillées (méthode, résolution, état, environnement, entités, organismes, PubMed).

4. **Alignement et raccommodage (`alignment.py`)** :
   - Aligne les modèles AlphaFold sur chaque structure expérimentale.
   - Raccommode les régions manquantes (gaps).
   - Gère les complexes homomériques et hétéromériques.
   - Sélectionne automatiquement la correspondance de chaînes la plus optimale (RMSD minimal).

5. **Orientation des structures (`oriented.py`)** :
   - Parmi toutes les structures disponibles pour un gène ou complexe, le pipeline identifie automatiquement celles contenant DUM (en analysant les petites molécules/ligands de chaque structure).
   - Pour chaque groupe (gène ou complexe), la première structure contenant DUM sert de référence d’orientation.
   - Toutes les autres structures expérimentales ou prédictives du groupe sont superposées sur cette référence par alignement structurale précis (généralement sur les atomes Cα), garantissant un repère spatial commun.
   - Le pipeline enregistre les valeurs de RMSD (avant/après alignement) et trace toutes les étapes pour assurer la reproductibilité.

6. **Compilation des résultats (`pipeline.py`)** :
   - Génère un rapport Excel multi-feuilles avec liens cliquables.
   - Résume les métriques d’alignement, patching, RMSD, gaps, etc.
   - Organise tous les outputs dans des dossiers standardisés.

---

## Utilisation :

Préparer un fichier Excel d’entrée (ex : `uniprot_template.xlsx`) avec les UniProt IDs en première colonne (en-tête optionnel).  
Exécuter :

```bash
python pipeline.py 
```

Ou pour des IDs directs :
```bash
python pipeline.py --ids P12345 P67890
```

Options :
- `--base PATH` : Répertoire de sortie (défaut : courant)
- `--window N` : Fenêtre autour des gaps (défaut : 5)
- `--input` : Fichier Excel d’entrée (défaut : "uniprot_template.xlsx")
- `--output` : Fichier Excel de sortie (défaut : "uniprot_results.xlsx")
- `--ids` : IDs à traiter (ignore l’Excel)
- `-v` : Log verbeux

---

## Sortie :

- `pipeline.log` : Log détaillé de l’exécution
- `pdb_structures_exp/` : Structures expérimentales, par `Gene_UniProt/`
- `pdb_structures_pred/` : Structures AlphaFold, par `Gene_UniProt/`
- `pdb_structures_af_aligned/` : Structures alignées/raccommodées par `Gene/`
- `pdb_structures_oriented/` : Structures orientées via DUM par `Gene/`
- `uniprot_results.xlsx` : Rapport multi-feuilles :
  - **UniProt Results** : Synthèse des entrées UniProt (liens, partenaires, cross-ref, variants)
  - **PDB_Details** : Infos structurales (méthode, résolution, état, environnement, entités, organisme, PubMed, ligands)
  - **Structural_Integrity** : Gaps/coils par structure (type, plage, contexte secondaire)
  - **Alignment_Summary** : Statuts de l’alignement, pourcentages patchés, chemins, erreurs
  - **Alignment_RMSD** : RMSD des structures orientées

Toutes les feuilles incluent des liens hypertextes.

---

## Notes :

- Les requêtes UniProt sont mises en cache dans `cache/uniprot/`.
- Utilise MDAnalysis et BioPython pour la manipulation structurale.
- Compatible mono-chaîne et multi-chaînes, gestion complexe des correspondances.
- Sorties nommées systématiquement pour la reproductibilité.
- Pipeline documenté via `pipeline.log` à chaque étape.
