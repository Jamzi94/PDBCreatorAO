#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Multi-reference alignment module for PDB structures
# ---------------------------------------------------------------------------
"""
oriented.py

Multi-reference structural orientation using OPM/PDBTM databases.

v2 Features:
- Dynamic OPM availability checking via API
- PDBTM fallback for orientation data
- RCSB Search API for finding homolog templates
- ADME family classification by gene prefix
- Hierarchical file organization: Family/Gene/Conformation/PDB_ID

Key features:
- Groups structures by gene name or UniProt ID
- Identifies OPM-oriented structures (DUM markers) as references
- Dynamically searches for orientation templates if none locally available
- Aligns all structures in a group to their common reference
- Tracks RMSD before and after alignment
- Saves results to Excel for quality control

The orientation process ensures that all structures of the same protein
share a common spatial reference frame, which is essential for comparative
analysis and visualization.
"""

import logging
import re
import sys
import shutil
import time
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Any

import requests
import numpy as np
import pandas as pd
import openpyxl
from Bio.PDB import PDBParser, Superimposer, PDBIO

# Try to import BeautifulSoup for HTML parsing (optional for PDBTM)
try:
    from bs4 import BeautifulSoup
    HAS_BS4 = True
except ImportError:
    HAS_BS4 = False

# --- Constants & regular expressions ---
# Regular expression to detect DUM (dummy atom) markers in small molecule lists
DUM_RE = re.compile(r"\bDUM\b", re.IGNORECASE)
# Regular expression to extract PDB IDs from Excel HYPERLINK formulas
HYP_RE = re.compile(
    r'=?HYPERLINK\(\s*"[^"]*"\s*[;,]\s*"(?P<id>[A-Za-z0-9\-]+)"\s*\)',
    re.IGNORECASE
)

# --- API URLs ---
OPM_DOWNLOAD_URL = "https://biomembhub.org/shared/opm-assets/pdb/{pdb_id}.pdb"
OPM_PROTEIN_PAGE = "https://opm.phar.umich.edu/proteins/{pdb_id}"
PDBTM_PAGE = "http://pdbtm.enzim.hu/?do=display&id={pdb_id}"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

# --- ADME Family Gene Prefixes ---
ADME_GENE_PREFIXES = {
    "ABC_efflux": ["ABCB", "ABCG", "ABCC", "MDR", "MRP", "BCRP", "P-GP"],
    "SLC22": ["SLC22", "OCT", "OAT"],
    "SLCO_OATP": ["SLCO", "OATP"],
    "SLC15_POT": ["SLC15", "PEPT", "PHT"],
    "SLC47_MATE": ["SLC47", "MATE"],
    "SLC10": ["SLC10", "NTCP", "ASBT", "SOAT"],
    "SLC29_ENT": ["SLC29", "ENT"],
}

# --- HTTP settings ---
HTTP_TIMEOUT = 10
HTTP_RETRY_DELAY = 1.0

# --- Logging configuration ---
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Utility functions for data extraction and processing
# ---------------------------------------------------------------------------

def extract_pdb_from_cell(cell) -> str:
    """
    Extract PDB identifier from an Excel cell containing hyperlink or plain text.
    
    Handles multiple cell formats:
    - Excel hyperlink objects
    - HYPERLINK formula strings
    - Plain text PDB IDs
    
    Parameters
    ----------
    cell : openpyxl.cell.Cell
        Excel cell to process
        
    Returns
    -------
    str
        Extracted PDB identifier in uppercase, or None if cell is empty
    """
    if cell.hyperlink and cell.hyperlink.target:
        return Path(cell.hyperlink.target).stem.upper()
    val = cell.value
    if isinstance(val, str):
        m = HYP_RE.match(val.strip())
        return (m.group("id") if m else val).upper()
    return None if val is None else str(val).strip().upper()


def load_sheet(excel_file: Path, name: str) -> pd.DataFrame:
    """
    Load a specific sheet from an Excel workbook and convert to DataFrame.
    
    Reads an Excel sheet, extracts PDB IDs from hyperlinks, and converts
    the data to a pandas DataFrame with normalized column names.
    
    Parameters
    ----------
    excel_file : Path
        Path to the Excel file
    name : str
        Name of the sheet to load
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the sheet data with normalized column names
        (lowercase, spaces replaced with underscores)
        
    Raises
    ------
    SystemExit
        If the Excel file or sheet does not exist
    """
    if not excel_file.is_file():
        logger.error(f"Missing Excel file: {excel_file}")
        sys.exit(1)
    wb = openpyxl.load_workbook(excel_file, data_only=False)
    if name not in wb.sheetnames:
        logger.error(f"Sheet '{name}' not found (available: {wb.sheetnames})")
        sys.exit(1)
    ws = wb[name]
    headers = [str(c.value).strip().lower().replace(" ", "_") for c in ws[1]]
    rows = []
    for row in ws.iter_rows(min_row=2):
        rec = {h: extract_pdb_from_cell(c) for h, c in zip(headers, row)}
        rows.append(rec)
    df = pd.DataFrame(rows)
    logger.info(f"[{name}] loaded: {df.shape[0]} rows")
    return df


def filter_and_cascade(df_det: pd.DataFrame, df_ali: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter structures and propagate filtering to alignment data.
    
    Ensures consistency between PDB_Details and Alignment_Summary sheets
    by keeping only structures that appear in both datasets.
    
    Parameters
    ----------
    df_det : pd.DataFrame
        DataFrame containing PDB structure details
    df_ali : pd.DataFrame
        DataFrame containing alignment information
        
    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Filtered DataFrames (details, alignment) with matching structure sets
    """
    # Cascade to Alignment_Summary
    valid = set(df_det['pdb_id'].dropna().str.upper())
    df_ali_proc = df_ali[df_ali['pdb_id'].str.upper().isin(valid)].reset_index(drop=True)
    logger.info(f"Processed: PDB_Details {df_det.shape[0]} rows, Alignment_Summary {df_ali_proc.shape[0]} rows")
    return df_det, df_ali_proc


# ---------------------------------------------------------------------------
# v2: Dynamic OPM/PDBTM/RCSB API functions
# ---------------------------------------------------------------------------

def check_opm_availability(pdb_id: str) -> Dict[str, Any]:
    """
    Check if an OPM-oriented structure is available for a PDB ID.
    
    Tries multiple methods:
    1. Direct download URL check
    2. OPM protein page scraping (if BeautifulSoup available)
    
    Parameters
    ----------
    pdb_id : str
        PDB identifier (4 characters)
        
    Returns
    -------
    dict
        {"available": bool, "url": str|None, "source": "opm"|"opm_page"|None}
    """
    pdb_id = pdb_id.lower().strip()
    
    # Try direct OPM download URL
    opm_url = OPM_DOWNLOAD_URL.format(pdb_id=pdb_id)
    try:
        resp = requests.head(opm_url, timeout=HTTP_TIMEOUT)
        if resp.ok:
            logger.debug(f"OPM available for {pdb_id.upper()} via direct URL")
            return {"available": True, "url": opm_url, "source": "opm"}
    except requests.RequestException as e:
        logger.debug(f"OPM head request failed for {pdb_id}: {e}")
    
    # Try OPM protein page (if BeautifulSoup available)
    if HAS_BS4:
        try:
            resp = requests.get(OPM_PROTEIN_PAGE.format(pdb_id=pdb_id), timeout=HTTP_TIMEOUT)
            if resp.ok:
                try:
                    soup = BeautifulSoup(resp.text, 'html.parser')
                    # Look for membrane-related content indicating orientation data
                    if soup.find(text=lambda t: t and 'membrane' in t.lower()):
                        logger.debug(f"OPM available for {pdb_id.upper()} via page scraping")
                        return {"available": True, "url": None, "source": "opm_page"}
                except Exception as e:
                    # BeautifulSoup may fail on Python 3.10+ due to collections.Callable deprecation
                    logger.debug(f"BeautifulSoup parsing failed for {pdb_id}: {e}")
                    # Fall back to simple text search
                    if 'membrane' in resp.text.lower():
                        logger.debug(f"OPM available for {pdb_id.upper()} via text search")
                        return {"available": True, "url": None, "source": "opm_text"}
        except requests.RequestException as e:
            logger.debug(f"OPM page request failed for {pdb_id}: {e}")
    
    return {"available": False, "url": None, "source": None}


def check_pdbtm_availability(pdb_id: str) -> bool:
    """
    Check if PDBTM has transmembrane annotation for a PDB ID.
    
    Parameters
    ----------
    pdb_id : str
        PDB identifier
        
    Returns
    -------
    bool
        True if PDBTM has orientation data for this structure
    """
    if not HAS_BS4:
        logger.debug("BeautifulSoup not available, skipping PDBTM check")
        return False
    
    try:
        resp = requests.get(PDBTM_PAGE.format(pdb_id=pdb_id.lower()), timeout=HTTP_TIMEOUT)
        if resp.ok:
            try:
                soup = BeautifulSoup(resp.text, 'html.parser')
                if soup.find(text=lambda t: t and 'Transmembrane' in t):
                    logger.debug(f"PDBTM available for {pdb_id.upper()}")
                    return True
            except Exception as e:
                # BeautifulSoup may fail on Python 3.10+ due to collections.Callable deprecation
                logger.debug(f"BeautifulSoup parsing failed for PDBTM {pdb_id}: {e}")
                # Fall back to simple text search
                if 'Transmembrane' in resp.text:
                    logger.debug(f"PDBTM available for {pdb_id.upper()} via text search")
                    return True
    except requests.RequestException as e:
        logger.debug(f"PDBTM check failed for {pdb_id}: {e}")
    
    return False


def search_pdb_by_uniprot(uniprot_id: str) -> List[str]:
    """
    Search RCSB PDB for all structures linked to a UniProt accession.
    
    Uses RCSB Search API v2 to find PDB entries.
    
    Parameters
    ----------
    uniprot_id : str
        UniProt accession number
        
    Returns
    -------
    list
        List of PDB IDs (uppercase) linked to this UniProt
    """
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id.upper()
            }
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    
    try:
        resp = requests.post(RCSB_SEARCH_URL, json=query, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        pdb_ids = [r['identifier'].upper() for r in data.get('result_set', [])]
        logger.debug(f"RCSB search for {uniprot_id}: found {len(pdb_ids)} PDBs")
        return sorted(set(pdb_ids))
    except requests.RequestException as e:
        logger.warning(f"RCSB search failed for {uniprot_id}: {e}")
        return []


def classify_adme_family(gene: str) -> str:
    """
    Classify a gene into an ADME transporter family by prefix matching.
    
    Uses gene name prefixes to determine family membership.
    
    Parameters
    ----------
    gene : str
        Gene name (from parse_uniprot_data()['gene'])
        
    Returns
    -------
    str
        ADME family name or "Other" if not matched
    """
    if not gene:
        return "Other"
    
    # Handle composite gene names (e.g., "ABCB1-ATP2A2")
    gene_upper = gene.upper().split("-")[0].strip()
    
    for family, prefixes in ADME_GENE_PREFIXES.items():
        for prefix in prefixes:
            if gene_upper.startswith(prefix.upper()):
                logger.debug(f"Gene {gene} classified as {family}")
                return family
    
    logger.debug(f"Gene {gene} not matched to any ADME family")
    return "Other"


def build_opm_mapping(uniprot_ids: List[str], 
                      pdb_ids_by_uniprot: Dict[str, List[str]],
                      gene_by_uniprot: Dict[str, str]) -> Dict[str, Dict[str, Any]]:
    """
    Build dynamic OPM availability mapping for UniProt IDs.
    
    For each UniProt ID, checks which PDB structures have OPM orientation
    data available and classifies the ADME family.
    
    Parameters
    ----------
    uniprot_ids : list
        List of UniProt accession numbers to process
    pdb_ids_by_uniprot : dict
        PDB IDs already retrieved by parse_uniprot_data() {upid: [pdb_ids]}
    gene_by_uniprot : dict
        Gene names from parse_uniprot_data() {upid: gene}
        
    Returns
    -------
    dict
        Mapping: {uniprot_id: {opm_pdbs, non_opm_pdbs, has_opm, family, endpoints_used}}
    """
    mapping = {}
    
    for upid in uniprot_ids:
        if not upid:
            continue
            
        logger.info(f"Building OPM mapping for {upid}")
        endpoints_used = []
        opm_pdbs = []
        non_opm_pdbs = []
        
        # Get PDB IDs already known from UniProt data (reuse existing data)
        known_pdbs = pdb_ids_by_uniprot.get(upid, [])
        
        # Optionally complement with RCSB Search (may find additional PDBs)
        rcsb_pdbs = search_pdb_by_uniprot(upid)
        if rcsb_pdbs:
            endpoints_used.append("RCSB Search")
        
        # Combine and deduplicate
        all_pdbs = sorted(set([p.upper() for p in known_pdbs] + rcsb_pdbs))
        
        # Check OPM/PDBTM availability for each PDB
        for pdb_id in all_pdbs:
            opm_info = check_opm_availability(pdb_id)
            
            if opm_info["available"]:
                opm_pdbs.append(pdb_id)
                endpoints_used.append(f"OPM ({opm_info['source']})")
            else:
                # Check PDBTM as fallback
                if check_pdbtm_availability(pdb_id):
                    opm_pdbs.append(pdb_id)
                    endpoints_used.append("PDBTM")
                else:
                    non_opm_pdbs.append(pdb_id)
            
            # Rate limiting to avoid hammering APIs
            time.sleep(HTTP_RETRY_DELAY * 0.2)
        
        # Get gene and classify family
        gene = gene_by_uniprot.get(upid, "")
        family = classify_adme_family(gene)
        
        mapping[upid] = {
            "opm_pdbs": opm_pdbs,
            "non_opm_pdbs": non_opm_pdbs,
            "has_opm": len(opm_pdbs) > 0,
            "family": family,
            "gene": gene,
            "template_candidates": [],  # Filled later if needed
            "endpoints_used": sorted(set(endpoints_used))
        }
        
        logger.info(f"  {upid}: {len(opm_pdbs)} OPM, {len(non_opm_pdbs)} non-OPM, family={family}")
    
    return mapping


def download_opm_structure(pdb_id: str, output_path: Path) -> bool:
    """
    Download an OPM-oriented structure.
    
    Parameters
    ----------
    pdb_id : str
        PDB identifier
    output_path : Path
        Where to save the downloaded file
        
    Returns
    -------
    bool
        True if download succeeded
    """
    opm_url = OPM_DOWNLOAD_URL.format(pdb_id=pdb_id.lower())
    
    try:
        resp = requests.get(opm_url, timeout=HTTP_TIMEOUT)
        resp.raise_for_status()
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "wb") as f:
            f.write(resp.content)
        
        logger.info(f"Downloaded OPM structure: {pdb_id} -> {output_path}")
        return True
        
    except requests.RequestException as e:
        logger.warning(f"Failed to download OPM structure {pdb_id}: {e}")
        return False


def find_opm_template_in_group(group_pdb_ids: List[str], 
                                opm_mapping: Dict[str, Dict[str, Any]],
                                upid: str) -> Optional[str]:
    """
    Find an OPM-oriented template within a group of PDB structures.
    
    Searches the group for structures that have OPM orientation data.
    
    Parameters
    ----------
    group_pdb_ids : list
        PDB IDs in the current group
    opm_mapping : dict
        OPM availability mapping from build_opm_mapping()
    upid : str
        UniProt ID for this group
        
    Returns
    -------
    str or None
        First OPM-available PDB ID found, or None if none available
    """
    if upid not in opm_mapping:
        return None
    
    opm_pdbs = set(opm_mapping[upid].get("opm_pdbs", []))
    
    for pdb_id in group_pdb_ids:
        if pdb_id.upper() in opm_pdbs:
            logger.info(f"Found OPM template {pdb_id} for group {upid}")
            return pdb_id.upper()
    
    return None


def superimpose_structures(ref_path: Path,
                           df_targets: pd.DataFrame,
                           out_root: Path) -> list[tuple[str, str, float, float]]:
    """
    Superimpose target structures onto a reference structure using CA atoms.
    
    Performs rigid-body superposition of multiple structures onto a common
    reference, typically one containing DUM orientation markers. Uses CA
    (alpha-carbon) atoms for alignment.
    
    Parameters
    ----------
    ref_path : Path
        Path to the reference structure
    df_targets : pd.DataFrame
        DataFrame containing information about target structures to align
        Must include columns: pdb_id, gene_name, uniprot_id, patched_pdb
    out_root : Path
        Root directory where oriented structures will be saved
        
    Returns
    -------
    list[tuple[str, str, float, float]]
        List of alignment results: (pdb_id, aligned_path, rmsd_before, rmsd_after)
        RMSD values are in Angstroms
    """
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure('ref', str(ref_path))
    ref_atoms = [a for a in ref_struct[0].get_atoms() if a.get_id() == 'CA']
    sup = Superimposer()
    io = PDBIO()

    results = []
    for idx, row in df_targets.iterrows():
        pdb_id = row['pdb_id']
        gene = row['gene_name']
        uniprot = row['uniprot_id']
        in_path = Path(row['patched_pdb'])
        # Prepare output path
        subdir = out_root / f"{gene}_{uniprot}"
        subdir.mkdir(parents=True, exist_ok=True)
        out_name = f"oriented_{gene}_{uniprot}_{pdb_id}.pdb"
        out_path = subdir / out_name
        # Load and align
        tgt_struct = parser.get_structure(pdb_id, str(in_path))
        tgt_atoms = [a for a in tgt_struct[0].get_atoms() if a.get_id() == 'CA']
        n = min(len(ref_atoms), len(tgt_atoms))
        ref_coords = np.array([a.get_coord() for a in ref_atoms[:n]])
        tgt_coords = np.array([a.get_coord() for a in tgt_atoms[:n]])
        # RMSD before alignment
        rms_before = np.sqrt(((ref_coords - tgt_coords)**2).sum(axis=1).mean())
        # Superimpose and get RMSD after
        sup.set_atoms(ref_atoms[:n], tgt_atoms[:n])
        rms_after = sup.rms
        sup.apply(tgt_struct.get_atoms())
        io.set_structure(tgt_struct)
        io.save(str(out_path))
        logger.info(f"{pdb_id}: before={rms_before:.3f}Å, after={rms_after:.3f}Å -> {out_path}")
        results.append((pdb_id, str(out_path), rms_before, rms_after))
    return results


# ---------------------------------------------------------------------------
# Functions for multi-reference alignment workflow
# ---------------------------------------------------------------------------

def group_structures_by_attribute(df_det: pd.DataFrame, attribute: str) -> Dict[str, List[str]]:
    """
    Group PDB structures by a specified attribute for batch processing.
    
    Creates groups of related structures (e.g., all structures for the same
    gene or protein) to enable group-wise orientation to a common reference.
    
    Parameters
    ----------
    df_det : pd.DataFrame
        DataFrame containing PDB structure details
    attribute : str
        Column name to group by ('gene_name' or 'uniprot_id')
        
    Returns
    -------
    Dict[str, List[str]]
        Dictionary mapping attribute values to lists of PDB IDs
        Example: {'ABCB1': ['6C0V', '6QEX', '7A69'], 'ABCB4': ['3G9P']}
    """
    if attribute not in df_det.columns:
        logger.error(f"Attribute '{attribute}' not found in DataFrame")
        return {}
    
    groups = {}
    for _, row in df_det.iterrows():
        attr_val = row[attribute]
        pdb_id = row['pdb_id']
        if attr_val and pdb_id:
            if attr_val not in groups:
                groups[attr_val] = []
            groups[attr_val].append(pdb_id)
    
    return groups


def find_dum_references(df_det: pd.DataFrame, groups: Dict[str, List[str]]) -> Dict[str, str]:
    """
    Identify DUM-containing structures to use as orientation references.
    
    For each group of structures, searches for structures containing DUM
    (dummy atom) markers. DUM atoms are artificial markers indicating membrane
    orientation in membrane protein structures. The first DUM-containing
    structure in each group is designated as the orientation reference.
    
    Parameters
    ----------
    df_det : pd.DataFrame
        DataFrame containing PDB structure details with 'small_molecules' column
    groups : Dict[str, List[str]]
        Dictionary of structure groups from group_structures_by_attribute
        
    Returns
    -------
    Dict[str, str]
        Dictionary mapping group names to reference PDB IDs
        Only groups with DUM-containing structures are included
    """
    references = {}
    
    # Check if 'small_molecules' column exists
    if 'small_molecules' not in df_det.columns:
        logger.warning("Column 'small_molecules' not found, cannot detect DUM")
        return references
    
    # For each group, find a reference with DUM
    for group, pdb_ids in groups.items():
        # Filter DataFrame for this group
        group_df = df_det[df_det['pdb_id'].isin(pdb_ids)].copy()
        
        # Search for structures with DUM
        sm = group_df['small_molecules'].fillna("").str.upper()
        dum_mask = sm.str.contains(DUM_RE)
        
        if dum_mask.any():
            # Take the first structure with DUM as reference
            ref_id = group_df.loc[dum_mask, 'pdb_id'].iloc[0]
            references[group] = ref_id
            logger.info(f"Group {group}: DUM reference found: {ref_id}")
        else:
            logger.warning(f"Group {group}: no DUM reference found")
    
    return references


def align_structures_by_groups(excel_file: Path, 
                              group_by: str = 'gene_name',
                              out_root: Optional[Path] = None) -> Dict[str, List[tuple]]:
    """
    Perform multi-reference structural alignment grouped by gene or protein.
    
    This is the main orchestrator function that:
    1. Loads structure and alignment data from Excel
    2. Groups structures by the specified attribute
    3. Identifies DUM-containing references for each group
    4. Aligns all structures in each group to their reference
    5. Saves oriented structures and RMSD metrics
    
    This approach ensures that all structures of the same protein share a
    common spatial orientation, which is critical for:
    - Comparative structural analysis
    - Visualization of conformational changes
    - Identification of conserved binding sites
    - Analysis of membrane protein orientation
    
    Parameters
    ----------
    excel_file : Path
        Path to the Excel file containing PDB_Details and Alignment_Summary sheets
    group_by : str, optional
        Attribute for grouping: 'gene_name' (default) or 'uniprot_id'
    out_root : Path, optional
        Root directory for oriented structures (default: 'pdb_structures_oriented')
        
    Returns
    -------
    Dict[str, List[tuple]]
        Results dictionary: {group: [(pdb_id, aligned_path, rmsd_before, rmsd_after)]}
        
    Notes
    -----
    - Structures without a DUM reference are not oriented
    - The reference structure itself is copied (not re-aligned)
    - Results are saved to a new 'Alignment_RMSD_MultiRef' sheet in the Excel file
    """
    if out_root is None:
        out_root = Path('pdb_structures_oriented')
    
    # Load Excel sheets
    df_det_all = load_sheet(excel_file, 'PDB_Details')
    df_ali_all = load_sheet(excel_file, 'Alignment_Summary')
    
    # Check that required columns exist
    required_cols = ['pdb_id', group_by, 'small_molecules']
    for col in required_cols:
        if col not in df_det_all.columns:
            logger.error(f"Required column '{col}' not found in PDB_Details")
            return {}
    
    if 'patched_pdb' not in df_ali_all.columns:
        logger.error("Column 'patched_pdb' not found in Alignment_Summary")
        return {}
    
    # Group structures by specified attribute
    groups = group_structures_by_attribute(df_det_all, group_by)
    logger.info(f"Structures grouped by {group_by}: {len(groups)} groups found")
    
    # Find DUM references for each group
    references = find_dum_references(df_det_all, groups)
    logger.info(f"DUM references found for {len(references)}/{len(groups)} groups")
    
    # Prepare results dictionary
    all_results = {}
    all_aligned_results = []
    
    # Process each group that has a DUM-containing reference
    for group, ref_id in references.items():
        logger.info(f"Processing group {group} with reference {ref_id}")
        
        # Locate the reference structure path
        ref_row = df_ali_all[df_ali_all['pdb_id'] == ref_id]
        if ref_row.empty:
            logger.warning(f"Reference {ref_id} not found in Alignment_Summary, skipping group {group}")
            continue
        
        ref_path = Path(ref_row.iloc[0]['patched_pdb'])
        if not ref_path.is_file():
            logger.error(f"Reference file not found: {ref_path}, skipping group {group}")
            continue
        
        # Identify target structures (all structures except the reference)
        group_pdb_ids = groups[group]
        target_pdb_ids = [pid for pid in group_pdb_ids if pid != ref_id]
        
        # Prepare DataFrame of target structures with necessary metadata
        df_targets = pd.merge(
            df_ali_all[df_ali_all['pdb_id'].isin(target_pdb_ids)],
            df_det_all[['pdb_id', 'gene_name', 'uniprot_id']],
            on='pdb_id', how='left'
        )
        
        if df_targets.empty:
            logger.warning(f"No target structures found for group {group}")
            continue
        
        # Perform structural superposition for all targets
        results = superimpose_structures(ref_path, df_targets, out_root)
        
        # Copy the reference structure to the output directory (already oriented)
        ref_det_row = df_det_all[df_det_all['pdb_id'] == ref_id].iloc[0]
        gene = ref_det_row['gene_name']
        uniprot = ref_det_row['uniprot_id']
        subdir = out_root / f"{gene}_{uniprot}"
        subdir.mkdir(parents=True, exist_ok=True)
        out_name = f"oriented_{gene}_{uniprot}_{ref_id}.pdb"
        out_path = subdir / out_name
        shutil.copy(str(ref_path), str(out_path))
        logger.info(f"Copied reference {ref_id} -> {out_path}")
        
        # Add reference to results (no RMSD since it wasn't aligned)
        results.append((ref_id, str(out_path), None, None))
        
        # Store results for this group and accumulate for Excel output
        all_results[group] = results
        all_aligned_results.extend(results)
    
    # Save consolidated results to Excel
    df_res = pd.DataFrame(all_aligned_results, columns=['pdb_id', 'aligned_path', 'rmsd_before', 'rmsd_after'])
    with pd.ExcelWriter(excel_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        df_res.to_excel(writer, sheet_name='Alignment_RMSD_MultiRef', index=False)
    logger.info(f"Results saved in 'Alignment_RMSD_MultiRef' in {excel_file}")
    
    return all_results


def align_structures_by_groups_v2(
    excel_file: Path,
    opm_mapping: Dict[str, Dict[str, Any]],
    conformation_by_pdb: Dict[str, str],
    group_by: str = 'gene_name',
    out_root: Optional[Path] = None
) -> Dict[str, List[Dict[str, Any]]]:
    """
    v2: Multi-reference structural alignment with dynamic OPM integration.
    
    Enhanced version that:
    1. Uses dynamic OPM availability mapping
    2. Downloads OPM templates on-demand if not locally available
    3. Organizes output by Family/Gene/Conformation/PDB_ID
    4. Saves metadata.json with orientation info for each structure
    
    Parameters
    ----------
    excel_file : Path
        Path to the Excel file containing PDB_Details and Alignment_Summary sheets
    opm_mapping : dict
        OPM availability mapping from build_opm_mapping()
    conformation_by_pdb : dict
        Conformation state per PDB ID from get_pdb_details()
    group_by : str, optional
        Attribute for grouping: 'gene_name' (default) or 'uniprot_id'
    out_root : Path, optional
        Root directory for oriented structures (default: 'pdb_structures_oriented_v2')
        
    Returns
    -------
    dict
        Results: {group: [result_dict for each structure]}
    """
    if out_root is None:
        out_root = Path('pdb_structures_oriented_v2')
    
    # Load Excel sheets
    df_det_all = load_sheet(excel_file, 'PDB_Details')
    df_ali_all = load_sheet(excel_file, 'Alignment_Summary')
    
    # Check required columns
    required_cols = ['pdb_id', group_by]
    for col in required_cols:
        if col not in df_det_all.columns:
            logger.error(f"Required column '{col}' not found in PDB_Details")
            return {}
    
    if 'patched_pdb' not in df_ali_all.columns:
        logger.error("Column 'patched_pdb' not found in Alignment_Summary")
        return {}
    
    # Group structures by specified attribute
    groups = group_structures_by_attribute(df_det_all, group_by)
    logger.info(f"v2: Structures grouped by {group_by}: {len(groups)} groups found")
    
    # Prepare results
    all_results = {}
    all_result_records = []
    
    # Process each group
    for group, pdb_ids in groups.items():
        logger.info(f"v2: Processing group {group} ({len(pdb_ids)} structures)")
        
        # Get UniProt ID for this group
        group_rows = df_det_all[df_det_all[group_by] == group]
        if group_rows.empty:
            continue
        
        upid_col = 'uniprot_id' if 'uniprot_id' in group_rows.columns else None
        upid = group_rows.iloc[0][upid_col] if upid_col else ""
        
        # Get family from OPM mapping (or classify from gene)
        if upid and upid in opm_mapping:
            family = opm_mapping[upid].get("family", "Other")
        else:
            family = classify_adme_family(group)
        
        # Find OPM template in this group
        ref_id = None
        ref_path = None
        orientation_method = "none"
        
        # First, check local DUM references (already downloaded from OPM)
        dum_refs = find_dum_references(df_det_all, {group: pdb_ids})
        if group in dum_refs:
            ref_id = dum_refs[group]
            ref_row = df_ali_all[df_ali_all['pdb_id'] == ref_id]
            if not ref_row.empty:
                ref_path = Path(ref_row.iloc[0]['patched_pdb'])
                if ref_path.is_file():
                    orientation_method = "local_dum"
                    logger.info(f"  Using local DUM reference: {ref_id}")
        
        # If no local DUM, check OPM mapping for available templates
        if ref_path is None and upid and upid in opm_mapping:
            opm_template = find_opm_template_in_group(pdb_ids, opm_mapping, upid)
            if opm_template:
                # Download OPM structure as reference
                opm_ref_path = out_root / "_opm_templates" / f"{opm_template}.pdb"
                if download_opm_structure(opm_template, opm_ref_path):
                    ref_id = opm_template
                    ref_path = opm_ref_path
                    orientation_method = "opm_download"
                    logger.info(f"  Downloaded OPM template: {ref_id}")
        
        # Process each structure in the group
        group_results = []
        parser = PDBParser(QUIET=True)
        io = PDBIO()
        sup = Superimposer()
        
        for pdb_id in pdb_ids:
            # Get structure info
            det_row = df_det_all[df_det_all['pdb_id'] == pdb_id]
            ali_row = df_ali_all[df_ali_all['pdb_id'] == pdb_id]
            
            if det_row.empty or ali_row.empty:
                continue
            
            gene = det_row.iloc[0].get('gene_name', group)
            upid_struct = det_row.iloc[0].get('uniprot_id', upid)
            patched_path = ali_row.iloc[0].get('patched_pdb')
            aligned_path = ali_row.iloc[0].get('aligned_pdb', None)
            
            if not patched_path or not Path(patched_path).is_file():
                continue
            
            # Get conformation
            conformation = conformation_by_pdb.get(pdb_id, "unspecified")
            if not conformation or conformation == "Not specified":
                conformation = "unspecified"
            conformation = conformation.replace(" ", "_").replace(",", "_").lower()
            
            # Create output directory: Family/Gene/Conformation/PDB_ID
            struct_dir = out_root / family / gene / conformation / pdb_id
            struct_dir.mkdir(parents=True, exist_ok=True)
            
            oriented_path = struct_dir / f"oriented_{pdb_id}.pdb"
            
            # Perform orientation
            result = {
                "pdb_id": pdb_id,
                "gene": gene,
                "upid": upid_struct,
                "family": family,
                "conformation": conformation,
                "method": orientation_method,
                "template_used": ref_id,
                "rmsd_before": None,
                "rmsd_after": None,
                "status": "ok",
                "notes": []
            }
            
            if ref_path and ref_path.is_file() and pdb_id != ref_id:
                # Align to reference
                try:
                    ref_struct = parser.get_structure('ref', str(ref_path))
                    tgt_struct = parser.get_structure(pdb_id, str(patched_path))
                    
                    ref_atoms = [a for a in ref_struct[0].get_atoms() if a.get_id() == 'CA']
                    tgt_atoms = [a for a in tgt_struct[0].get_atoms() if a.get_id() == 'CA']
                    
                    n = min(len(ref_atoms), len(tgt_atoms))
                    if n >= 10:
                        ref_coords = np.array([a.get_coord() for a in ref_atoms[:n]])
                        tgt_coords = np.array([a.get_coord() for a in tgt_atoms[:n]])
                        
                        rms_before = float(np.sqrt(((ref_coords - tgt_coords)**2).sum(axis=1).mean()))
                        
                        sup.set_atoms(ref_atoms[:n], tgt_atoms[:n])
                        rms_after = float(sup.rms)
                        sup.apply(tgt_struct.get_atoms())
                        
                        io.set_structure(tgt_struct)
                        io.save(str(oriented_path))
                        
                        result["rmsd_before"] = round(rms_before, 3)
                        result["rmsd_after"] = round(rms_after, 3)
                        result["notes"].append(f"Aligned to {ref_id}")
                        
                        logger.info(f"  {pdb_id}: RMSD {rms_before:.3f}Å -> {rms_after:.3f}Å")
                    else:
                        # Not enough atoms, copy without alignment
                        shutil.copy(str(patched_path), str(oriented_path))
                        result["status"] = "partial"
                        result["notes"].append("Insufficient CA atoms for alignment")
                        
                except Exception as e:
                    logger.warning(f"  {pdb_id}: alignment failed: {e}")
                    shutil.copy(str(patched_path), str(oriented_path))
                    result["status"] = "error"
                    result["notes"].append(f"Alignment error: {str(e)}")
                    
            elif pdb_id == ref_id:
                # This is the reference, copy it
                if ref_path and ref_path.is_file():
                    shutil.copy(str(ref_path), str(oriented_path))
                else:
                    shutil.copy(str(patched_path), str(oriented_path))
                result["notes"].append("Reference structure (already oriented)")
                
            else:
                # No reference available, copy without orientation
                shutil.copy(str(patched_path), str(oriented_path))
                result["method"] = "none"
                result["status"] = "needs_manual_review"
                result["notes"].append("No OPM template available")
            
            # Copy other versions
            shutil.copy(str(patched_path), struct_dir / f"patched_{pdb_id}.pdb")
            if aligned_path and Path(aligned_path).is_file():
                shutil.copy(str(aligned_path), struct_dir / f"aligned_{pdb_id}.pdb")
            
            # Save metadata.json
            result["oriented_path"] = str(oriented_path)
            metadata = {k: v for k, v in result.items()}
            with open(struct_dir / "metadata.json", "w", encoding="utf-8") as f:
                json.dump(metadata, f, indent=2, ensure_ascii=False)
            
            group_results.append(result)
            all_result_records.append(result)
        
        all_results[group] = group_results
    
    # Save to Excel
    if all_result_records:
        df_v2 = pd.DataFrame(all_result_records)
        with pd.ExcelWriter(excel_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            df_v2.to_excel(writer, sheet_name='Orientation_v2', index=False)
        logger.info(f"v2 results saved in 'Orientation_v2' sheet")
    
    return all_results


# ---------------------------------------------------------------------------
# Main function for command-line usage
# ---------------------------------------------------------------------------

def main():
    """
    Main function for command-line usage.
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Align and orient PDB files by groups on references containing DUM"
    )
    parser.add_argument('-e', '--excel', required=True,
                        help='Excel file containing the data')
    parser.add_argument('-g', '--group-by', choices=['gene_name', 'uniprot_id'], default='gene_name',
                        help='Attribute for grouping (gene_name or uniprot_id)')
    parser.add_argument('-o', '--output', default='pdb_structures_oriented',
                        help='Root directory for aligned structures')
    
    args = parser.parse_args()
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    
    # Execute alignment
    excel_file = Path(args.excel)
    out_root = Path(args.output)
    
    align_structures_by_groups(excel_file, args.group_by, out_root)


if __name__ == '__main__':
    main()
