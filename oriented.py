#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Multi-reference alignment module for PDB structures
# ---------------------------------------------------------------------------
"""
oriented.py

Multi-reference structural orientation using DUM markers.

This module performs systematic orientation of protein structures using
DUM (dummy atom) ligands as spatial references. DUM atoms are artificial
markers placed in membrane protein structures to indicate membrane orientation.

Key features:
- Groups structures by gene name or UniProt ID
- Identifies DUM-containing structures as orientation references
- Aligns all structures in a group to their common DUM reference
- Handles cases where no DUM reference exists
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
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Any

import numpy as np
import pandas as pd
import openpyxl
from Bio.PDB import PDBParser, Superimposer, PDBIO

# --- Constants & regular expressions ---
# Regular expression to detect DUM (dummy atom) markers in small molecule lists
DUM_RE = re.compile(r"\bDUM\b", re.IGNORECASE)
# Regular expression to extract PDB IDs from Excel HYPERLINK formulas
# Regular expression to extract PDB IDs from Excel HYPERLINK formulas
HYP_RE = re.compile(
    r'=?HYPERLINK\(\s*"[^"]*"\s*[;,]\s*"(?P<id>[A-Za-z0-9\-]+)"\s*\)',
    re.IGNORECASE
)

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
