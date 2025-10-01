#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Multi-reference alignment module for PDB structures
# ---------------------------------------------------------------------------

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
DUM_RE = re.compile(r"\bDUM\b", re.IGNORECASE)
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
    Extract PDB identifier from an Excel cell, whether as hyperlink or text.
    
    Parameters
    ----------
    cell : openpyxl.cell.Cell
        Excel cell to process
        
    Returns
    -------
    str
        Extracted PDB identifier, or None if the cell is empty
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
    Load a specific sheet from an Excel file and convert it to DataFrame.
    
    Parameters
    ----------
    excel_file : Path
        Path to the Excel file
    name : str
        Name of the sheet to load
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the sheet data
        
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
    Filter structures and cascade filtering to alignment data.
    
    Parameters
    ----------
    df_det : pd.DataFrame
        DataFrame containing PDB structure details
    df_ali : pd.DataFrame
        DataFrame containing alignment information
        
    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Tuple containing filtered DataFrames (details, alignment)
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
    Superimpose target structures onto a reference structure.
    
    Parameters
    ----------
    ref_path : Path
        Path to the reference structure
    df_targets : pd.DataFrame
        DataFrame containing information about target structures
    out_root : Path
        Root directory for aligned structures
        
    Returns
    -------
    list[tuple[str, str, float, float]]
        List of tuples (pdb_id, aligned_path, rmsd_before, rmsd_after)
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
# New functions for multi-reference alignment
# ---------------------------------------------------------------------------

def group_structures_by_attribute(df_det: pd.DataFrame, attribute: str) -> Dict[str, List[str]]:
    """
    Group PDB structures by specified attribute (gene_name or uniprot_id).
    
    Parameters
    ----------
    df_det : pd.DataFrame
        DataFrame containing PDB structure details
    attribute : str
        Attribute for grouping ('gene_name' or 'uniprot_id')
        
    Returns
    -------
    Dict[str, List[str]]
        Dictionary {attribute_value: [list_of_pdb_ids]}
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
    Find a reference containing DUM for each group.
    
    Parameters
    ----------
    df_det : pd.DataFrame
        DataFrame containing PDB structure details
    groups : Dict[str, List[str]]
        Dictionary of structure groups
        
    Returns
    -------
    Dict[str, str]
        Dictionary {group: reference_pdb_id}
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
    Align PDB structures by groups, using a DUM reference for each group.
    
    Parameters
    ----------
    excel_file : Path
        Path to the Excel file containing the data
    group_by : str, optional
        Attribute for grouping ('gene_name' or 'uniprot_id'), default 'gene_name'
    out_root : Path, optional
        Root directory for aligned structures, default 'pdb_structures_oriented'
        
    Returns
    -------
    Dict[str, List[tuple]]
        Dictionary {group: [(pdb_id, aligned_path, rmsd_before, rmsd_after)]}
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
    
    # For each group with a reference
    for group, ref_id in references.items():
        logger.info(f"Processing group {group} with reference {ref_id}")
        
        # Get the reference path
        ref_row = df_ali_all[df_ali_all['pdb_id'] == ref_id]
        if ref_row.empty:
            logger.warning(f"Reference {ref_id} not found in Alignment_Summary, skipping group {group}")
            continue
        
        ref_path = Path(ref_row.iloc[0]['patched_pdb'])
        if not ref_path.is_file():
            logger.error(f"Reference file not found: {ref_path}, skipping group {group}")
            continue
        
        # Filter group structures (except reference)
        group_pdb_ids = groups[group]
        target_pdb_ids = [pid for pid in group_pdb_ids if pid != ref_id]
        
        # Prepare target DataFrame
        df_targets = pd.merge(
            df_ali_all[df_ali_all['pdb_id'].isin(target_pdb_ids)],
            df_det_all[['pdb_id', 'gene_name', 'uniprot_id']],
            on='pdb_id', how='left'
        )
        
        if df_targets.empty:
            logger.warning(f"No target structures found for group {group}")
            continue
        
        # Align structures to reference
        results = superimpose_structures(ref_path, df_targets, out_root)
        
        # Copy reference (already oriented)
        ref_det_row = df_det_all[df_det_all['pdb_id'] == ref_id].iloc[0]
        gene = ref_det_row['gene_name']
        uniprot = ref_det_row['uniprot_id']
        subdir = out_root / f"{gene}_{uniprot}"
        subdir.mkdir(parents=True, exist_ok=True)
        out_name = f"oriented_{gene}_{uniprot}_{ref_id}.pdb"
        out_path = subdir / out_name
        shutil.copy(str(ref_path), str(out_path))
        logger.info(f"Copied reference {ref_id} -> {out_path}")
        
        # Add reference to results
        results.append((ref_id, str(out_path), None, None))
        
        # Store results for this group
        all_results[group] = results
        all_aligned_results.extend(results)
    
    # Save all results to Excel
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
