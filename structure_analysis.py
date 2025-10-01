"""
structure_analysis.py

Parse PDB files to extract structural features and quality metrics.

This module provides utilities for analyzing protein structure files:
 - Secondary structure annotations (HELIX/SHEET records)
 - SEQRES declared sequence length
 - Missing-residue gaps and low-confidence coil regions
 - Small-molecule ligands (HETATM entries)
 
These analyses help assess structure quality and identify regions that
may need patching or special attention during structural alignment.
"""

import re
from typing import List, Tuple

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa


def parse_helix_sheet_annotations(pdb_file: str, chain_id: str = 'A') -> dict:
    """
    Parse HELIX and SHEET records from a PDB file for secondary structure annotation.
    
    Reads the PDB header to identify helical and sheet regions for the specified chain,
    creating a mapping of residue positions to their secondary structure type.
    
    Parameters
    ----------
    pdb_file : str
        Path to the PDB file to parse
    chain_id : str, optional
        Chain identifier to extract annotations for, default 'A'
        
    Returns
    -------
    dict
        Dictionary mapping residue numbers to secondary structure type:
        'H' for helix, 'E' for beta-sheet (extended strand)
    """
    ss_map = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("HELIX") and len(line) >= 38 and line[19] == chain_id:
                start = int(line[21:25].strip())
                end = int(line[33:37].strip())
                for res_id in range(start, end + 1):
                    ss_map[res_id] = 'H'
            elif line.startswith("SHEET") and len(line) >= 37 and line[21] == chain_id:
                start = int(line[22:26].strip())
                end = int(line[33:37].strip())
                for res_id in range(start, end + 1):
                    ss_map[res_id] = 'E'
    return ss_map


def get_seqres_length(pdb_file: str, chain_id: str = 'A') -> int:
    """
    Get the declared sequence length from SEQRES records in a PDB file.
    
    SEQRES records contain the complete amino acid sequence as defined in the
    source database, which may differ from the observed residues in the ATOM records.
    
    Parameters
    ----------
    pdb_file : str
        Path to the PDB file to parse
    chain_id : str, optional
        Chain identifier to extract length for, default 'A'
        
    Returns
    -------
    int
        Total number of residues declared in SEQRES for the specified chain
    """
    length = 0
    with open(pdb_file, 'r') as f:
        for ln in f:
            if ln.startswith("SEQRES") and ln[11] == chain_id:
                # count residue names on that line
                length += len(ln[19:].split())
    return length


def group_segments(numbers: List[int]) -> List[Tuple[int, int]]:
    """
    Group a sorted list of integers into contiguous segments.
    
    Identifies consecutive integers and groups them into (start, end) pairs,
    useful for representing residue ranges or gap regions compactly.
    
    Parameters
    ----------
    numbers : List[int]
        Sorted list of integers to group
        
    Returns
    -------
    List[Tuple[int, int]]
        List of (start, end) tuples representing contiguous segments
        
    Examples
    --------
    >>> group_segments([1, 2, 3, 5, 6, 10])
    [(1, 3), (5, 6), (10, 10)]
    """
    from itertools import groupby, count
    segments = []
    for _, grp in groupby(numbers, key=lambda n, c=count(): n - next(c)):
        grp = list(grp)
        segments.append((grp[0], grp[-1]))
    return segments


def detect_gaps_and_discard_coils(
    pdb_file: str,
    chain_id: str,
    plddt_threshold: float,
    coil_length_threshold: int
) -> List[Tuple[str, int, int]]:
    """
    Detect structural issues including missing residues and low-confidence regions.
    
    Identifies two types of structural problems:
    1. Missing residues (gaps) detected from REMARK 465 records or residue numbering
    2. Low-confidence coil segments where average B-factor is below threshold
    
    Gaps are classified as N-terminal, C-terminal, or internal based on their position.
    Low-confidence regions are only reported if they exceed the minimum length threshold.
    
    Parameters
    ----------
    pdb_file : str
        Path to the PDB file to analyze
    chain_id : str
        Chain identifier to analyze
    plddt_threshold : float
        B-factor/pLDDT threshold below which residues are considered low-confidence
    coil_length_threshold : int
        Minimum number of consecutive low-confidence residues to report
        
    Returns
    -------
    List[Tuple[str, int, int]]
        List of issues as (issue_type, start_residue, end_residue) tuples
        Issue types: "GAP_NTER" (N-terminal gap), "GAP" (internal gap),
                    "GAP_CTER" (C-terminal gap), "Coil_rejeté" (rejected coil)
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("struct", pdb_file)
    except Exception:
        return []

    # find chain
    chain = next((c for c in structure[0] if c.id == chain_id), None)
    if chain is None:
        return []

    present_ids = []
    seq_plddt = []

    # record which residues are present and their avg B-factor
    for res in chain:
        hetflag, resnum, icode = res.id
        if hetflag == ' ' and is_aa(res):
            present_ids.append(resnum)
            bvals = [a.bfactor for a in res if a.bfactor is not None]
            avg = sum(bvals) / len(bvals) if bvals else 0.0
            seq_plddt.append((resnum, avg))

    if not present_ids:
        return []

    min_pres, max_pres = min(present_ids), max(present_ids)
    missing = []

    # parse REMARK 465 lines
    pattern = re.compile(rf"^REMARK 465\s+\w{{3}}\s+{chain_id}\s+(\d+)")
    with open(pdb_file, 'r') as fh:
        for ln in fh:
            m = pattern.match(ln)
            if m:
                missing.append(int(m.group(1)))

    issues = []
    if missing:
        for start, end in group_segments(sorted(missing)):
            if end < min_pres:
                t = "GAP_NTER"
            elif start > max_pres:
                t = "GAP_CTER"
            else:
                t = "GAP"
            issues.append((t, start, end))
    else:
        # infer gaps by numbering discontinuities
        sorted_ids = sorted(present_ids)
        for a, b in zip(sorted_ids, sorted_ids[1:]):
            if b != a + 1:
                issues.append(("GAP", a + 1, b - 1))
        # N‐terminal
        if min_pres > 1:
            issues.append(("GAP_NTER", 1, min_pres - 1))
        # C‐terminal
        seq_len = get_seqres_length(pdb_file, chain_id) or max_pres
        if max_pres < seq_len:
            issues.append(("GAP_CTER", max_pres + 1, seq_len))

    # now detect coil (low‐confidence) segments
    coil_regions = []
    current = []
    for resnum, avg in seq_plddt:
        if avg < plddt_threshold:
            current.append(resnum)
        else:
            if len(current) >= coil_length_threshold:
                coil_regions.append(tuple(current))
            current = []
    if len(current) >= coil_length_threshold:
        coil_regions.append(tuple(current))

    for region in coil_regions:
        issues.append(("Coil_rejeté", region[0], region[-1]))

    return issues


def parse_small_molecules(pdb_file: str) -> List[str]:
    """
    Extract small molecule ligand names from a PDB file.
    
    Parses both HET header records and HETATM coordinate lines to identify
    all small molecules (ligands, cofactors, substrates) present in the structure.
    Water molecules (HOH, WAT) are automatically excluded from the results.
    
    Parameters
    ----------
    pdb_file : str
        Path to the PDB file to parse
        
    Returns
    -------
    List[str]
        Sorted list of unique small molecule residue names (3-letter codes)
    """
    ligands = set()
    with open(pdb_file, 'r') as f:
        for ln in f:
            # header HET records
            if ln.startswith("HET   "):
                resn = ln[7:10].strip()
                if resn and resn not in ("HOH", "WAT"):
                    ligands.add(resn)
            # actual coordinates
            elif ln.startswith("HETATM"):
                resn = ln[17:20].strip()
                if resn and resn not in ("HOH", "WAT"):
                    ligands.add(resn)
    return sorted(ligands)

