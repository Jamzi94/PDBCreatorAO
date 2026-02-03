"""
alignment.py

Structural alignment and gap patching for protein structures.

This module provides comprehensive functionality for aligning predicted protein
structures (e.g., from AlphaFold) onto experimental structures and patching
missing regions (gaps). It handles:

- Multi-chain complexes (both homomeric and heteromeric)
- Optimal chain-to-chain correspondence using sequence and structural alignment
- Gap detection and filling using predicted structure coordinates
- Fusion protein detection and handling
- RMSD calculation before and after alignment
- Temporary file management with automatic cleanup

The alignment process uses both BioPython for sequence alignment and MDAnalysis
for structural superposition to find the best chain correspondence.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Imports & logging
# ---------------------------------------------------------------------------
import glob
import logging
import os
import tempfile
import shutil
from openpyxl.styles import Font
from typing import Dict, List, Tuple, Optional, Any, Union
import warnings
import MDAnalysis as mda
import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Select, Superimposer
from Bio.PDB.Polypeptide import PPBuilder
from Bio import pairwise2
from MDAnalysis.analysis import align, rms
from contextlib import contextmanager
import re
from pathlib import Path

from structure_analysis import group_segments 


# Suppress specific warnings commonly encountered with MDAnalysis
# These warnings are expected and don't indicate problems with the analysis
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", message="Reader has no dt information")
warnings.filterwarnings("ignore", message="Unit cell dimensions not found*")
warnings.filterwarnings("ignore", message="Found no information for attr: 'formalcharges'")
warnings.filterwarnings("ignore", message="Unknown element .*")

# Reduce MDAnalysis logger level if needed
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)

log = logging.getLogger(__name__)

# Regular expression pattern for detecting DUM (dummy atom) ligands
# DUM atoms are used as orientation references in membrane protein structures
DUM_RE = re.compile(r"\bDUM\b", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Helper: highest residue number containing a Cα atom for a given chain
# ---------------------------------------------------------------------------

def _get_max_resid(src: str, chain_id: str, parser: Optional[PDBParser] = None) -> int:
    """
    Get the maximum residue number containing a Cα atom for a given chain.
    
    Parameters
    ----------
    src : str
        Path to the PDB file
    chain_id : str
        Chain identifier
    parser : PDBParser, optional
        PDB parser to reuse, default None (creates a new parser)
        
    Returns
    -------
    int
        Maximum residue number
        
    Raises
    ------
    RuntimeError
        If the chain is not found
    """
    if parser is None:
        parser = PDBParser(QUIET=True)
    struct = parser.get_structure("tmp", src)
    for ch in struct.get_chains():
        if ch.id == chain_id:
            return max(res.id[1] for res in ch if res.id[0] == " " and "CA" in res)
    raise RuntimeError(f"Chain {chain_id} not found in {src!r}")

# ---------------------------------------------------------------------------
# Temporary file management
# ---------------------------------------------------------------------------

def create_persistent_temp_file(suffix=".pdb"):
    """
    Create a persistent temporary file that won't be automatically deleted.
    
    Creates a temporary file that must be explicitly removed. This is necessary
    because some operations require the file to persist beyond the scope of
    the function that creates it.
    
    Parameters
    ----------
    suffix : str, optional
        File extension for the temporary file, default ".pdb"
        
    Returns
    -------
    str
        Absolute path to the created temporary file
    """
    fd, tmp_path = tempfile.mkstemp(suffix=suffix)
    os.close(fd)
    return tmp_path

def cleanup_temp_files(file_list):
    """
    Clean up a list of temporary files.
    
    Safely removes temporary files, logging warnings if deletion fails.
    This ensures proper cleanup even when the pipeline encounters errors.
    
    Parameters
    ----------
    file_list : list
        List of file paths to remove
    """
    for file_path in file_list:
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
            except Exception as e:
                log.warning(f"Failed to remove temporary file {file_path}: {e}")

# ---------------------------------------------------------------------------
# Chain-pairing utilities
# ---------------------------------------------------------------------------

def _find_best_chain_pair(pred_dir: str, exp_path: str) -> Dict[str, Tuple[str, Dict]]:
    """
    Pair predicted chains with experimental chains using optimal RMSD alignment.

    Uses both sequence alignment and structural superposition to find the best
    matching between predicted and experimental chains. This handles complex cases:
    - Multi-chain structures
    - Different chain orderings
    - Partial overlaps
    
    The function uses a two-pass approach:
    1. Sequence alignment to find similar chains
    2. Structural alignment (RMSD) to confirm and rank matches
    
    Parameters
    ----------
    pred_dir : str
        Directory containing predicted PDB files
    exp_path : str
        Path to experimental PDB file
        
    Returns
    -------
    mapping : dict
        Dictionary mapping "<pred_file>:<pred_chain>" to (best_exp_chain, details)
        where details contains:
        - rmsd: Best RMSD value
        - pairs: List of all tested chain pairs
        - source: Absolute path to predicted PDB
        - pred_chain: Chain ID in predicted structure
    """
    # Use both Biopython's sequence alignment and structural alignment
    parser = PDBParser(QUIET=True)
    pp_builder = PPBuilder()  # For extracting sequences

    # Pre-extract CA coordinates and sequences for experimental chains
    # This amortizes the cost across all predicted chains
    exp_struct = parser.get_structure("exp", exp_path)
    exp_chains = {}
    exp_sequences = {}

    for ch in exp_struct.get_chains():
        # Extract CA coordinates for structural alignment
        exp_chains[ch.id] = {res.id[1]: res["CA"].get_coord() for res in ch if res.id[0] == " " and "CA" in res}
        # Extract amino acid sequence for sequence alignment
        seq = "".join(str(pp.get_sequence()) for pp in pp_builder.build_peptides(ch))
        if seq:  # Only store if sequence is not empty
            exp_sequences[ch.id] = seq

    results: Dict[str, Tuple[str, Dict]] = {}
    # Track used experimental chains to avoid duplicate assignments in multimers
    used_exp_chains = set()

    for pred_path in sorted(glob.glob(os.path.join(pred_dir, "*.pdb"))):
        pred_name = os.path.basename(pred_path)
        pred_struct = parser.get_structure(pred_name, pred_path)

        # Iterate over all chains in the predicted structure
        for p_chain in pred_struct.get_chains():
            # Extract CA coordinates for this predicted chain
            p_ca = {res.id[1]: res["CA"].get_coord() for res in p_chain if res.id[0] == " " and "CA" in res}
            
            # Get sequence for this predicted chain
            pred_seq = "".join(str(pp.get_sequence()) for pp in pp_builder.build_peptides(p_chain))
            if not pred_seq:
                continue  # Skip chains with no sequence
                
            best_rmsd = float("inf")
            best_exp_chain = None
            best_score = float("-inf")
            pairs = []

            # First pass: use sequence alignment to find the best matching chain
            for exp_id, exp_seq in exp_sequences.items():
                if exp_id in used_exp_chains:
                    continue  # Skip already assigned chains for better multimer handling
                    
                # Perform sequence alignment
                alignment = pairwise2.align.globalxx(pred_seq, exp_seq, one_alignment_only=True)
                score = alignment[0].score if alignment else 0
                
                # Track best sequence match
                if score > best_score:
                    best_score = score
                    
                    # Calculate RMSD for chains with common residues
                    if exp_id in exp_chains and p_ca:
                        shared = sorted(p_ca.keys() & exp_chains[exp_id].keys())
                        if shared:
                            coords_pred = np.array([p_ca[i] for i in shared])
                            coords_exp = np.array([exp_chains[exp_id][i] for i in shared])
                            rmsd_val = rms.rmsd(coords_pred, coords_exp, superposition=True)
                            
                            pairs.append({
                                "pred_chain": p_chain.id,
                                "exp_chain": exp_id,
                                "rmsd": float(rmsd_val),
                                "seq_score": score
                            })
                            
                            if rmsd_val < best_rmsd:
                                best_rmsd = rmsd_val
                                best_exp_chain = exp_id

            # Second pass: if sequence alignment didn't work, try pure structural alignment
            if best_exp_chain is None:
                for exp_id, exp_ca in exp_chains.items():
                    if exp_id in used_exp_chains:
                        continue
                        
                    shared = sorted(p_ca.keys() & exp_ca.keys())
                    if not shared:
                        continue
                        
                    coords_pred = np.array([p_ca[i] for i in shared])
                    coords_exp = np.array([exp_ca[i] for i in shared])
                    rmsd_val = rms.rmsd(coords_pred, coords_exp, superposition=True)
                    
                    pairs.append({
                        "pred_chain": p_chain.id,
                        "exp_chain": exp_id,
                        "rmsd": float(rmsd_val)
                    })
                    
                    if rmsd_val < best_rmsd:
                        best_rmsd = rmsd_val
                        best_exp_chain = exp_id

            if best_exp_chain is not None:
                # Mark this experimental chain as assigned to prevent duplicate usage
                used_exp_chains.add(best_exp_chain)
                
                key = f"{pred_name}:{p_chain.id}"  # Unique identifier for this chain pair
                results[key] = (
                    best_exp_chain,
                    {
                        "rmsd": best_rmsd,
                        "pairs": pairs,
                        "source": pred_path,      # Absolute path for later use
                        "pred_chain": p_chain.id,
                    },
                )

    return results

# ---------------------------------------------------------------------------
# Strict PDB reorder utility 
# ---------------------------------------------------------------------------

def reorder_residues_in_pdb(pdb_input: str, pdb_output: str) -> None:
    """
    Reorder and renumber atoms in a PDB file for consistency.
    
    Ensures proper ordering of atoms by:
    1. Grouping atoms by chain and residue
    2. Sorting by chain ID, residue number, and insertion code
    3. Renumbering atoms sequentially
    4. Inserting TER records between chains
    5. Adding proper END record
    
    This normalization is important for downstream processing and
    ensures compatibility with various structure analysis tools.
    
    Parameters
    ----------
    pdb_input : str
        Path to input PDB file
    pdb_output : str
        Path to output reordered PDB file
    """
    from collections import defaultdict

    with open(pdb_input, "r") as f:
        lines = f.readlines()

    atom_lines = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    other_lines = [l for l in lines if not l.startswith(("ATOM", "HETATM", "TER"))]

    # Group atoms per (chain, resid, insertion_code, resname)
    residues = defaultdict(list)
    for line in atom_lines:
        chain_id = line[21]
        res_seq = int(line[22:26])
        i_code = line[26] if line[26] != ' ' else ''
        resname = line[17:20]
        key = (chain_id, res_seq, i_code, resname)
        residues[key].append(line)

    # Sort keys by chain, resid, insertion code
    sorted_keys = sorted(residues.keys(), key=lambda x: (x[0], x[1], x[2]))

    new_lines = []
    serial = 1
    prev_chain = None
    for key in sorted_keys:
        chain_id = key[0]
        if prev_chain and chain_id != prev_chain:
            new_lines.append("TER\n")
        for line in sorted(residues[key], key=lambda l: l[12:16]):
            new_line = f"{line[:6]}{serial:5d}{line[11:]}"[:80] + "\n"
            new_lines.append(new_line)
            serial += 1
        prev_chain = chain_id

    new_lines.append("TER\n")
    new_lines += [l for l in other_lines if not l.strip().startswith("END")]
    new_lines.append("END\n")

    with open(pdb_output, "w") as f:
        f.writelines(new_lines)

# ---------------------------------------------------------------------------
# PDB manipulation helpers
# ---------------------------------------------------------------------------

class _ChainSelector(Select):
    """
    BioPython selector to extract a specific chain up to a maximum residue number.
    
    Used with PDBIO to write only the desired chain and residue range to a file.
    This is essential for isolating individual chains in multi-chain structures.
    """

    def __init__(self, wanted: str, maxr: int):
        """
        Initialize the chain selector.
        
        Parameters
        ----------
        wanted : str
            Chain ID to select
        maxr : int
            Maximum residue number (inclusive)
        """
        super().__init__()
        self.wanted = wanted; self.max = maxr

    def accept_residue(self, res):
        """Determine if a residue should be included in the output."""
        return res.get_parent().id == self.wanted and res.id[0] == " " and res.id[1] <= self.max


def _clean_chain(src: str, chain: str, maxr: int, parser: Optional[PDBParser] = None) -> str:
    """
    Extract a specific chain from a PDB file and write to a temporary file.
    
    Isolates a single chain up to the specified maximum residue number,
    which is useful for processing multi-chain structures one chain at a time.
    
    Parameters
    ----------
    src : str
        Path to source PDB file
    chain : str
        Chain ID to extract
    maxr : int
        Maximum residue number to include
    parser : PDBParser, optional
        PDB parser to reuse (default: create new parser)
        
    Returns
    -------
    str
        Path to temporary file containing the extracted chain
    """
    if parser is None:
        parser = PDBParser(QUIET=True)
    
    # Create a persistent temporary file
    tmp = create_persistent_temp_file(suffix=".pdb")
    io = PDBIO()
    io.set_structure(parser.get_structure("s", src))
    io.save(tmp, _ChainSelector(chain, maxr))
    return tmp


def _pad_pdb(src: str) -> str:
    """
    Normalize PDB line format and add proper chain separators.
    
    Ensures all ATOM/HETATM lines are exactly 80 characters and adds TER
    records between chains. Removes ANISOU, TER, and END records for clean
    processing. This normalization is required for consistent parsing.
    
    Parameters
    ----------
    src : str
        Path to source PDB file
        
    Returns
    -------
    str
        Path to temporary file with normalized format
    """
    # Create a persistent temporary file
    tmp = create_persistent_temp_file(suffix=".pdb")
    prev_chain = None
    with open(src) as fin, open(tmp, "w") as fout:
        for ln in fin:
            if ln.startswith(("ATOM  ", "HETATM")):
                chain = ln[21]
                if prev_chain and chain != prev_chain:
                    fout.write("TER\n")
                prev_chain = chain
                fout.write((ln.rstrip("\n") + " " * 80)[:80] + "\n")
            elif ln.startswith(("ANISOU", "TER", "END")):
                continue
            else:
                fout.write(ln)
    return tmp


def _trim_chain_to_universe(src: str, chain_id: str, parser: Optional[PDBParser] = None) -> Tuple[mda.Universe, List[str]]:
    """
    Create an MDAnalysis Universe for a specific chain with cleanup tracking.
    
    Extracts a single chain from a PDB file and creates an MDAnalysis Universe
    object for structural analysis. Tracks temporary files for later cleanup.
    
    Parameters
    ----------
    src : str
        Path to the PDB file
    chain_id : str
        Chain identifier to extract
    parser : PDBParser, optional
        PDB parser to reuse (default: create new parser)
        
    Returns
    -------
    universe : MDAnalysis.Universe
        Universe containing only the specified chain
    temp_files : list
        List of temporary files that should be cleaned up later
    """
    if parser is None:
        parser = PDBParser(QUIET=True)
    
    # Create temporary files that will persist until explicitly cleaned up
    temp_files = []
    
    # Clean the chain
    tmp = _clean_chain(src, chain_id, _get_max_resid(src, chain_id, parser), parser)
    temp_files.append(tmp)
    
    # Pad the PDB
    padded = _pad_pdb(tmp)
    temp_files.append(padded)
    
    # Create the universe
    u = mda.Universe(padded)
    
    return u, temp_files

# ---------------------------------------------------------------------------
# Reliable RMSD calculation (altloc-aware)
# ---------------------------------------------------------------------------

def _safe_rmsd(u1: mda.Universe, u2: mda.Universe, sel: str) -> Tuple[float, int]:
    """
    Calculate RMSD between two structures handling alternate locations.
    
    Computes root-mean-square deviation of CA atoms between two structures,
    properly handling alternate location indicators (altloc) in PDB files.
    Prioritizes atoms with no altloc or altloc 'A'.
    
    Parameters
    ----------
    u1 : MDAnalysis.Universe
        First structure (mobile)
    u2 : MDAnalysis.Universe
        Second structure (reference)
    sel : str
        Atom selection string
        
    Returns
    -------
    rmsd : float
        RMSD value in Angstroms
    n_atoms : int
        Number of common CA atoms used in calculation
        
    Raises
    ------
    RuntimeError
        If no common residues are found
    """
    def _pick(u):
        ag = u.select_atoms(sel + ' and altloc ""')
        if not ag.n_atoms:
            ag = u.select_atoms(sel + ' and altloc "A"')
        return ag if ag.n_atoms else u.select_atoms(sel)

    mob, ref = _pick(u1), _pick(u2)
    mob_map = {a.resid: a.position for a in mob if a.name == "CA"}
    ref_map = {a.resid: a.position for a in ref if a.name == "CA"}
    common = sorted(mob_map.keys() & ref_map.keys())
    if not common:
        raise RuntimeError("No common residues for RMSD")
    mob_pos = np.vstack([mob_map[r] for r in common])
    ref_pos = np.vstack([ref_map[r] for r in common])
    return float(rms.rmsd(mob_pos, ref_pos, superposition=False)), len(common)

# ---------------------------------------------------------------------------
# Helper – ensure patched segment has correct chain IDs
# ---------------------------------------------------------------------------

def _reassign_chain(seg: mda.AtomGroup, chain_id: str) -> None:
    """
    Assign chain ID to all atoms in a segment for proper structure assembly.
    
    Sets both the chainID (atom-level) and segid (segment-level) attributes
    to ensure the chain assignment is recognized across different MDAnalysis
    versions and PDB writing operations.
    
    Parameters
    ----------
    seg : MDAnalysis.AtomGroup
        Atoms to reassign
    chain_id : str
        Target chain identifier
    """
    # atoms: chainID attribute exists in topology
    for atom in seg.atoms:
        atom.chainID = chain_id
    # segments: proper place for segid
    for segment in seg.segments:
        segment.segid = chain_id


# ---------------------------------------------------------------------------
# Enhanced patching helpers for fusion proteins
# ---------------------------------------------------------------------------

def _is_fusion_protein(pred_chain, exp_chain) -> bool:
    """
    Detect if the experimental chain is a fusion protein (longer than predicted).
    
    Fusion proteins contain additional sequences beyond the target protein,
    often used for purification or stabilization. This function identifies
    such cases by comparing maximum residue numbers.
    
    Parameters
    ----------
    pred_chain : Bio.PDB.Chain
        Predicted chain from AlphaFold
    exp_chain : Bio.PDB.Chain
        Experimental chain from crystal structure
        
    Returns
    -------
    bool
        True if experimental chain extends beyond the predicted chain
    """
    if not pred_chain or not exp_chain:
        return False
        
    # Get maximum residue ID for each chain
    try:
        pred_max = max(res.id[1] for res in pred_chain.get_residues() if res.has_id('CA'))
        exp_max = max(res.id[1] for res in exp_chain.get_residues() if res.has_id('CA'))
        return exp_max > pred_max
    except ValueError:  # Empty chain
        return False


def _truncate_fusion_protein(exp_chain, max_resid) -> None:
    """
    Truncate fusion protein experimental chain to match predicted length.
    
    Removes extra residues beyond the target protein sequence to ensure
    proper alignment with the predicted structure.
    
    Parameters
    ----------
    exp_chain : Bio.PDB.Chain
        Experimental chain to truncate
    max_resid : int
        Maximum residue number to keep (from predicted structure)
    """
    # Remove residues beyond the maximum residue ID of the predicted chain
    for res in list(exp_chain.get_residues()):
        if res.id[1] > max_resid:
            exp_chain.detach_child(res.id)

# ---------------------------------------------------------------------------
# Alignment with reference structure for orientation
# ---------------------------------------------------------------------------

def align_with_reference(structure_path: str, reference_path: str, output_path: str) -> Tuple[float, float]:
    """
    Align a structure to a reference structure using CA atom superposition.
    
    Performs rigid-body superposition of two protein structures based on their
    CA (alpha-carbon) atoms. This is used for orienting structures to a common
    reference frame, typically using a structure containing DUM orientation markers.
    
    Parameters
    ----------
    structure_path : str
        Path to the structure to be aligned
    reference_path : str
        Path to the reference structure
    output_path : str
        Path where the aligned structure will be saved
        
    Returns
    -------
    tuple
        (RMSD before alignment, RMSD after alignment) both in Angstroms
        
    Raises
    ------
    ValueError
        If no CA atoms are found in either structure
    """
    parser = PDBParser(QUIET=True)
    
    # Load structures
    ref_struct = parser.get_structure('ref', reference_path)
    tgt_struct = parser.get_structure('target', structure_path)
    
    # Get CA atoms
    ref_atoms = [a for a in ref_struct[0].get_atoms() if a.get_id() == 'CA']
    tgt_atoms = [a for a in tgt_struct[0].get_atoms() if a.get_id() == 'CA']
    
    # Use the minimum number of atoms
    n = min(len(ref_atoms), len(tgt_atoms))
    
    if n == 0:
        raise ValueError("No CA atoms found in one or both structures")
    
    # Calculate RMSD before alignment
    ref_coords = np.array([a.get_coord() for a in ref_atoms[:n]])
    tgt_coords = np.array([a.get_coord() for a in tgt_atoms[:n]])
    rms_before = np.sqrt(((ref_coords - tgt_coords)**2).sum(axis=1).mean())
    
    # Superimpose
    sup = Superimposer()
    sup.set_atoms(ref_atoms[:n], tgt_atoms[:n])
    rms_after = sup.rms
    
    # Apply transformation
    sup.apply(tgt_struct.get_atoms())
    
    # Save aligned structure
    io = PDBIO()
    io.set_structure(tgt_struct)
    io.save(output_path)
    
    log.info(f"Aligned structure: RMSD before={rms_before:.3f}Å, after={rms_after:.3f}Å -> {output_path}")
    
    return rms_before, rms_after

# ---------------------------------------------------------------------------
# Main alignment and patching function
# ---------------------------------------------------------------------------

def align_and_patch(
    af_path: str,
    exp_path: str,
    gene: str,
    prot_len: int,  # Kept for API compatibility
    window: int,    # Kept for API compatibility
    outdir: str,
    partner_ids: Optional[List[str]] = None,
    uniprot_id: Optional[str] = None,
    reference_paths: Optional[Dict[str, str]] = None
) -> dict:
    """
    Align predicted structures onto experimental structure and patch missing regions.
    
    This is the main function that orchestrates the complete alignment and patching
    workflow. It:
    
    1. Finds optimal chain correspondence between predicted and experimental structures
    2. Performs structural superposition using MDAnalysis
    3. Detects gaps (missing regions) in the experimental structure
    4. Patches gaps using coordinates from the predicted structure
    5. Handles both single-chain and multi-chain complexes (homo/hetero-oligomers)
    6. Detects and handles fusion proteins
    7. Calculates quality metrics (RMSD, gap statistics)
    8. Manages temporary files with automatic cleanup
    
    The function uses a sophisticated chain-pairing algorithm that combines
    sequence similarity and structural RMSD to find the best correspondence,
    which is essential for multi-chain complexes.
    
    Parameters
    ----------
    af_path : str
        Directory containing AlphaFold predicted PDB files
    exp_path : str
        Path to experimental PDB file
    gene : str
        Gene name for output file naming
    prot_len : int
        Protein length (retained for API compatibility, not actively used)
    window : int
        Window size (retained for API compatibility, not actively used)
    outdir : str
        Output directory for aligned and patched structures
    partner_ids : Optional[List[str]], optional
        List of UniProt IDs for interaction partners (for hetero-oligomers)
    uniprot_id : Optional[str], optional
        UniProt ID of the main protein
    reference_paths : Optional[Dict[str, str]], optional
        Dictionary mapping UniProt IDs to reference structure paths
        for subsequent orientation
        
    Returns
    -------
    dict
        Comprehensive results dictionary containing:
        - gene: Gene name
        - pdb_id: PDB ID with hyperlink
        - aligned_pdb: Path to aligned structure
        - patched_pdb: Path to patched structure
        - patched_with_hetatm_pdb: Patched structure with ligands
        - oriented_pdb: Path to oriented structure (if reference provided)
        - fusion_protein: Boolean indicating fusion protein detection
        - is_heteromer: Boolean for hetero-oligomer
        - interaction_partners: List of partner IDs
        - n_CA_common_before_patch: Number of aligned CA atoms before patching
        - n_CA_common_after_patch: Number of aligned CA atoms after patching
        - gaps_found: String listing all detected gaps
        - gaps_used: String listing gaps that were patched
        - gap_count: Total number of gaps
        - nb_atoms_to_patch: Total atoms in gap regions
        - patched_atoms: Number of atoms successfully patched
        - %_gaps_patched: Percentage of gaps patched
        - rmsd_exp: RMSD between predicted and experimental (Å)
        - rmsd_patched: RMSD between predicted and patched (Å)
        - rmsd_oriented: RMSD after orientation (Å)
        - rmsd_pairs: Per-chain RMSD values
        - status: 'ok' or 'error'
        
    Raises
    ------
    RuntimeError
        If no PDB files found in af_path or if no chains can be paired
    """
    # List to track all temporary files for cleanup at the end
    all_temp_files = []

    try:
        # ---------------------------------------------------------------------
        # Discover and validate prediction files
        # ---------------------------------------------------------------------
        if not glob.glob(os.path.join(af_path, "*.pdb")):
            raise RuntimeError(f"No PDB files found in {af_path!r}")

        mapping = _find_best_chain_pair(af_path, exp_path)
        valid = [(k, v) for k, v in mapping.items() if v[0] is not None]
        if not valid:
            raise RuntimeError("No chain could be paired with experimental model")

        # Extract PDB ID and create hyperlink for Excel output
        pdbid_raw = os.path.splitext(os.path.basename(exp_path))[0].split("_")[-1].upper()
        pdbid = f'=HYPERLINK("https://www.rcsb.org/structure/{pdbid_raw}", "{pdbid_raw}")'
        
        # Initialize results dictionary with all expected fields
        res: Dict = {
            "gene": gene,
            "pdb_id": pdbid,
            "aligned_pdb": "",
            "patched_pdb": "",
            "patched_with_hetatm_pdb": "",
            "oriented_pdb": "",  # New field for the oriented structure
            "fusion_protein": False,
            "is_heteromer": False,
            "interaction_partners": "",
            "n_CA_common_before_patch": 0,
            "n_CA_common_after_patch": 0,
            "gaps_found": "",
            "gaps_used": "",
            "gap_count": 0,
            "nb_atoms_to_patch": 0,
            "patched_atoms": 0,
            "%_gaps_patched": 0.0,
            "rmsd_exp": None,
            "rmsd_patched": None,
            "rmsd_oriented": None,  # New field for orientation RMSD
            "rmsd_pairs": {},
            "status": "error"
        }

        # ---------------------------------------------------------------------
        # Analyze experimental structure and determine oligomeric state
        # ---------------------------------------------------------------------
        # Parse experimental structure to identify chains
        exp_parser = PDBParser(QUIET=True)
        exp_structure = exp_parser.get_structure("exp", exp_path)
        exp_chains = list(exp_structure[0].get_chains())
        exp_chain_ids = [ch.id for ch in exp_chains]

        log.info(f"Experimental structure: {pdbid} with chains {exp_chain_ids}")
        log.info(f"Partner IDs: {partner_ids} (type: {type(partner_ids).__name__})")
        if uniprot_id:
            log.info(f"UniProt ID: {uniprot_id}")

        # Normalize partner_ids to list format regardless of chain count
        if partner_ids is None:
            partner_ids = []
        elif isinstance(partner_ids, str):
            partner_ids = [partner_ids]
        
        # Determine oligomeric state and log chain assignments
        if len(exp_chains) > 1:
            # Multi-chain structure: log chain correspondences for transparency
            log.info("Chain-protein correspondences:")
            for key, (exp_chain, details) in mapping.items():
                pred_file = os.path.basename(details["source"])
                pred_chain = details["pred_chain"]
                log.info(f"• {pred_file}:{pred_chain} → {pdbid}:{exp_chain} (RMSD: {details.get('rmsd', 'N/A')})")
            
            # Map partners to chains if information is available
            if partner_ids and len(partner_ids) <= len(exp_chain_ids):
                for i, (partner, chain) in enumerate(zip(partner_ids, exp_chain_ids)):
                    log.info(f"• UniProt {partner} → Chain {chain}")
        
        # Set oligomeric state metadata
        res["interaction_partners"] = partner_ids
        res["is_heteromer"] = bool(partner_ids)
        log.info(f"Structure defined as: {'HETEROMER' if res['is_heteromer'] else 'HOMOMER'}")
        
        # ---------------------------------------------------------------------
        # Load experimental structure and prepare for gap patching
        # ---------------------------------------------------------------------
        # Create normalized experimental structure (protein only)
        padded_exp = _pad_pdb(exp_path)
        all_temp_files.append(padded_exp)
        
        hybrid = mda.Universe(padded_exp)

        # Collector for per-chain universes and RMSDs
        predicted: List[Tuple[str, str, mda.Universe, List[str]]] = []

        # ---------------------------------------------------------------------
        # Process each chain pair: alignment and fusion detection
        # ---------------------------------------------------------------------
        for key, (exp_chain, details) in valid:
            src_pred = details["source"]
            pred_chain = details["pred_chain"]

            # Check for fusion protein and truncate if necessary
            pred_struct = exp_parser.get_structure("pred", src_pred)
            for p_chain in pred_struct.get_chains():
                if p_chain.id == pred_chain:
                    for e_chain in exp_chains:
                        if e_chain.id == exp_chain:
                            if _is_fusion_protein(p_chain, e_chain):
                                res["fusion_protein"] = True
                                pred_max = max(res.id[1] for res in p_chain.get_residues() if res.has_id('CA'))
                                _truncate_fusion_protein(e_chain, pred_max)
                            break
                    break

            # Use the corrected _trim_chain_to_universe function that returns temp files
            uA, temp_files_A = _trim_chain_to_universe(src_pred, pred_chain, exp_parser)
            uE, temp_files_E = _trim_chain_to_universe(exp_path, exp_chain, exp_parser)
            
            # Add temporary files to the tracking list
            all_temp_files.extend(temp_files_A)
            all_temp_files.extend(temp_files_E)

            # Calculate common residues once and reuse
            ca_resids_A = set(uA.select_atoms("name CA").resids)
            ca_resids_E = set(uE.select_atoms("name CA").resids)
            common = sorted(ca_resids_A & ca_resids_E)
            
            if not common:
                continue  # nothing to align – rare but safe
                
            sel = "protein and name CA and resid " + " ".join(map(str, common))
            
            try:
                align.AlignTraj(uA, uE, select=sel, in_memory=True, weights="mass").run()
                rms_val, n_common = _safe_rmsd(uA, uE, sel)
                rms_val = round(rms_val, 3)     
                res["rmsd_pairs"][exp_chain] = rms_val
                predicted.append((pred_chain, exp_chain, uA, temp_files_A))
                res["n_CA_common_before_patch"] += len(common)
            except (ValueError, RuntimeError) as e:
                log.warning(f"Failed alignment for chain {pred_chain}-{exp_chain}: {e}")
                continue
        
        if predicted:
            gene_dir = os.path.join(outdir, gene)
            os.makedirs(gene_dir, exist_ok=True)
            aligned_path = os.path.join(
                gene_dir, f"aligned_on_{os.path.basename(exp_path)}"
            )

            ag_pred = [uA.select_atoms("protein") for pred_chain, exp_chain, uA, _ in predicted]

            merged = mda.Merge(*ag_pred)      
            merged.atoms.write(aligned_path) 

            res["aligned_pdb"] = aligned_path


        # ---------------------------------------------------------------------
        # Detect and patch gaps per chain
        # ---------------------------------------------------------------------
        for pred_chain, exp_chain, uA, _ in predicted:
            log.info(f"Checking gaps for predicted chain {pred_chain} against experimental chain {exp_chain}")

            ids_A = set(uA.select_atoms("name CA").resids)
            ids_H = set(hybrid.select_atoms(f"chainID {exp_chain} and name CA").resids)

            extras = sorted(ids_A - ids_H)
            gaps = group_segments(extras) if extras else []

            res["gap_count"] += len(gaps)
            if gaps:
                res["gaps_found"] += ";".join(f"{exp_chain}:{s}-{e}" for s, e in gaps) + ";"

            patched_ranges = []  # For building the dynamic list per pair
            
            # Pre-calculate total atoms to patch
            total_atoms_to_patch = 0
            for start, end in gaps:
                seg = uA.select_atoms(f"protein and resid {start}:{end}")
                if seg.n_atoms:
                    total_atoms_to_patch += seg.n_atoms
            
            res["nb_atoms_to_patch"] += total_atoms_to_patch

            for start, end in gaps:
                seg = uA.select_atoms(f"protein and resid {start}:{end}")
                if not seg.n_atoms:
                    continue
                _reassign_chain(seg, exp_chain)
                hybrid = mda.Merge(hybrid.select_atoms("protein"), seg.atoms)
                res["patched_atoms"] += seg.n_atoms
                patched_ranges.append(f"{exp_chain}:{start}-{end}")

            if patched_ranges:
                res["gaps_used"] += ";".join(patched_ranges) + ";"

        # ---------------------------------------------------------------------
        # Write patched structure
        # ---------------------------------------------------------------------
        patched_path = os.path.join(gene_dir, f"patched_{os.path.basename(exp_path)}")
        
        # CORRECTION: For homomers, keep only one chain
        if not res["is_heteromer"] and predicted:
            # For homomers, keep only the first experimental chain
            first_exp_chain = predicted[0][1]  # Get the first experimental chain
            log.info(f"Homomer detected: keeping only chain {first_exp_chain}")
            # Select only atoms from the first chain
            hybrid.select_atoms(f"protein and chainID {first_exp_chain}").write(patched_path)
        else:
            # Original behavior for heteromers
            hybrid.select_atoms("protein").write(patched_path)

        # Reorder and rewrite the same file directly
        reorder_residues_in_pdb(patched_path, patched_path)
        res["patched_pdb"] = patched_path

        # ---------------------------------------------------------------------
        # Unified RMSD calculation (experimental and patched)
        # ---------------------------------------------------------------------
        rmsd_exp_vals = []
        rmsd_patch_vals = []

        for pred_chain, exp_chain, uA, _ in predicted:
            try:
                # Experimental comparison
                uE, temp_files_E = _trim_chain_to_universe(exp_path, exp_chain, exp_parser)
                all_temp_files.extend(temp_files_E)
                
                sel_common = "name CA"
                mob = uA.select_atoms(sel_common)
                ref = uE.select_atoms(sel_common)
                mob_map = {a.resid: a.position for a in mob}
                ref_map = {a.resid: a.position for a in ref}
                common = sorted(mob_map.keys() & ref_map.keys())
                if common:
                    rmsd = float(rms.rmsd(
                        np.vstack([mob_map[r] for r in common]),
                        np.vstack([ref_map[r] for r in common]),
                        superposition=False
                    ))
                    rmsd_exp_vals.append((rmsd, len(common)))
                    log.info(f"→ RMSD EXP: pred={pred_chain}, exp={exp_chain}, n_CA={len(common)}, rmsd={rmsd:.3f}")

                # Patched comparison
                ref = hybrid.select_atoms(f"chainID {exp_chain} and name CA")
                ref_map = {a.resid: a.position for a in ref}
                common = sorted(mob_map.keys() & ref_map.keys())
                if common:
                    rmsd = float(rms.rmsd(
                        np.vstack([mob_map[r] for r in common]),
                        np.vstack([ref_map[r] for r in common]),
                        superposition=False
                    ))
                    rmsd_patch_vals.append((rmsd, len(common)))
                    res["n_CA_common_after_patch"] += len(common)
                    log.info(f"→ RMSD PATCHED: pred={pred_chain}, patched, n_CA={len(common)}, rmsd={rmsd:.3f}")
            except (ValueError, RuntimeError) as e:
                log.warning(f"Failed RMSD computation for chain {pred_chain}-{exp_chain}: {e}")
            except Exception as e:
                log.error(f"Unexpected error during RMSD computation: {e}")
                # Don't re-raise to maintain compatibility

        # Final RMSD aggregation
        if rmsd_exp_vals:
            total = sum(n for _, n in rmsd_exp_vals)
            res["rmsd_exp"] = round(sum(r * n for r, n in rmsd_exp_vals) / total, 3)
            log.info(f"✓ Global experimental RMSD: {res['rmsd_exp']:.3f}")

        if rmsd_patch_vals:
            total = sum(n for _, n in rmsd_patch_vals)
            res["rmsd_patched"] = round(sum(r * n for r, n in rmsd_patch_vals) / total, 3)
            log.info(f"✓ Global patched RMSD: {res['rmsd_patched']:.3f}")

        # ---------------------------------------------------------------------
        # HETATM reintegration (optional)
        # ---------------------------------------------------------------------
        with open(exp_path) as f:
            het_atoms = [ln for ln in f if ln.startswith("HETATM")]
            
        if het_atoms:
            hetout = os.path.join(gene_dir, f"patched_with_hetatm_{os.path.basename(exp_path)}")
            with open(patched_path) as p_in, open(hetout, "w") as p_out:
                p_out.writelines([ln for ln in p_in if not ln.startswith("END")] + het_atoms + ["END\n"])
            res["patched_with_hetatm_pdb"] = hetout
        # ---------------------------------------------------------------------
        # Final metrics
        # ---------------------------------------------------------------------
        if res["gap_count"]:
            done = len(res["gaps_used"].rstrip(";").split(";")) if res["gaps_used"] else 0
            res["%_gaps_patched"] = round(100 * done / res["gap_count"], 1)
        res["status"] = "ok"

        return res
    
    finally:
        # Clean up all temporary files
        cleanup_temp_files(all_temp_files)
