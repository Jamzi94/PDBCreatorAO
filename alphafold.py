"""
alphafold.py

Interact with the AlphaFold database API to retrieve predicted protein structures.

This module provides functions to:
 - Query the AlphaFold database for available predictions
 - Download AlphaFold predicted PDB structures
 - Extract prediction metadata and quality metrics

The AlphaFold database provides AI-predicted 3D structures for millions of proteins,
which serve as complementary data to experimental structures.
"""

import os
import requests
import logging

log = logging.getLogger(__name__)

# If you have an API key, set it as environment variable ALPHAFOLD_API_KEY;
# otherwise you can omit the key parameter.
API_KEY = os.getenv("AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94", "")
API_PRED_TEMPLATE = (
    "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    + ("?key={api_key}" if API_KEY else "")
)

def _fetch_predictions(uniprot_id: str) -> list:
    """
    Fetch the list of AlphaFold predictions for a UniProt accession.
    
    Internal helper function that queries the AlphaFold API to retrieve
    all available predictions for the specified protein.
    
    Parameters
    ----------
    uniprot_id : str
        UniProt accession ID
        
    Returns
    -------
    list
        List of prediction dictionaries (empty list if none found or error occurs)
    """
    url = API_PRED_TEMPLATE.format(uniprot_id=uniprot_id, api_key=API_KEY)
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        if not isinstance(data, list):
            log.error(f"Unexpected AlphaFold response format for {uniprot_id}")
            return []
        return data
    except Exception as e:
        log.error(f"Error fetching AlphaFold prediction for {uniprot_id}: {e}")
        return []

def get_alphafold_pdb_info(uniprot_id: str) -> str:
    """
    Get the AlphaFold entry identifier for a protein.
    
    Retrieves the AlphaFold database entry ID (e.g., 'AF-P14859-F1') which
    uniquely identifies the predicted structure.
    
    Parameters
    ----------
    uniprot_id : str
        UniProt accession ID
        
    Returns
    -------
    str
        AlphaFold entry identifier, or empty string if no prediction exists
    """
    preds = _fetch_predictions(uniprot_id)
    if not preds:
        log.info(f"No AlphaFold prediction available for {uniprot_id}")
        return ""
    entry_id = preds[0].get("entryId", "")
    return entry_id

def download_alphafold_pdb(uniprot_id: str, output_path: str) -> bool:
    """
    Download the AlphaFold predicted structure file in PDB format.
    
    Retrieves the PDB file for the AlphaFold prediction and saves it to
    the specified path. The file contains atomic coordinates and confidence
    scores (pLDDT) for the predicted structure.
    
    Parameters
    ----------
    uniprot_id : str
        UniProt accession ID
    output_path : str
        Path where the PDB file should be saved
        
    Returns
    -------
    bool
        True if download succeeded, False otherwise
    """
    preds = _fetch_predictions(uniprot_id)
    if not preds or not preds[0].get("pdbUrl"):
        log.warning(f"No PDB URL for AlphaFold model of {uniprot_id}")
        return False
    pdb_url = preds[0]["pdbUrl"]
    try:
        resp = requests.get(pdb_url, timeout=10)
        resp.raise_for_status()
        with open(output_path, "wb") as f:
            f.write(resp.content)
        log.info(f"Downloaded AlphaFold PDB for {uniprot_id} to {output_path}")
        return True
    except Exception as e:
        log.error(f"Error downloading AlphaFold PDB for {uniprot_id}: {e}")
        return False

def get_alphafold_details(uniprot_id: str, prot_name: str) -> dict:
    """
    Retrieve summary metadata for an AlphaFold predicted structure.
    
    Gathers basic information about the AlphaFold prediction to provide
    consistent metadata fields comparable to experimental structures.
    Since AlphaFold models are computational predictions, some fields
    (resolution, conformation, environment) are marked as not applicable.
    
    Parameters
    ----------
    uniprot_id : str
        UniProt accession ID
    prot_name : str
        Protein name to use in the description
        
    Returns
    -------
    dict
        Dictionary containing:
        - method: 'Prediction'
        - resolution: 'N/A' (not applicable for predictions)
        - methodology: 'AlphaFold'
        - conformation: 'Not specified'
        - environment: 'NA'
        - count: 1 if prediction exists, 0 otherwise
        - descs: Protein name/description
        - chains: 'A' (AlphaFold models are single-chain)
        - orgs: Source organism scientific name
    """
    preds = _fetch_predictions(uniprot_id)
    if not preds:
        return {
            "method": "Prediction",
            "resolution": "N/A",
            "methodology": "AlphaFold",
            "conformation": "Not specified",
            "environment": "NA",
            "count": 0,
            "descs": prot_name,
            "chains": "",
            "orgs": ""
        }

    p = preds[0]
    method = "Prediction"
    resolution = "N/A"
    methodology = "AlphaFold"
    conformation = "Not specified"
    environment = "NA"

    # AlphaFold predictions don't have explicit conformational states
    # Conformation could potentially be inferred from structural analysis
    entry_id = p.get("entryId", "")
    count = 1 if entry_id else 0

    # Use the provided protein name as the description
    descs = prot_name

    # AlphaFold models are single-chain structures with chain ID 'A'
    chains = "A"

    # Organism name from the returned JSON
    orgs = p.get("organismScientificName", "")

    return {
        "method": method,
        "resolution": resolution,
        "methodology": methodology,
        "conformation": conformation,
        "environment": environment,
        "count": count,
        "descs": descs,
        "chains": chains,
        "orgs": orgs
    }
