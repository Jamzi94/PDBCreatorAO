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
    Internal: fetch the list of AlphaFold predictions for a UniProt ID.
    Returns a list (possibly empty) of prediction dicts.
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
    Return the AlphaFold 'entryId' (e.g. 'AF-P14859-F1') for the given UniProt ID,
    or an empty string if none is available.
    """
    preds = _fetch_predictions(uniprot_id)
    if not preds:
        log.info(f"No AlphaFold prediction available for {uniprot_id}")
        return ""
    entry_id = preds[0].get("entryId", "")
    return entry_id

def download_alphafold_pdb(uniprot_id: str, output_path: str) -> bool:
    """
    Download the AlphaFold PDB file for uniprot_id to the given path.
    Returns True on success, False otherwise.
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
    Retrieve summary details of the AlphaFold model for the given UniProt ID.
    Returns a dict with keys:
      method, resolution, methodology, conformation,
      count, descs, chains, orgs
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

    # You could decide to infer conformation from templateModelQuality or other fields
    # but by default AlphaFold predictions have no explicit 'conformation'.
    entry_id = p.get("entryId", "")
    count = 1 if entry_id else 0

    # Description: use the protein name passed in
    descs = prot_name

    # Chains: AlphaFold models are single‐chain, chain ID 'A'
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
