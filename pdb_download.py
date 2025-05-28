import os
import requests
import logging

log = logging.getLogger(__name__)

def download_pdb_structure(pdb_id: str, output_dir: str, filename: str = None) -> str:
    """
    Download the PDB file for pdb_id into output_dir, first trying OPM,
    then falling back to the RCSB files server. If filename is given,
    use it; otherwise use '{pdb_id}.pdb'. Returns the path to the
    downloaded file, or empty string on failure.
    """
    pdb_id = pdb_id.lower()
    os.makedirs(output_dir, exist_ok=True)
    if not filename:
        filename = f"{pdb_id}.pdb"
    output_path = os.path.join(output_dir, filename)

    # 1) Try OPM mirror
    opm_url = f"https://biomembhub.org/shared/opm-assets/pdb/{pdb_id}.pdb"
    try:
        resp = requests.get(opm_url, timeout=10)
        if resp.ok:
            with open(output_path, "wb") as f:
                f.write(resp.content)
            log.info(f"Downloaded PDB {pdb_id} from OPM to {output_path}")
            return output_path
        else:
            log.debug(f"OPM returned status {resp.status_code} for {pdb_id}")
    except Exception as e:
        log.warning(f"OPM download failed for {pdb_id}: {e}")

    # 2) Fallback to RCSB files server
    rcsb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        resp = requests.get(rcsb_url, timeout=10)
        resp.raise_for_status()
        with open(output_path, "wb") as f:
            f.write(resp.content)
        log.info(f"Downloaded PDB {pdb_id} from RCSB to {output_path}")
        return output_path
    except Exception as e:
        log.error(f"Failed to download PDB {pdb_id} from RCSB: {e}")
        return ""
