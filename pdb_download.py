import os
import requests
import logging

log = logging.getLogger(__name__)

def download_pdb_structure(pdb_id: str, output_dir: str, filename: str = None) -> str:
    """Downloads a PDB structure from OPM or RCSB.

    This function attempts to download a PDB file from two sources. It first
    tries the OPM (Orientations of Proteins in Membranes) database, which is
    specialized for membrane proteins. If the OPM download fails (e.g., the
    protein is not in OPM or the server is down), it falls back to the main
    RCSB PDB data server.

    The downloaded file is saved in the specified output directory.

    Args:
        pdb_id: The 4-character PDB accession ID (e.g., "1A2B").
        output_dir: The directory where the PDB file will be saved.
        filename: An optional name for the output file. If not provided,
            it defaults to "{pdb_id}.pdb".

    Returns:
        The full path to the downloaded PDB file if successful, otherwise an
        empty string.
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
