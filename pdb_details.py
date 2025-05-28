import logging
import requests
import re

# Conformational state keywords
CONF_TERMS = [
    "inward-facing","inward facing", "inward-open",  "inward open",
    "outward-facing","outward facing", "outward-open", "outward open",
    "apo", "occluded", "closed", "open" 
]
# Experimental environment keywords
ENV_TERMS = [
    "nanodisc", "lipid", "micelle", "detergent", "bicelle",
    "liposome", "nanodiscs"
]

# Pre-compile regex patterns (word boundaries + case-insensitive)
CONF_PATTERN = re.compile(
    r'\b(' + '|'.join(map(re.escape, CONF_TERMS)) + r')\b',
    flags=re.IGNORECASE
)
ENV_PATTERN = re.compile(
    r'\b(' + '|'.join(map(re.escape, ENV_TERMS)) + r')\b',
    flags=re.IGNORECASE
)

log = logging.getLogger(__name__)

RCSB_CORE = "https://data.rcsb.org/rest/v1/core"

def get_pdb_details(pdb_id: str) -> dict:
    """
    Retrieve experimental details for pdb_id from RCSB. Returns a dict:
    method, resolution, methodology, conformation, environment,
    count, descs, chains, orgs, pubmed_id.
    """
    pdb_id = pdb_id.lower()
    try:
        # Fetch entry JSON (for PubMed ID and free-text keywords)
        url_entry = f"{RCSB_CORE}/entry/{pdb_id}"
        resp_e = requests.get(url_entry, timeout=10)
        resp_e.raise_for_status()
        entry = resp_e.json()
        container = entry.get("rcsb_entry_container_identifiers", {})
        pubmed_id = container.get("pubmed_id", "")
        # Free-text keywords for environment detection
        entry_keywords = entry.get("struct_keywords", {}).get("text", "")
    except Exception as e:
        log.error(f"Error retrieving RCSB entry for {pdb_id}: {e}")
        return {
            "method": "N/A",
            "resolution": "N/A",
            "methodology": "N/A",
            "conformation": "Not specified",
            "environment": "Not specified",
            "count": 0,
            "descs": "",
            "chains": "",
            "orgs": "",
            "pubmed_id": ""
        }

    try:
        # Fetch struct JSON (for title and controlled keywords)
        url_struct = f"{RCSB_CORE}/struct/{pdb_id}"
        resp_s = requests.get(url_struct, timeout=10)
        resp_s.raise_for_status()
        struct = resp_s.json()
    except Exception:
        struct = {}

    info = entry.get("rcsb_entry_info", {})
    method = info.get("experimental_method", "N/A")
    resolution = info.get("resolution_combined", ["N/A"])[0]
    methodology = info.get("structure_determination_methodology", "N/A")

    # Combine title, controlled keywords, free-text keywords, and citation for scanning
    combined_text = " ".join([
        struct.get("struct", {}).get("title", ""),
#       struct.get("struct_keywords", {}).get("pdbx_keywords", ""), <- pas le end point
        entry_keywords,
        entry.get("struct", {}).get("title", "")
    ]).lower()

    # Extract conformational states via regex
    raw_confs = CONF_PATTERN.findall(combined_text)
    conformation = ", ".join(dict.fromkeys(raw_confs)) if raw_confs else "Not specified"

    # Extract experimental environment terms via regex
    raw_envs = ENV_PATTERN.findall(combined_text)
    environment = ", ".join(dict.fromkeys(raw_envs)) if raw_envs else "Not specified"

    # Retrieve polymer entity descriptions
    count = info.get("polymer_entity_count", 0)
    descs = []
    chains = []
    orgs = []
    for i in range(1, count + 1):
        try:
            url_poly = f"{RCSB_CORE}/polymer_entity/{pdb_id}/{i}"
            resp_p = requests.get(url_poly, timeout=10)
            resp_p.raise_for_status()
            poly = resp_p.json()
            d = (
                poly.get("entity", {}).get("pdbx_description")
                or poly.get("rcsb_polymer_entity", {}).get("pdbx_description")
                or "Inconnue"
            )
            ch_list = poly.get("rcsb_polymer_entity_container_identifiers", {}) \
                          .get("auth_asym_ids", [])
            o = poly.get("rcsb_entity_source_organism", [{}])[0] \
                   .get("scientific_name", "N/A")
            descs.append(d)
            chains.append(",".join(ch_list))
            orgs.append(o)
        except Exception as e:
            log.warning(f"Failed to retrieve polymer_entity {i} for {pdb_id}: {e}")
            continue

    return {
        "method": method,
        "resolution": resolution,
        "methodology": methodology,
        "conformation": conformation,
        "environment": environment,
        "count": len(descs),
        "descs": "; ".join(descs),
        "chains": "; ".join(chains),
        "orgs": "; ".join(orgs),
        "pubmed_id": pubmed_id or ""
    }
