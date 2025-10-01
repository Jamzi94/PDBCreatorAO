import os
import requests
import json
import logging
from typing import List, Dict

BASE_URL = "https://rest.uniprot.org/uniprotkb"
FIELDS = ",".join([
    "gene_primary", "protein_name", "xref_pdb", "ft_variant",
    "cc_pharmaceutical", "xref_tcdb", "xref_drugcentral",
    "cc_tissue_specificity", "annotation_score", "cc_cofactor",
    "cc_polymorphism", "sequence", "cc_interaction"
])

log = logging.getLogger(__name__)

# Directory for caching UniProt JSON results
CACHE_DIR = os.path.join(os.getcwd(), "cache", "uniprot")
os.makedirs(CACHE_DIR, exist_ok=True)


def fetch_uniprot_json(uniprot_id: str) -> dict:
    """
    Fetch UniProt entry JSON for the given UniProt ID, with caching.
    """
    uniprot_id = uniprot_id.strip()
    cache_file = os.path.join(CACHE_DIR, f"{uniprot_id}.json")
    if os.path.isfile(cache_file):
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                data = json.load(f)
                return data
        except Exception as e:
            log.warning(f"Failed to load cache for {uniprot_id}: {e}. Fetching fresh.")
    url = f"{BASE_URL}/{uniprot_id}"
    resp = requests.get(url, params={"fields": FIELDS}, headers={"Accept": "application/json"})
    resp.raise_for_status()
    data = resp.json()
    try:
        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(data, f)
    except Exception as e:
        log.warning(f"Failed to write cache for {uniprot_id}: {e}")
    return data


def get_first_comment(comments: List[dict], comment_type: str) -> str:
    """
    Extract the first comment of a given type from UniProt JSON comments.
    """
    for c in comments:
        if c.get("commentType") == comment_type:
            texts = c.get("texts", [])
            if texts:
                val = texts[0]
                if isinstance(val, dict) and "value" in val:
                    return val["value"]
                elif isinstance(val, str):
                    return val
    return ""


def get_interactions(uniprot_id: str) -> dict:
    """
    Retrieve only the cc_interaction block.
    """
    url = f"{BASE_URL}/{uniprot_id}"
    resp = requests.get(url,
                        params={"fields": "cc_interaction"},
                        headers={"accept": "application/json"})
    resp.raise_for_status()
    return resp.json()

from typing import List, Dict
import logging

log = logging.getLogger(__name__)

def find_related_uniprots(entry_uniprot: str, data: dict) -> List[Dict[str, str]]:
    """
    Extract UniProt partners whose ID starts with the same prefix as entry_uniprot.
    """
    related = []
    prefix = entry_uniprot[:5]
    for c in data.get("comments", []):
        if c.get("commentType") != "INTERACTION":
            continue
        for inter in c.get("interactions", []):
            p = inter.get("interactantTwo", {})
            pid = p.get("uniProtKBAccession", "")
            pg  = p.get("geneName", "")
            # Only add entries with both UniProt ID and gene name
            if not pid or "-" in pid or not pg:
                continue
            if pid.startswith(prefix):
                related.append({"uniprot_id": pid, "gene_name": pg})
    return related


def parse_uniprot_data(data: dict) -> dict:
    """
    Parse UniProt JSON data and extract standard fields plus interactions.
    """
    # Gene name
    gene = ""
    try:
        gene = data.get("genes", [{}])[0].get("geneName", {}).get("value", "")
    except Exception:
        log.debug("No gene name found in UniProt data")
    # Protein name
    prot = ""
    try:
        prot = data.get("proteinDescription", {}) \
                   .get("recommendedName", {}) \
                   .get("fullName", {}) \
                   .get("value", "")
    except Exception:
        log.debug("No recommended protein name found")
    # PDB cross-references
    pdb_ids = [
        x.get("id")
        for x in data.get("uniProtKBCrossReferences", [])
        if x.get("database") == "PDB" and x.get("id")
    ]
    # Natural variant (SNP) features
    variants_list = []
    for feat in data.get("features", []):
        if feat.get("type") == "Natural variant":
            vid  = feat.get("featureId", "")
            loc  = feat.get("location", {}).get("start", {}).get("value", "")
            alt  = feat.get("alternativeSequence", {})
            orig = alt.get("originalSequence", "")
            new  = ""
            if alt.get("alternativeSequences"):
                new = alt["alternativeSequences"][0]
            if vid and loc:
                variants_list.append(f"{vid}({loc}{orig}>{new})")
    variants = "; ".join(variants_list) if variants_list else ""
    # Comments
    pharm        = get_first_comment(data.get("comments", []), "PHARMACEUTICAL")
    tissue       = get_first_comment(data.get("comments", []), "TISSUE SPECIFICITY")
    cofactor     = get_first_comment(data.get("comments", []), "COFACTOR")
    polymorphism = get_first_comment(data.get("comments", []), "POLYMORPHISM")
    # Cross-refs TCDB / DrugCentral
    tcdb_id, drugcentral_id = "", ""
    for x in data.get("uniProtKBCrossReferences", []):
        if x.get("database") == "TCDB":
            tcdb_id = x.get("id", "")
        if x.get("database") == "DrugCentral":
            drugcentral_id = x.get("id", "")
    # Annotation score & sequence
    annotation_score = data.get("annotationScore", "")
    sequence         = data.get("sequence", {}).get("value", "")
    # Add interactions
    related = find_related_uniprots(
        data.get("primaryAccession", ""), 
        data
    )
    related_ids = ""
    related_genes = ""
    if related:
        # If partners found, concatenate gene names
        related_genes = "-".join([gene] + [r["gene_name"] for r in related])
        related_ids = ";".join([r["uniprot_id"] for r in related])
    else:
        related_genes = gene
    return {
        "gene": related_genes,
        "protein_name": prot,
        "pdb_ids": pdb_ids,
        "variants": variants,
        "pharmaceutical": pharm,
        "tcdb_id": tcdb_id,
        "drugcentral_id": drugcentral_id,
        "tissue": tissue,
        "annotation_score": annotation_score,
        "cofactor": cofactor,
        "polymorphism": polymorphism,
        "sequence": sequence,
        "interaction_partners": related_ids,
    }
