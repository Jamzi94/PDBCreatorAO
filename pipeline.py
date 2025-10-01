import os, sys, argparse, logging
import pandas as pd
from openpyxl import load_workbook
from openpyxl.workbook import Workbook
from openpyxl.styles import PatternFill, Font
from uniprot import fetch_uniprot_json, parse_uniprot_data, find_related_uniprots
from alphafold import (get_alphafold_pdb_info, download_alphafold_pdb,
                       get_alphafold_details)
from pdb_download import download_pdb_structure
from pdb_details import get_pdb_details
from structure_analysis import (detect_gaps_and_discard_coils,
                                parse_small_molecules,
                                parse_helix_sheet_annotations)
from alignment import align_and_patch
from oriented import align_structures_by_groups
from pathlib import Path


# --------------------------------------------------------------------------- #

def setup_logging(logfile: str, verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(logfile, "w", "utf-8"),
                  logging.StreamHandler(sys.stdout)]
    )

# --------------------------------------------------------------------------- #

def style_link(cell, link_font):
    if cell.value and str(cell.value).startswith("=HYPERLINK"):
        cell.font = link_font

# --------------------------------------------------------------------------- #
def _add_af_gap_row(ws, upid, gene, prot, pdb_id,
                    start, end, ss_map, link_font):
    ws.append([
        f'=HYPERLINK("https://www.uniprot.org/uniprotkb/{upid}/entry","{upid}")',
        gene, prot,
        f'=HYPERLINK("https://www.rcsb.org/structure/{pdb_id}", "{pdb_id}")',
        "GAP_AF", f"{start}–{end}", end - start + 1,
        " ".join(ss_map.get(r, "C") for r in range(start, end + 1))
    ])
    for col in (1, 3):
        style_link(ws.cell(ws.max_row, col), link_font)

# --------------------------------------------------------------------------- #
def create_excel_sheets(path: str):
    wb = Workbook()
    ws1 = wb.active; ws1.title = "UniProt Results"
    ws1.append([
        "UniProt ID", "Gene Name", "Protein Name",
        "PDB Entries", "AlphaFoldDB Entries",
        "SNPs", "Pharmaceutical Use", "TCDB ID", "DrugCentral ID",
        "Tissue Specificity", "Annotation Score", "Cofactors",
        "Polymorphism", "FASTA Sequence"
    ])

    ws2 = wb.create_sheet("PDB_Details")
    ws2.append([
        "UniProt ID", "Gene Name", "Protein Name", "PDB ID",
        "Experimental Method", "Resolution (Å)", "Determination Methodology",
        "Conformation", "Environment", "Polymer Entity Count",
        "Entity Descriptions", "Chains", "Entity Organisms",
        "PubMed ID", "Small Molecules"
    ])

    ws3 = wb.create_sheet("Structural_Integrity")
    ws3.append([
        "UniProt ID", "Gene Name", "Protein Name", "PDB ID",
        "Type", "Residus", "Length", "Secondary Structure"
    ])

    wb.save(path)
    return wb, ws1, ws2, ws3

# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser("Protein structure pipeline")
    ap.add_argument("--base", default=os.getcwd())
    ap.add_argument("--input", default="uniprot_template.xlsx")
    ap.add_argument("--output", default="uniprot_results.xlsx")
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("--ids", nargs="+")
    args = ap.parse_args()

    BASE = os.path.abspath(args.base)
    input_xlsx = os.path.abspath(args.input)
    output_xlsx = os.path.join(BASE, args.output)
    EXP_DIR = os.path.join(BASE, "pdb_structures_exp")
    PRED_DIR = os.path.join(BASE, "pdb_structures_pred")
    ALIGNED_DIR = os.path.join(BASE, "pdb_structures_af_aligned")
    for d in (EXP_DIR, PRED_DIR, ALIGNED_DIR):
        os.makedirs(d, exist_ok=True)

    setup_logging(os.path.join(BASE, "pipeline.log"), args.verbose)
    log = logging.getLogger(__name__)
    log.info("Starting pipeline")

    # ---------- read UniProt IDs ------------------------------------------
    if args.ids:
        uniprot_ids = [uid.strip() for uid in args.ids]
    else:
        if not os.path.isfile(input_xlsx):
            Workbook().save(input_xlsx)
            log.info(f"Template '{input_xlsx}' created – fill and rerun.")
            sys.exit()
        ws_in = load_workbook(input_xlsx).active
        uniprot_ids = [str(r[0]).strip()
                       for r in ws_in.iter_rows(min_row=2, values_only=True)
                       if r[0]]

    wb_out, ws_uni, ws_det, ws_struct = create_excel_sheets(output_xlsx)
    link_font = Font(color="0000FF", underline="single")

    # prot_len map per folder "Gene_UPID"
    prot_len_map = {}
    partner_ids_map = {}

    # ---------- UniProt & experimental PDB part ---------------------------
    for upid in uniprot_ids:
        if not upid:
            continue
        log.info(f"Processing {upid}")
        try:
            up_json = fetch_uniprot_json(upid)
        except Exception as e:
            log.error(f"UniProt fetch failed for {upid}: {e}")
            continue

        info = parse_uniprot_data(up_json)
        partners = find_related_uniprots(upid, up_json)
        info["interaction_partners"] = [p["uniprot_id"] for p in partners]
        partner_ids_map[f"{info['gene']}_{upid}"] = info["interaction_partners"]
        gene = info["gene"] or "UnknownGene"
        prot = info["protein_name"]
        prot_len_map[f"{gene}_{upid}"] = len(info["sequence"])
        exp_ids = info["pdb_ids"]
        pred_id = get_alphafold_pdb_info(upid)
        pred_ids = [pred_id] if pred_id else []

        # --- write UniProt sheet ------------------------------------------
        ws_uni.append([
            f'=HYPERLINK("https://www.uniprot.org/uniprotkb/{upid}/entry","{upid}")',
            f'=HYPERLINK("https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}", "{gene}")',
            prot, ", ".join(exp_ids), ", ".join(pred_ids),
            info["variants"], info["pharmaceutical"], info["tcdb_id"],
            info["drugcentral_id"], info["tissue"], info["annotation_score"],
            info["cofactor"], info["polymorphism"], info["sequence"]
        ])
        for col in (1, 2):
            style_link(ws_uni.cell(ws_uni.max_row, col), link_font)

        # --- download experimental PDBs, fill Structural_Integrity --------
        for pdb_id in exp_ids:
            out_dir = os.path.join(EXP_DIR, f"{gene}_{upid}")
            pdb_path = download_pdb_structure(
                pdb_id, out_dir, f"{gene}_{upid}_{pdb_id}.pdb"
            )
            if not pdb_path:
                continue
            issues = detect_gaps_and_discard_coils(
                pdb_path, chain_id="A", plddt_threshold=60.0,
                coil_length_threshold=4
            )
            ss_map = parse_helix_sheet_annotations(pdb_path, chain_id="A")
            for typ, start, end in issues:
                ws_struct.append([
                    f'=HYPERLINK("https://www.uniprot.org/uniprotkb/{upid}/entry","{upid}")',
                    gene, prot,
                    f'=HYPERLINK("https://www.rcsb.org/structure/{pdb_id}", "{pdb_id}")',
                    typ, f"{start}–{end}", end - start + 1,
                    " ".join(ss_map.get(r, "C") for r in range(start, end + 1))
                ])
                for c in (1, 4):
                    style_link(ws_struct.cell(ws_struct.max_row, c), link_font)

        # --- download AlphaFold PDB ---------------------------------------
        if pred_ids:
            pred_dir = os.path.join(PRED_DIR, f"{gene}_{upid}")
            os.makedirs(pred_dir, exist_ok=True)

            download_alphafold_pdb(
                upid,
                os.path.join(pred_dir, f"{gene}_{upid}_{pred_id}.pdb")
            )

            partners = info["interaction_partners"]
            if isinstance(partners, str):
                partner_ids = [p.strip() for p in partners.split(",") if p.strip()]
            else:
                partner_ids = partners

            log.info(f"▶ DEBUG: gene={gene!r}, upid={upid!r}, pred_id={pred_id!r}")
            log.info(f"▶ DEBUG: raw partners={partners!r}")
            log.info(f"▶ DEBUG: parsed partner_ids={partner_ids!r}")
            log.info(f"▶ DEBUG: pred_dir={pred_dir!r}")

            for partner_id in partner_ids:
                if not partner_id:
                    continue
                partner_json = fetch_uniprot_json(partner_id)
                partner_info = parse_uniprot_data(partner_json)
                partner_gene = partner_info["gene"] or partner_id
                prot_len_map[f"{partner_gene}_{partner_id}"] = len(partner_info["sequence"])
                partner_pred_id = get_alphafold_pdb_info(partner_id)
                if partner_pred_id:
                    download_alphafold_pdb(
                        partner_id,
                        os.path.join(pred_dir, f"{gene}_{partner_id}_{partner_pred_id}.pdb")
                    )

        # --- write PDB_Details (experimental) ----------------------------
        for pdb_id in exp_ids:
            det = get_pdb_details(pdb_id)
            pdb_file = os.path.join(EXP_DIR, f"{gene}_{upid}",
                                    f"{gene}_{upid}_{pdb_id}.pdb")
            smalls = parse_small_molecules(pdb_file) if os.path.isfile(pdb_file) else []
            pmid = det.get("pubmed_id", "")
            pmcell = (f'=HYPERLINK("https://pubmed.ncbi.nlm.nih.gov/{pmid}", "{pmid}")'
                      if pmid else "")
            ws_det.append([
                f'=HYPERLINK("https://www.uniprot.org/uniprotkb/{upid}/entry","{upid}")',
                gene, prot,
                f'=HYPERLINK("https://www.rcsb.org/structure/{pdb_id}", "{pdb_id}")',
                det["method"], det["resolution"], det["methodology"],
                det["conformation"], det["environment"], det["count"],
                det["descs"], det["chains"], det["orgs"], pmcell,
                ", ".join(smalls)
            ])
            for col in (1, 4, 14):
                style_link(ws_det.cell(ws_det.max_row, col), link_font)

        # --- PDB details for AlphaFold -----------------------------------
        if pred_ids:
            pred_det = get_alphafold_details(upid, prot)
            ws_det.append([
                f'=HYPERLINK("https://www.uniprot.org/uniprotkb/{upid}/entry","{upid}")',
                gene, prot,
                f'=HYPERLINK("https://alphafold.ebi.ac.uk/entry/{upid}", "{pred_id}")',
                pred_det["method"], pred_det["resolution"], pred_det["methodology"],
                pred_det["conformation"], pred_det["environment"], pred_det["count"],
                pred_det["descs"], pred_det["chains"], pred_det["orgs"], "", ""
            ])
            for col in (1, 4):
                style_link(ws_det.cell(ws_det.max_row, col), link_font)

    wb_out.save(output_xlsx)

    # ---------- Alignment & patching -------------------------------------
    summary = []
    for gene_dir in os.listdir(PRED_DIR):
        af_dir = os.path.join(PRED_DIR, gene_dir)
        if not os.path.isdir(af_dir):
            continue

        gene_name, upid = gene_dir.split("_", 1)
        prot_len = prot_len_map.get(gene_dir)
        if prot_len is None:
            log.warning(f"No prot_len for {gene_dir} – skipped.")
            continue

        exp_dir = os.path.join(EXP_DIR, gene_dir)
        if not os.path.isdir(exp_dir):
            continue
        partner_ids = partner_ids_map.get(gene_dir, [])
        for exp_file in os.listdir(exp_dir):
            if not exp_file.endswith(".pdb"):
                continue
            exp_path = os.path.join(exp_dir, exp_file)
            partner_ids = partner_ids_map.get(gene_dir, [])
            res = align_and_patch(
                af_path=af_dir,
                exp_path=exp_path,
                gene=gene_name,
                prot_len=prot_len,
                window=args.window,
                outdir=ALIGNED_DIR,
                partner_ids=partner_ids
            )
            summary.append(res)

    if summary:
        df = pd.DataFrame(summary)
        wb = load_workbook(output_xlsx)
        if "Alignment_Summary" in wb.sheetnames:
            del wb["Alignment_Summary"]
            wb.save(output_xlsx)

        with pd.ExcelWriter(output_xlsx, mode="a",
                            engine="openpyxl", if_sheet_exists="overlay") as wr:
            df.to_excel(wr, sheet_name="Alignment_Summary", index=False)

        wb = load_workbook(output_xlsx)
        ws = wb["Alignment_Summary"]
        col_stat = [i for i, c in enumerate(ws[1], 1) if c.value == "status"][0]
        for row in ws.iter_rows(min_row=2, min_col=col_stat, max_col=col_stat):
            cell = row[0]
            cell.fill = PatternFill("solid",
                                    fgColor="C6EFCE" if cell.value == "ok" else "FFC7CE")
        
        col_pdb = [i for i, c in enumerate(ws[1], 1) if c.value == "pdb_id"][0]
        for row_idx in range(2, ws.max_row + 1):
            style_link(ws.cell(row_idx, col_pdb), link_font)
            wb.save(output_xlsx)

    # ---------- Multi-reference orientation (align_and_patch) ---------------
    log.info("Running multi-reference alignment and orientation")
    try:
        align_structures_by_groups(
            excel_file=Path(output_xlsx),  # <-- convert to Path
            group_by='gene_name',
            out_root=Path(os.path.join(BASE, "pdb_structures_oriented"))
        )
        log.info("Multi-reference alignment completed")
    except Exception as e:
        log.error(f"Multi-reference alignment failed: {e}")


        log.info("Pipeline completed")

if __name__ == "__main__":
    main()
