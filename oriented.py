#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Module d'alignement multi-référence pour structures PDB
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

# --- Constantes & expressions régulières ---
DUM_RE = re.compile(r"\bDUM\b", re.IGNORECASE)
HYP_RE = re.compile(
    r'=?HYPERLINK\(\s*"[^"]*"\s*[;,]\s*"(?P<id>[A-Za-z0-9\-]+)"\s*\)',
    re.IGNORECASE
)

# --- Configuration du logging ---
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Fonctions utilitaires pour l'extraction et le traitement des données
# ---------------------------------------------------------------------------

def extract_pdb_from_cell(cell) -> str:
    """
    Extrait l'identifiant PDB d'une cellule Excel, qu'il soit sous forme de lien hypertexte ou de texte.
    
    Parameters
    ----------
    cell : openpyxl.cell.Cell
        Cellule Excel à traiter
        
    Returns
    -------
    str
        Identifiant PDB extrait, ou None si la cellule est vide
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
    Charge une feuille spécifique d'un fichier Excel et la convertit en DataFrame.
    
    Parameters
    ----------
    excel_file : Path
        Chemin vers le fichier Excel
    name : str
        Nom de la feuille à charger
        
    Returns
    -------
    pd.DataFrame
        DataFrame contenant les données de la feuille
        
    Raises
    ------
    SystemExit
        Si le fichier Excel ou la feuille n'existe pas
    """
    if not excel_file.is_file():
        logger.error(f"Excel manquant: {excel_file}")
        sys.exit(1)
    wb = openpyxl.load_workbook(excel_file, data_only=False)
    if name not in wb.sheetnames:
        logger.error(f"Feuille '{name}' manquante ({wb.sheetnames})")
        sys.exit(1)
    ws = wb[name]
    headers = [str(c.value).strip().lower().replace(" ", "_") for c in ws[1]]
    rows = []
    for row in ws.iter_rows(min_row=2):
        rec = {h: extract_pdb_from_cell(c) for h, c in zip(headers, row)}
        rows.append(rec)
    df = pd.DataFrame(rows)
    logger.info(f"[{name}] chargée: {df.shape[0]} lignes")
    return df


def filter_and_cascade(df_det: pd.DataFrame, df_ali: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filtre les structures et cascade le filtrage aux données d'alignement.
    
    Parameters
    ----------
    df_det : pd.DataFrame
        DataFrame contenant les détails des structures PDB
    df_ali : pd.DataFrame
        DataFrame contenant les informations d'alignement
        
    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Tuple contenant les DataFrames filtrés (détails, alignement)
    """
    # Cascade vers Alignment_Summary
    valid = set(df_det['pdb_id'].dropna().str.upper())
    df_ali_proc = df_ali[df_ali['pdb_id'].str.upper().isin(valid)].reset_index(drop=True)
    logger.info(f"Traité: PDB_Details {df_det.shape[0]} lignes, Alignment_Summary {df_ali_proc.shape[0]} lignes")
    return df_det, df_ali_proc


def superimpose_structures(ref_path: Path,
                           df_targets: pd.DataFrame,
                           out_root: Path) -> list[tuple[str, str, float, float]]:
    """
    Superpose les structures cibles sur une structure de référence.
    
    Parameters
    ----------
    ref_path : Path
        Chemin vers la structure de référence
    df_targets : pd.DataFrame
        DataFrame contenant les informations sur les structures cibles
    out_root : Path
        Répertoire racine pour les structures alignées
        
    Returns
    -------
    list[tuple[str, str, float, float]]
        Liste de tuples (pdb_id, chemin_aligné, rmsd_avant, rmsd_après)
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
        # Prépare le chemin de sortie
        subdir = out_root / f"{gene}_{uniprot}"
        subdir.mkdir(parents=True, exist_ok=True)
        out_name = f"oriented_{gene}_{uniprot}_{pdb_id}.pdb"
        out_path = subdir / out_name
        # Charge et aligne
        tgt_struct = parser.get_structure(pdb_id, str(in_path))
        tgt_atoms = [a for a in tgt_struct[0].get_atoms() if a.get_id() == 'CA']
        n = min(len(ref_atoms), len(tgt_atoms))
        ref_coords = np.array([a.get_coord() for a in ref_atoms[:n]])
        tgt_coords = np.array([a.get_coord() for a in tgt_atoms[:n]])
        # RMSD avant
        rms_before = np.sqrt(((ref_coords - tgt_coords)**2).sum(axis=1).mean())
        # Superpose et RMSD après
        sup.set_atoms(ref_atoms[:n], tgt_atoms[:n])
        rms_after = sup.rms
        sup.apply(tgt_struct.get_atoms())
        io.set_structure(tgt_struct)
        io.save(str(out_path))
        logger.info(f"{pdb_id}: avant={rms_before:.3f}Å, après={rms_after:.3f}Å -> {out_path}")
        results.append((pdb_id, str(out_path), rms_before, rms_after))
    return results


# ---------------------------------------------------------------------------
# Nouvelles fonctions pour l'alignement multi-référence
# ---------------------------------------------------------------------------

def group_structures_by_attribute(df_det: pd.DataFrame, attribute: str) -> Dict[str, List[str]]:
    """
    Regroupe les structures PDB par attribut spécifié (gene_name ou uniprot_id).
    
    Parameters
    ----------
    df_det : pd.DataFrame
        DataFrame contenant les détails des structures PDB
    attribute : str
        Attribut pour le regroupement ('gene_name' ou 'uniprot_id')
        
    Returns
    -------
    Dict[str, List[str]]
        Dictionnaire {valeur_attribut: [liste_pdb_ids]}
    """
    if attribute not in df_det.columns:
        logger.error(f"Attribut '{attribute}' non trouvé dans le DataFrame")
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
    Trouve une référence contenant DUM pour chaque groupe.
    
    Parameters
    ----------
    df_det : pd.DataFrame
        DataFrame contenant les détails des structures PDB
    groups : Dict[str, List[str]]
        Dictionnaire des groupes de structures
        
    Returns
    -------
    Dict[str, str]
        Dictionnaire {groupe: pdb_id_reference}
    """
    references = {}
    
    # Vérifie si 'small_molecules' est dans les colonnes
    if 'small_molecules' not in df_det.columns:
        logger.warning("Colonne 'small_molecules' non trouvée, impossible de détecter DUM")
        return references
    
    # Pour chaque groupe, trouve une référence avec DUM
    for group, pdb_ids in groups.items():
        # Filtre le DataFrame pour ce groupe
        group_df = df_det[df_det['pdb_id'].isin(pdb_ids)].copy()
        
        # Cherche les structures avec DUM
        sm = group_df['small_molecules'].fillna("").str.upper()
        dum_mask = sm.str.contains(DUM_RE)
        
        if dum_mask.any():
            # Prend la première structure avec DUM comme référence
            ref_id = group_df.loc[dum_mask, 'pdb_id'].iloc[0]
            references[group] = ref_id
            logger.info(f"Groupe {group}: référence avec DUM trouvée: {ref_id}")
        else:
            logger.warning(f"Groupe {group}: aucune référence avec DUM trouvée")
    
    return references


def align_structures_by_groups(excel_file: Path, 
                              group_by: str = 'gene_name',
                              out_root: Optional[Path] = None) -> Dict[str, List[tuple]]:
    """
    Aligne les structures PDB par groupes, en utilisant une référence avec DUM pour chaque groupe.
    
    Parameters
    ----------
    excel_file : Path
        Chemin vers le fichier Excel contenant les données
    group_by : str, optional
        Attribut pour le regroupement ('gene_name' ou 'uniprot_id'), par défaut 'gene_name'
    out_root : Path, optional
        Répertoire racine pour les structures alignées, par défaut 'pdb_structures_oriented'
        
    Returns
    -------
    Dict[str, List[tuple]]
        Dictionnaire {groupe: [(pdb_id, chemin_aligné, rmsd_avant, rmsd_après)]}
    """
    if out_root is None:
        out_root = Path('pdb_structures_oriented')
    
    # Charge les feuilles Excel
    df_det_all = load_sheet(excel_file, 'PDB_Details')
    df_ali_all = load_sheet(excel_file, 'Alignment_Summary')
    
    # Vérifie que les colonnes nécessaires existent
    required_cols = ['pdb_id', group_by, 'small_molecules']
    for col in required_cols:
        if col not in df_det_all.columns:
            logger.error(f"Colonne requise '{col}' non trouvée dans PDB_Details")
            return {}
    
    if 'patched_pdb' not in df_ali_all.columns:
        logger.error("Colonne 'patched_pdb' non trouvée dans Alignment_Summary")
        return {}
    
    # Groupe les structures par attribut spécifié
    groups = group_structures_by_attribute(df_det_all, group_by)
    logger.info(f"Structures regroupées par {group_by}: {len(groups)} groupes trouvés")
    
    # Trouve les références avec DUM pour chaque groupe
    references = find_dum_references(df_det_all, groups)
    logger.info(f"Références avec DUM trouvées pour {len(references)}/{len(groups)} groupes")
    
    # Prépare le dictionnaire de résultats
    all_results = {}
    all_aligned_results = []
    
    # Pour chaque groupe avec une référence
    for group, ref_id in references.items():
        logger.info(f"Traitement du groupe {group} avec référence {ref_id}")
        
        # Obtient le chemin de la référence
        ref_row = df_ali_all[df_ali_all['pdb_id'] == ref_id]
        if ref_row.empty:
            logger.warning(f"Référence {ref_id} non trouvée dans Alignment_Summary, groupe {group} ignoré")
            continue
        
        ref_path = Path(ref_row.iloc[0]['patched_pdb'])
        if not ref_path.is_file():
            logger.error(f"Fichier de référence non trouvé: {ref_path}, groupe {group} ignoré")
            continue
        
        # Filtre les structures du groupe (sauf la référence)
        group_pdb_ids = groups[group]
        target_pdb_ids = [pid for pid in group_pdb_ids if pid != ref_id]
        
        # Prépare le DataFrame des cibles
        df_targets = pd.merge(
            df_ali_all[df_ali_all['pdb_id'].isin(target_pdb_ids)],
            df_det_all[['pdb_id', 'gene_name', 'uniprot_id']],
            on='pdb_id', how='left'
        )
        
        if df_targets.empty:
            logger.warning(f"Aucune structure cible trouvée pour le groupe {group}")
            continue
        
        # Aligne les structures sur la référence
        results = superimpose_structures(ref_path, df_targets, out_root)
        
        # Copie la référence (déjà orientée)
        ref_det_row = df_det_all[df_det_all['pdb_id'] == ref_id].iloc[0]
        gene = ref_det_row['gene_name']
        uniprot = ref_det_row['uniprot_id']
        subdir = out_root / f"{gene}_{uniprot}"
        subdir.mkdir(parents=True, exist_ok=True)
        out_name = f"oriented_{gene}_{uniprot}_{ref_id}.pdb"
        out_path = subdir / out_name
        shutil.copy(str(ref_path), str(out_path))
        logger.info(f"Copié référence {ref_id} -> {out_path}")
        
        # Ajoute la référence aux résultats
        results.append((ref_id, str(out_path), None, None))
        
        # Stocke les résultats pour ce groupe
        all_results[group] = results
        all_aligned_results.extend(results)
    
    # Sauvegarde tous les résultats dans Excel
    df_res = pd.DataFrame(all_aligned_results, columns=['pdb_id', 'aligned_path', 'rmsd_before', 'rmsd_after'])
    with pd.ExcelWriter(excel_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        df_res.to_excel(writer, sheet_name='Alignment_RMSD_MultiRef', index=False)
    logger.info(f"Résultats sauvegardés dans 'Alignment_RMSD_MultiRef' dans {excel_file}")
    
    return all_results


# ---------------------------------------------------------------------------
# Fonction principale pour l'utilisation en ligne de commande
# ---------------------------------------------------------------------------

def main():
    """
    Fonction principale pour l'utilisation en ligne de commande.
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Aligne et oriente les fichiers PDB par groupes sur des références contenant DUM"
    )
    parser.add_argument('-e', '--excel', required=True,
                        help='Fichier Excel contenant les données')
    parser.add_argument('-g', '--group-by', choices=['gene_name', 'uniprot_id'], default='gene_name',
                        help='Attribut pour le regroupement (gene_name ou uniprot_id)')
    parser.add_argument('-o', '--output', default='pdb_structures_oriented',
                        help='Répertoire racine pour les structures alignées')
    
    args = parser.parse_args()
    
    # Configure le logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    
    # Exécute l'alignement
    excel_file = Path(args.excel)
    out_root = Path(args.output)
    
    align_structures_by_groups(excel_file, args.group_by, out_root)


if __name__ == '__main__':
    main()
