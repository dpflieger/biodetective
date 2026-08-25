#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_genome_stats.py — Taille du génome et nombre de chromosomes par espèce.

Source : les rapports NCBI GENOME_REPORTS, bien plus légers que d'interroger
les assemblages un par un.

    python3 build_genome_stats.py            # utilise genome_reports/
    python3 build_genome_stats.py --download # récupère les fichiers d'abord

Produit genome_stats.json, relu par data_pipeline.py. Les rapports bruts
(plus de 200 Mo) ne sont pas versionnés ; le JSON, compact, l'est.

Une espèce a souvent des dizaines d'assemblages. On retient en priorité celui
dont les réplicons portent des accessions RefSeq (NC_…) : c'est le génome de
référence, et c'est aussi lui qu'on retrouvera dans nt.
"""

import argparse
import json
import os
import re
import sqlite3
import sys
import urllib.request

REPORTS_DIR = "genome_reports"
FILES = ["eukaryotes.txt", "prokaryotes.txt", "viruses.txt"]
BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS"
OUT = "genome_stats.json"

# « chromosome 1:NC_003070.9/CP002684.1; chromosome 2:… »
REPLICON_RE = re.compile(r"chromosome\s+([^:;]+):([^;/]+)", re.I)


def download():
    os.makedirs(REPORTS_DIR, exist_ok=True)
    for f in FILES:
        dest = os.path.join(REPORTS_DIR, f)
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"  {f} déjà présent")
            continue
        print(f"  téléchargement de {f}…", flush=True)
        urllib.request.urlretrieve(f"{BASE}/{f}", dest)


def parse_replicons(field):
    """[(nom du chromosome, accession)] — sert aussi à situer un hit."""
    out = []
    for name, acc in REPLICON_RE.findall(field or ""):
        name = name.strip()
        acc = acc.strip()
        if acc and acc.lower() not in ("-", "na"):
            out.append((name, acc))
    return out


def score(chroms, size, has_refseq, has_accession):
    """Plus c'est haut, meilleur c'est.

    L'accession d'abord : certaines lignes n'en portent pas et donnent des
    tailles fantaisistes — 11 Mb pour Escherichia coli, qui en fait 4,6.
    Puis RefSeq, qui désigne le génome de référence.
    """
    return (1 if has_accession else 0,
            1 if has_refseq else 0,
            len(chroms), size or 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--db", default="biodetective.db")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    if args.download:
        download()

    conn = sqlite3.connect(args.db)
    wanted = {str(r[0]) for r in
              conn.execute("SELECT taxonomy_id FROM organisms")}
    conn.close()
    print(f"{len(wanted):,} espèces à renseigner")

    best = {}
    seen_rows = 0
    for f in FILES:
        path = os.path.join(REPORTS_DIR, f)
        if not os.path.exists(path):
            print(f"  {path} absent, ignoré")
            continue
        with open(path, encoding="utf-8", errors="replace") as fh:
            header = next(fh, "").lstrip("#").rstrip("\n").split("\t")
            # Les trois rapports n'ont PAS la même disposition : dans
            # prokaryotes.txt l'accession est en 19e colonne et les réplicons
            # en 9e, l'inverse d'eukaryotes.txt. Et viruses.txt donne des Kb.
            # On repère donc les colonnes par leur intitulé.
            idx = {name.strip().lower(): i for i, name in enumerate(header)}
            i_tax = idx.get("taxid")
            i_rep = idx.get("replicons")
            i_asm = idx.get("assembly accession")
            i_mb = idx.get("size (mb)")
            i_kb = idx.get("size (kb)")
            if i_tax is None:
                print(f"  {f} : colonne TaxID introuvable, ignoré")
                continue
            print(f"  {f} : taille={'Mb' if i_mb is not None else 'Kb'}, "
                  f"réplicons={'oui' if i_rep is not None else 'non'}")

            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) <= i_tax:
                    continue
                taxid = c[i_tax].strip()
                if taxid not in wanted:
                    continue
                seen_rows += 1
                size = None
                if i_mb is not None and len(c) > i_mb:
                    try:
                        size = float(c[i_mb])
                    except ValueError:
                        pass
                elif i_kb is not None and len(c) > i_kb:
                    try:
                        size = float(c[i_kb]) / 1000.0
                    except ValueError:
                        pass
                chroms = parse_replicons(
                    c[i_rep] if i_rep is not None and len(c) > i_rep else "")
                accession = (c[i_asm].strip()
                             if i_asm is not None and len(c) > i_asm else "")
                has_accession = bool(accession) and accession != "-"
                has_refseq = any(a.startswith(("NC_", "NZ_"))
                                 for _n, a in chroms)
                s = score(chroms, size, has_refseq, has_accession)
                if taxid not in best or s > best[taxid][0]:
                    best[taxid] = (s, {
                        "genome_size_mb": round(size, 2) if size else None,
                        "chromosomes": len(chroms),
                        "assembly": accession if has_accession else None,
                        # nom -> accession : permettra de situer un hit sur
                        # le bon chromosome, et de le dessiner un jour.
                        "replicons": {n: a for n, a in chroms} or None,
                    })

    stats = {t: v for t, (_s, v) in best.items()}
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, ensure_ascii=False, sort_keys=True)
        fh.write("\n")

    with_size = sum(1 for v in stats.values() if v["genome_size_mb"])
    with_chr = sum(1 for v in stats.values() if v["chromosomes"])
    print(f"\n{seen_rows:,} lignes retenues -> {len(stats):,} espèces")
    print(f"  taille de génome : {with_size:,} "
          f"({100*with_size/max(len(wanted),1):.1f} % du total)")
    print(f"  chromosomes      : {with_chr:,}")
    print(f"  {args.out} : {os.path.getsize(args.out)/1e6:.1f} Mo")
    return 0


if __name__ == "__main__":
    sys.exit(main())
