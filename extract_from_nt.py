#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_from_nt.py — À LANCER SUR LE SERVEUR QUI HÉBERGE nt.

Extrait de nt une séquence représentative par espèce, pour les taxons de
BioDetective qui possèdent une image. Produit un FASTA d'une trentaine de
mégaoctets, à rapatrier sur la machine de démonstration.

    python3 extract_from_nt.py --db /chemin/vers/nt --taxids taxids.txt

Ne nécessite que blastdbcmd et python3. Rien d'autre du projet.

⚠️ La TAILLE de la base finale est le paramètre critique. La E-value est
proportionnelle à la taille de la banque : les briques du Brickopore font
16 bases, ce qui donne E = 0,01 contre 16 Mpb, mais E = 1 contre 1,6 Gpb,
c'est-à-dire plus rien d'exploitable. D'où une seule séquence par espèce et
un plafond de longueur : surtout pas les génomes complets.
"""

import argparse
import re
import subprocess
import sys
from collections import defaultdict

# Fourchette utile : un marqueur, pas un chromosome.
MIN_LEN, MAX_LEN = 400, 3000

# Marqueurs classiques de code-barres ADN, par ordre de préférence.
PREFERRED = [
    r"cytochrome c oxidase subunit (?:1|I)\b", r"\bCOI\b", r"\bCOX1\b",
    r"ribulose.*carboxylase", r"\brbcL\b", r"\bmatK\b",
    r"internal transcribed spacer", r"\bITS\b",
    r"16S ribosomal RNA", r"18S ribosomal RNA", r"28S ribosomal RNA",
    r"cytochrome b\b", r"\bcytb\b",
]
PREFERRED_RE = [re.compile(p, re.I) for p in PREFERRED]

# À écarter : tout ce qui est génomique en vrac.
REJECT_RE = re.compile(
    r"complete genome|chromosome|whole genome shotgun|scaffold|contig|"
    r"unplaced|patch|assembly|clone library|predicted|hypothetical",
    re.I)


def rank(title, length):
    """Note une séquence : plus c'est bas, mieux c'est."""
    if REJECT_RE.search(title):
        return (9, 0)
    for i, rx in enumerate(PREFERRED_RE):
        if rx.search(title):
            # à marqueur égal, on préfère une longueur médiane
            return (i, abs(length - 800))
    return (len(PREFERRED_RE), abs(length - 800))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="chemin de la banque nt")
    ap.add_argument("--taxids", default="taxids.txt")
    ap.add_argument("--out", default="biodetective_subset.fasta")
    ap.add_argument("--per-species", type=int, default=1)
    args = ap.parse_args()

    wanted = {l.strip() for l in open(args.taxids) if l.strip()}
    print(f"{len(wanted)} taxons demandés", flush=True)

    # 1) Catalogue : accession, taxid, longueur, titre. Sortie texte, rapide.
    print("Inventaire des séquences disponibles (peut prendre plusieurs "
          "minutes)…", flush=True)
    cmd = [ "blastdbcmd", "-db", args.db, "-taxidlist", args.taxids,
            "-outfmt", "%a\t%T\t%l\t%t" ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, errors="replace")

    best = defaultdict(list)
    seen = 0
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t", 3)
        if len(parts) < 4:
            continue
        acc, taxid, length, title = parts
        seen += 1
        if seen % 500000 == 0:
            print(f"  {seen:,} séquences examinées, "
                  f"{len(best):,} espèces couvertes", flush=True)
        try:
            length = int(length)
        except ValueError:
            continue
        if not (MIN_LEN <= length <= MAX_LEN):
            continue
        if taxid not in wanted:
            continue
        r = rank(title, length)
        if r[0] >= 9:
            continue
        best[taxid].append((r, acc, length, title))

    proc.wait()
    err = proc.stderr.read()
    if proc.returncode != 0:
        print(f"blastdbcmd a échoué : {err[:500]}", file=sys.stderr)
        return 1

    print(f"{seen:,} séquences examinées, {len(best):,} espèces retenues",
          flush=True)

    chosen = []
    for taxid, cands in best.items():
        cands.sort()
        for _r, acc, _l, _t in cands[:args.per_species]:
            chosen.append((acc, taxid))

    with open("chosen_accessions.txt", "w") as fh:
        for acc, _ in chosen:
            fh.write(acc + "\n")
    print(f"{len(chosen)} séquences choisies", flush=True)

    # 2) Extraction des seules séquences retenues.
    print("Extraction du FASTA…", flush=True)
    raw = subprocess.run(
        ["blastdbcmd", "-db", args.db, "-entry_batch", "chosen_accessions.txt",
         "-outfmt", "%a\t%T\t%s"],
        capture_output=True, text=True, errors="replace")
    if raw.returncode != 0:
        print(f"échec : {raw.stderr[:500]}", file=sys.stderr)
        return 1

    # Le taxid est écrit dans l'identifiant : la machine de démonstration n'a
    # alors besoin d'aucun fichier de taxonomie pour relier un hit à sa fiche.
    n = 0
    total = 0
    with open(args.out, "w") as fh:
        for line in raw.stdout.splitlines():
            p = line.split("\t")
            if len(p) < 3:
                continue
            acc, taxid, seq = p[0], p[1], p[2]
            fh.write(f">{taxid}|{acc}\n{seq}\n")
            n += 1
            total += len(seq)

    print("-" * 55)
    print(f"{n} séquences écrites dans {args.out}")
    print(f"{total/1e6:.1f} Mpb  ->  E-value attendue pour 16 pb : "
          f"~{0.010 * (total/16.2e6):.3f}")
    print("\nÀ rapatrier sur la machine de démonstration, puis :")
    print(f"    makeblastdb -in {args.out} -dbtype nucl -out biodetective_db")
    return 0


if __name__ == "__main__":
    sys.exit(main())
