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

# Longueur minimale utile. Il n'y a PLUS de plafond par défaut : les
# chromosomes et génomes d'organites sont les bienvenus, c'est le budget par
# espèce qui borne le volume.
MIN_LEN = 200

# Budget par espèce, en bases. On prend les séquences par ordre de préférence
# jusqu'à l'épuiser. C'est ce qui remplace l'ancien « une seule séquence » :
# ce plafond-là rendait invisible tout ce qui n'était pas le marqueur, et une
# séquence prise ailleurs dans le génome ne donnait aucun hit.
DEFAULT_BUDGET = 20_000
PRIORITY_BUDGET = 300_000_000   # pour les organismes de la démonstration

# Marqueurs classiques de code-barres ADN, par ordre de préférence.
PREFERRED = [
    r"cytochrome c oxidase subunit (?:1|I)\b", r"\bCOI\b", r"\bCOX1\b",
    r"ribulose.*carboxylase", r"\brbcL\b", r"\bmatK\b",
    r"internal transcribed spacer", r"\bITS\b",
    r"16S ribosomal RNA", r"18S ribosomal RNA", r"28S ribosomal RNA",
    r"cytochrome b\b", r"\bcytb\b",
]
PREFERRED_RE = [re.compile(p, re.I) for p in PREFERRED]

# On n'écarte plus le génomique : c'est précisément ce qui manquait. Seuls
# restent bannis les enregistrements sans intérêt pour une identification.
REJECT_RE = re.compile(r"unverified|clone library", re.I)


# Après les marqueurs, on privilégie les organites : mitochondrie et
# chloroplaste sont petits, très séquencés, et couvrent beaucoup de terrain.
ORGANELLE_RE = re.compile(r"mitochondri|chloroplast|plastid", re.I)


def rank(title, length):
    """Note une séquence : plus c'est bas, mieux c'est."""
    if REJECT_RE.search(title):
        return (99, 0)
    for i, rx in enumerate(PREFERRED_RE):
        if rx.search(title):
            return (i, abs(length - 800))
    n = len(PREFERRED_RE)
    if ORGANELLE_RE.search(title):
        return (n, -length)          # organite : le plus complet d'abord
    return (n + 1, -length)          # le reste : le plus long d'abord


def main():
    ap = argparse.ArgumentParser(
        description="Extrait de nt les séquences des espèces imagées")
    ap.add_argument("--db", required=True, help="chemin de la banque nt")
    ap.add_argument("--taxids", default="taxids.txt")
    ap.add_argument("--priority", default="priority_taxids.txt",
                    help="taxons de la démonstration, servis largement")
    ap.add_argument("--out", default="biodetective_subset.fasta")
    ap.add_argument("--budget", type=int, default=DEFAULT_BUDGET,
                    help="bases par espèce ordinaire")
    ap.add_argument("--priority-budget", type=int, default=PRIORITY_BUDGET,
                    help="bases par espèce de la démonstration")
    ap.add_argument("--survey", action="store_true",
                    help="ne rien extraire, seulement mesurer le volume")
    args = ap.parse_args()

    wanted = {l.strip() for l in open(args.taxids) if l.strip()}
    priority = set()
    if os.path.exists(args.priority):
        priority = {l.strip() for l in open(args.priority) if l.strip()}
    print(f"{len(wanted)} taxons demandés, dont {len(priority)} prioritaires",
          flush=True)
    print(f"budget : {args.budget:,} b par espèce, "
          f"{args.priority_budget:,} b pour les prioritaires", flush=True)

    print("Inventaire des séquences disponibles (plusieurs minutes)…",
          flush=True)
    cmd = ["blastdbcmd", "-db", args.db, "-taxidlist", args.taxids,
           "-outfmt", "%a\t%T\t%l\t%t"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True,
                            errors="replace")

    cands = defaultdict(list)
    seen = available = 0
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t", 3)
        if len(parts) < 4:
            continue
        acc, taxid, length, title = parts
        seen += 1
        if seen % 1000000 == 0:
            print(f"  {seen:,} séquences examinées, "
                  f"{len(cands):,} espèces vues", flush=True)
        try:
            length = int(length)
        except ValueError:
            continue
        if length < MIN_LEN or taxid not in wanted:
            continue
        r = rank(title, length)
        if r[0] >= 99:
            continue
        available += length
        cands[taxid].append((r, acc, length))

    proc.wait()
    if proc.returncode != 0:
        print(f"blastdbcmd a échoué : {proc.stderr.read()[:500]}",
              file=sys.stderr)
        return 1
    print(f"{seen:,} séquences examinées ; {available/1e9:.1f} Gpb "
          f"disponibles pour {len(cands):,} espèces", flush=True)

    # Sélection sous budget, séquences les mieux notées d'abord.
    chosen, total = [], 0
    for taxid, lst in cands.items():
        lst.sort()
        budget = args.priority_budget if taxid in priority else args.budget
        used = 0
        for _r, acc, length in lst:
            if used and used + length > budget:
                continue
            chosen.append(acc)
            used += length
            if used >= budget:
                break
        total += used

    print(f"\n{len(chosen):,} séquences retenues, {total/1e6:,.0f} Mpb")
    print(f"E-value attendue pour 24 pb : ~{1.06e-6 * (total/32e6):.1e}")
    if args.survey:
        print("\n(--survey : rien n'a été extrait)")
        return 0

    with open("chosen_accessions.txt", "w") as fh:
        fh.write("\n".join(chosen) + "\n")

    print("\nExtraction du FASTA…", flush=True)
    n = written = 0
    with open(args.out, "w") as out:
        p2 = subprocess.Popen(
            ["blastdbcmd", "-db", args.db, "-entry_batch",
             "chosen_accessions.txt", "-outfmt", "%a\t%T\t%s"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            errors="replace")
        for line in p2.stdout:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            acc, taxid, seq = p[0], p[1], p[2]
            out.write(f">{taxid}|{acc}\n{seq}\n")
            n += 1
            written += len(seq)
            if n % 20000 == 0:
                print(f"  {n:,} séquences écrites", flush=True)
        p2.wait()

    print("-" * 55)
    print(f"{n:,} séquences, {written/1e6:,.0f} Mpb -> {args.out}")
    print(f"\nÀ rapatrier, puis : python3 prepare_blastdb.py "
          f"--fasta {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
