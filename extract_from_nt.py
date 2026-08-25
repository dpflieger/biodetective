#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_from_nt.py — À LANCER SUR LE SERVEUR QUI HÉBERGE nt.

Extrait de nt **toutes** les séquences des taxons de BioDetective, c'est-à-dire
ceux qui possèdent une image. Autrement dit : un petit nt restreint à notre
liste de taxid. Aucun plafond, aucune sélection.

    python3 extract_from_nt.py --db /chemin/vers/nt
    python3 extract_from_nt.py --db /chemin/vers/nt --survey   # mesurer d'abord

Ne nécessite que blastdbcmd et python3.

Le taxid est écrit dans l'identifiant FASTA (`>taxid|accession description`),
pour que la machine de démonstration n'ait besoin d'aucun fichier de taxonomie
pour relier un hit à sa fiche. La description est conservée : c'est elle qui
permet d'annoncer « chromosome 1 » plutôt qu'un numéro d'accession nu.

Note de dimensionnement : la E-value est proportionnelle à la taille de la
banque, mais la marge est confortable. Mesuré avec -dbsize, un brin de
24 briques donne E = 1e-6 contre 32 Mpb, 7,6e-5 contre 3,2 Gpb et 0,004 contre
320 Gpb. Il n'y a donc pas lieu de rogner.
"""

import argparse
import os
import subprocess
import sys


def human(n):
    for unit in ("pb", "kpb", "Mpb", "Gpb"):
        if n < 1000:
            return f"{n:,.1f} {unit}"
        n /= 1000
    return f"{n:,.1f} Tpb"


def survey(db, taxids):
    """Mesure le volume disponible sans rien extraire."""
    print("Inventaire (plusieurs minutes sur nt)…", flush=True)
    proc = subprocess.Popen(
        ["blastdbcmd", "-db", db, "-taxidlist", taxids, "-outfmt", "%T\t%l"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        errors="replace")
    n = total = 0
    species = set()
    longest = 0
    for line in proc.stdout:
        p = line.split("\t")
        if len(p) < 2:
            continue
        try:
            length = int(p[1])
        except ValueError:
            continue
        n += 1
        total += length
        longest = max(longest, length)
        species.add(p[0])
        if n % 1000000 == 0:
            print(f"  {n:,} séquences, {human(total)}", flush=True)
    proc.wait()
    if proc.returncode != 0:
        print(f"blastdbcmd a échoué : {proc.stderr.read()[:500]}",
              file=sys.stderr)
        return 1
    print("-" * 55)
    print(f"{n:,} séquences")
    print(f"{len(species):,} taxons")
    print(f"{human(total)} au total")
    print(f"séquence la plus longue : {human(longest)}")
    print(f"\nFASTA attendu : ~{total/1e9:.1f} Go sur disque")
    print(f"E-value attendue pour 24 pb : ~{1.06e-6 * (total/32e6):.1e}")
    return 0


def extract(db, taxids, out, min_length):
    """Écrit toutes les séquences des taxons demandés, sans plafond."""
    print(f"Extraction depuis {db}…", flush=True)
    proc = subprocess.Popen(
        ["blastdbcmd", "-db", db, "-taxidlist", taxids,
         "-outfmt", "%a\t%T\t%t\t%s"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        errors="replace")

    n = skipped = total = 0
    species = set()
    with open(out, "w") as fh:
        for line in proc.stdout:
            p = line.rstrip("\n").split("\t", 3)
            if len(p) < 4:
                continue
            acc, taxid, title, seq = p
            if len(seq) < min_length:
                skipped += 1
                continue
            # Le titre est conservé : c'est lui qui dira « chromosome 1 »
            # plutôt qu'un simple numéro d'accession, et sans lui la position
            # rendue par BLAST ne se rapporte à rien de nommable.
            title = title.replace("\t", " ").strip()
            fh.write(f">{taxid}|{acc} {title}\n{seq}\n")
            n += 1
            total += len(seq)
            species.add(taxid)
            if n % 100000 == 0:
                print(f"  {n:,} séquences, {human(total)}", flush=True)
    proc.wait()
    if proc.returncode != 0:
        print(f"blastdbcmd a échoué : {proc.stderr.read()[:500]}",
              file=sys.stderr)
        return 1

    size = os.path.getsize(out)
    print("-" * 55)
    print(f"{n:,} séquences, {len(species):,} taxons, {human(total)}")
    if skipped:
        print(f"{skipped:,} séquences ignorées (< {min_length} pb)")
    print(f"{out} : {size/1e9:.2f} Go")
    print(f"E-value attendue pour 24 pb : ~{1.06e-6 * (total/32e6):.1e}")
    print(f"\nÀ rapatrier, puis sur la machine de démonstration :")
    print(f"    python3 prepare_blastdb.py --fasta {os.path.basename(out)}")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="Extrait de nt toutes les séquences des taxons imagés")
    ap.add_argument("--db", required=True, help="chemin de la banque nt")
    ap.add_argument("--taxids", default="taxids.txt")
    ap.add_argument("--out", default="biodetective_subset.fasta")
    ap.add_argument("--min-length", type=int, default=0,
                    help="ignorer les séquences plus courtes (0 = tout garder)")
    ap.add_argument("--survey", action="store_true",
                    help="mesurer le volume sans rien extraire")
    args = ap.parse_args()

    if not os.path.exists(args.taxids):
        print(f"{args.taxids} introuvable.", file=sys.stderr)
        return 1
    n = sum(1 for l in open(args.taxids) if l.strip())
    print(f"{n:,} taxons demandés")

    if args.survey:
        return survey(args.db, args.taxids)
    return extract(args.db, args.taxids, args.out, args.min_length)


if __name__ == "__main__":
    sys.exit(main())
