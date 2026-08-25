#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_strips.py — Fabrique les brins préparés du Brickopore.

Depuis que l'identification passe par un vrai BLAST, une séquence inventée
ne sert plus à rien : les brins de la démonstration doivent être de vrais
fragments d'ADN, découpés dans les séquences de référence de la banque.

Chaque candidat est vérifié PAR BLAST : on ne garde un fragment que s'il
place le bon organisme en tête. Beaucoup de fragments tombent dans des
régions conservées et désignent une famille entière plutôt qu'une espèce.

    python3 make_strips.py --length 24
"""

import argparse
import json
import os
import random
import sqlite3
import sys

import blast_search as bs

FASTA_CANDIDATES = ["biodetective_subset.fasta", "dev_subset.fasta"]
OUT = "sequences.json"
LEGO_MAX_RUN = 2          # jamais 3 briques identiques d'affilée
TRIES_PER_SPECIES = 250


def load_fasta(path):
    seqs, taxid = {}, None
    for line in open(path):
        if line.startswith(">"):
            taxid = line[1:].strip().split("|")[0]
        elif taxid:
            seqs.setdefault(taxid, "")
            seqs[taxid] += line.strip().upper()
    return seqs


def longest_run(s):
    best = run = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def usable(frag):
    """Écarte ce qu'un enfant ne saurait pas assembler ni compter."""
    if set(frag) - set("ATCG"):
        return False
    if longest_run(frag) > LEGO_MAX_RUN:
        return False
    gc = (frag.count("G") + frag.count("C")) / len(frag)
    return 0.30 <= gc <= 0.70


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--length", type=int, default=24,
                    help="nombre de briques par brin")
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--tries", type=int, default=TRIES_PER_SPECIES,
                    help="positions testées au maximum par espèce")
    args = ap.parse_args()

    fasta = args.fasta or next((f for f in FASTA_CANDIDATES
                                if os.path.exists(f)), None)
    if not fasta:
        print("Aucun FASTA de référence. Lancer extract_from_nt.py d'abord.",
              file=sys.stderr)
        return 1
    if not bs.db_available():
        print(f"Banque BLAST absente ({bs.BLAST_DB}).", file=sys.stderr)
        return 1
    print(f"Références : {fasta}")
    print(f"Banque     : {bs.BLAST_DB}")

    seqs = load_fasta(fasta)
    names = {int(k): v for k, v in
             json.load(open("common_names_fr.json", encoding="utf-8")).items()}

    conn = sqlite3.connect("biodetective.db")
    conn.row_factory = sqlite3.Row

    rng = random.Random(20261003)
    entries, skipped = [], []

    for taxid, nom_fr in sorted(names.items(), key=lambda kv: kv[1]):
        ref = seqs.get(str(taxid))
        row = conn.execute(
            "SELECT scientific_name, organism_type, image_path "
            "FROM organisms WHERE taxonomy_id=?", (taxid,)).fetchone()
        if not ref or len(ref) < args.length + 20 or row is None:
            skipped.append((nom_fr, "pas de séquence de référence"))
            continue
        if not row["image_path"]:
            skipped.append((nom_fr, "pas d'image"))
            continue

        # Balayage systématique de toutes les positions utilisables, dans un
        # ordre mélangé. Le tirage aléatoire avec remise repassait sur les
        # mêmes positions et échouait sur les espèces à proches parents :
        # chez la cigogne ou la vigne, un 24-mer discriminant est rare.
        positions = [i for i in range(len(ref) - args.length + 1)
                     if usable(ref[i:i + args.length])]
        rng.shuffle(positions)

        chosen = None
        for start in positions[:args.tries]:
            frag = ref[start:start + args.length]
            # Le juge de paix : ce fragment désigne-t-il bien cette espèce ?
            hits = bs.search_sync(frag)
            if not hits:
                continue
            top = hits[0]
            if top["taxonomy_id"] == taxid and bs.is_confident(top):
                chosen = (frag, top, start)
                break

        if chosen is None:
            skipped.append((nom_fr, f"aucun fragment discriminant "
                                    f"({len(positions)} positions testées)"))
            continue

        frag, top, start = chosen
        entries.append({
            "sequence": frag,
            "taxonomy_id": taxid,
            "nom_fr": nom_fr,
            "nom_scientifique": row["scientific_name"],
            "type": row["organism_type"],
            "evalue_attendue": float(f"{top['evalue']:.2g}"),
            "position_reference": start + 1,
        })
        print(f"  {nom_fr:<26} {frag}  E={top['evalue']:.2g}")

    doc = {
        "_comment": [
            "Brins préparés du Brickopore : de VRAIS fragments d'ADN,",
            "découpés dans les séquences de référence de la banque BLAST.",
            f"{args.length} briques, jamais 3 identiques d'affilée.",
            "Chaque brin a été vérifié par BLAST : il place bien son",
            "organisme en tête. Régénérer avec make_strips.py après tout",
            "changement de banque.",
        ],
        "length": args.length,
        "sequences": entries,
    }
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, ensure_ascii=False, indent=2)
        fh.write("\n")

    print("-" * 60)
    print(f"{len(entries)} brins écrits dans {args.out}")
    if skipped:
        print(f"{len(skipped)} organismes écartés :")
        for n, why in skipped:
            print(f"    {n:<26} {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
