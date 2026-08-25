#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_blastdb.py — Nettoie le FASTA issu de nt puis construit la banque.

    python3 prepare_blastdb.py --fasta taxids.fasta

Deux corrections indispensables avant makeblastdb :

1. `blastdbcmd -taxidlist` descend dans la hiérarchie : le fichier contient des
   sous-espèces et des variétés dont le taxid n'existe pas dans
   biodetective.db. Un hit sur l'une d'elles ne trouverait aucune fiche et
   serait écarté, alors que c'est souvent la bonne réponse. On les rattache
   donc à leur ancêtre présent dans la base, via nodes.dmp.

2. Ce qui ne se rattache à rien est retiré : une séquence qu'on ne saurait ni
   nommer ni illustrer n'a rien à faire dans la banque, et alourdirait les
   E-values pour rien.
"""

import argparse
import os
import sqlite3
import subprocess
import sys

NODES = os.path.join("new_taxdump", "nodes.dmp")
SEP = "\t|\t"


def load_parents(path):
    """taxid -> taxid parent, depuis nodes.dmp."""
    parents = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            cols = line.split(SEP, 3)
            if len(cols) >= 2:
                parents[cols[0].strip()] = cols[1].strip()
    return parents


def resolve(taxid, wanted, parents, cache):
    """Remonte jusqu'au premier ancêtre présent dans la base. None sinon."""
    if taxid in wanted:
        return taxid
    if taxid in cache:
        return cache[taxid]
    chain, cur = [], taxid
    seen = set()
    while cur and cur not in seen:
        seen.add(cur)
        if cur in wanted:
            break
        if cur in cache:
            cur = cache[cur]
            break
        chain.append(cur)
        nxt = parents.get(cur)
        if nxt is None or nxt == cur:      # racine atteinte
            cur = None
            break
        cur = nxt
    result = cur if (cur and cur in wanted) else None
    for t in chain:
        cache[t] = result
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", default="taxids.fasta")
    ap.add_argument("--out", default="biodetective_subset.fasta")
    ap.add_argument("--db", default="blastdb/biodetective")
    ap.add_argument("--skip-makeblastdb", action="store_true",
                    help="écrire le FASTA nettoyé au lieu de construire "
                         "la banque (occupe le disque en double)")
    args = ap.parse_args()

    conn = sqlite3.connect("biodetective.db")
    wanted = {str(r[0]) for r in
              conn.execute("SELECT taxonomy_id FROM organisms")}
    conn.close()
    print(f"{len(wanted):,} espèces avec fiche")

    if not os.path.exists(NODES):
        print(f"{NODES} introuvable : impossible de rattacher les "
              f"sous-espèces.", file=sys.stderr)
        return 1
    parents = load_parents(NODES)
    print(f"{len(parents):,} relations de parenté chargées")

    cache = {}
    kept = remapped = dropped = 0
    total_bp = 0
    species = set()

    # Le FASTA d'entrée peut peser plusieurs centaines de gigaoctets. On
    # n'en écrit pas une seconde copie : les en-têtes sont réécrits à la
    # volée et poussés directement dans makeblastdb, qui lit stdin. Sans
    # cela il faudrait deux fois la place sur le disque.
    mkdb = None
    sink = None
    if args.skip_makeblastdb:
        sink = open(args.out, "w")
    else:
        os.makedirs(os.path.dirname(args.db) or ".", exist_ok=True)
        sys.path.insert(0, ".")
        import blast_search
        exe = os.path.join(os.path.dirname(blast_search.BLASTN or ""),
                           "makeblastdb") if blast_search.BLASTN \
            else "makeblastdb"
        print(f"Construction de la banque en flux avec {exe}…", flush=True)
        mkdb = subprocess.Popen(
            [exe, "-in", "-", "-dbtype", "nucl", "-out", args.db,
             "-title", "BioDetective"],
            stdin=subprocess.PIPE, text=True)
        sink = mkdb.stdin

    try:
        with open(args.fasta) as src:
            emit = False
            for line in src:
                if line.startswith(">"):
                    raw = line[1:].rstrip("\n")
                    taxid, _, rest = raw.partition("|")
                    acc, _, title = rest.partition(" ")
                    target = resolve(taxid, wanted, parents, cache)
                    if target is None:
                        emit = False
                        dropped += 1
                        continue
                    emit = True
                    kept += 1
                    if target != taxid:
                        remapped += 1
                    species.add(target)
                    sink.write(f">{target}|{acc}"
                               + (f" {title}" if title else "") + "\n")
                    if kept % 500000 == 0:
                        print(f"  {kept:,} séquences, {total_bp/1e9:.1f} Gpb",
                              flush=True)
                elif emit:
                    sink.write(line)
                    total_bp += len(line.strip())
    finally:
        if mkdb is not None:
            sink.close()
            mkdb.wait()
        else:
            sink.close()

    print(f"\n{kept:,} séquences conservées")
    print(f"  dont {remapped:,} rattachées à l'espèce parente")
    print(f"{dropped:,} écartées (aucun ancêtre en base)")
    print(f"{len(species):,} espèces représentées")
    print(f"{total_bp/1e6:,.0f} Mpb  ->  E-value attendue pour 24 pb : "
          f"~{1.06e-6*(total_bp/32e6):.1e}")
    if mkdb is not None:
        return mkdb.returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())
