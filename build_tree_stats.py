#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_tree_stats.py — Combien d'espèces décrites dans chaque grand groupe.

L'écran « arbre du vivant » montrait les effectifs de notre collection, donc
la forme de ce qui a été photographié plutôt que celle du vivant. Ce script
compte les **espèces réellement décrites** dans la taxonomie du NCBI, tous
groupes confondus, images ou non.

    python3 build_tree_stats.py     # -> tree_stats.json

Ne comptent que les taxons de rang « species » : la taxonomie contient aussi
les genres, familles et sous-espèces, et les additionner gonflerait les
chiffres sans rien vouloir dire.
"""

import json
import os
import sys
import time
from collections import Counter

TAXDUMP = "new_taxdump"
SEP = "\t|\t"
OUT = "tree_stats.json"

# Suffixes de rang viraux imposés par l'ICTV (voir data_pipeline.py).
VIRAL_SUFFIXES = ("viria", "virae", "viricota", "viricetes", "virales",
                  "viridae", "viriformidae", "satellitidae", "viroidae")


def load_species_ranks(path):
    """Les taxid de rang « species », et eux seuls."""
    species = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            cols = line.split(SEP, 3)
            if len(cols) >= 3 and cols[2].strip() == "species":
                species.add(cols[0].strip())
    return species


def main():
    nodes = os.path.join(TAXDUMP, "nodes.dmp")
    ranked = os.path.join(TAXDUMP, "rankedlineage.dmp")
    for p in (nodes, ranked):
        if not os.path.exists(p):
            print(f"{p} introuvable.", file=sys.stderr)
            return 1

    t0 = time.time()
    print("Lecture des rangs…", flush=True)
    species = load_species_ranks(nodes)
    print(f"  {len(species):,} taxons de rang espèce", flush=True)

    by_super = Counter()
    by_kingdom = Counter()
    viruses = 0
    seen = 0

    print("Lecture des lignées…", flush=True)
    with open(ranked, encoding="utf-8") as fh:
        for line in fh:
            taxid = line.partition(SEP)[0]
            if taxid not in species:
                continue
            cols = [c.strip() for c in
                    line.rstrip("\n").rstrip("\t|").split(SEP)]
            cols += [""] * (10 - len(cols))
            name, kingdom, superkingdom = cols[1], cols[8], cols[9]
            seen += 1
            if seen % 500000 == 0:
                print(f"  {seen:,} espèces classées", flush=True)

            candidates = [name, superkingdom, kingdom, cols[7], cols[5], cols[4]]
            if superkingdom == "Viruses" or any(
                    c.endswith(VIRAL_SUFFIXES) for c in candidates if c):
                viruses += 1
                continue
            if superkingdom:
                by_super[superkingdom] += 1
            if superkingdom == "Eukaryota":
                by_kingdom[kingdom or "(autres eucaryotes)"] += 1

    euk = by_super.get("Eukaryota", 0)
    stats = {
        "source": "NCBI Taxonomy (new_taxdump), taxons de rang species",
        "total_species": seen,
        "bacteria": by_super.get("Bacteria", 0),
        "archaea": by_super.get("Archaea", 0),
        "eukaryota": euk,
        "metazoa": by_kingdom.get("Metazoa", 0),
        "viridiplantae": by_kingdom.get("Viridiplantae", 0),
        "fungi": by_kingdom.get("Fungi", 0),
        "autres_eucaryotes": euk - sum(by_kingdom.get(k, 0) for k in
                                       ("Metazoa", "Viridiplantae", "Fungi")),
        "viruses": viruses,
    }
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, ensure_ascii=False, indent=2)
        fh.write("\n")

    print("-" * 52)
    for k, v in stats.items():
        if isinstance(v, int):
            print(f"  {k:<22} {v:>10,}")
    print(f"\n{OUT} écrit en {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
