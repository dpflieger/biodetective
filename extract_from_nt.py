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


# Plafonds proposés lors de l'inventaire, en bases par espèce.
SURVEY_CAPS = [50_000, 200_000, 1_000_000, 5_000_000, 20_000_000,
               100_000_000, 500_000_000, None]


def survey(db, taxids, priority):
    """Mesure le volume disponible et simule plusieurs plafonds.

    Ne lit que les métadonnées, jamais les séquences : bien plus rapide que
    l'extraction, et suffisant pour choisir un plafond en connaissance de
    cause plutôt qu'au jugé.
    """
    print("Inventaire des métadonnées (plusieurs minutes sur nt)…",
          flush=True)
    proc = subprocess.Popen(
        ["blastdbcmd", "-db", db, "-taxidlist", taxids, "-outfmt", "%T\t%l"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        errors="replace")

    per_taxid = {}
    n = total = longest = 0
    for line in proc.stdout:
        p = line.split("\t")
        if len(p) < 2:
            continue
        try:
            length = int(p[1])
        except ValueError:
            continue
        taxid = p[0]
        n += 1
        total += length
        longest = max(longest, length)
        per_taxid[taxid] = per_taxid.get(taxid, 0) + length
        if n % 2000000 == 0:
            print(f"  {n:,} séquences, {human(total)}", flush=True)
    proc.wait()
    if proc.returncode != 0:
        print(f"blastdbcmd a échoué : {proc.stderr.read()[:500]}",
              file=sys.stderr)
        return 1

    print("-" * 62)
    print(f"{n:,} séquences, {len(per_taxid):,} taxons, {human(total)}")
    print(f"séquence la plus longue : {human(longest)}")

    # Les quelques espèces qui pèsent le plus lourd : ce sont elles qui
    # décident du volume final.
    top = sorted(per_taxid.items(), key=lambda kv: -kv[1])[:8]
    print("\nEspèces les plus volumineuses :")
    for taxid, v in top:
        print(f"    taxid {taxid:<10} {human(v)}")

    print("\nVolume selon le plafond par espèce :")
    print(f"    {'plafond':>14}  {'total':>12}  {'E (24 pb)':>11}  "
          f"{'FASTA':>9}")
    for cap in SURVEY_CAPS:
        vol = sum(min(v, cap) if cap else v for v in per_taxid.values())
        if priority and cap:
            # les taxons prioritaires ne sont pas plafonnés
            vol += sum(max(0, per_taxid.get(t, 0) - cap) for t in priority)
        e = 1.06e-6 * (vol / 32e6)
        label = human(cap) if cap else "aucun"
        print(f"    {label:>14}  {human(vol):>12}  {e:>11.1e}  "
              f"{vol/1e9:>7.1f} Go")
    if priority:
        print(f"\n({len(priority)} taxons prioritaires jamais plafonnés)")
    print("\nChoisir ensuite : --budget <bases>")
    return 0


def extract(db, taxids, out, min_length, budget, priority):
    """Écrit les séquences des taxons demandés, sous plafond par espèce.

    Le plafond porte sur le total de bases par espèce, pas sur le nombre de
    séquences : sans lui, quelques organismes très séquencés — l'humain, le
    maïs, le blé — pèsent à eux seuls des centaines de gigaoctets.

    Les taxons prioritaires (ceux de la démonstration) n'ont pas de plafond :
    il faut que leur génome entier soit présent pour qu'une séquence prise
    n'importe où les retrouve.
    """
    print(f"Extraction depuis {db}…", flush=True)
    if budget:
        print(f"plafond : {human(budget)} par espèce, "
              f"{len(priority)} taxons exemptés", flush=True)
    proc = subprocess.Popen(
        ["blastdbcmd", "-db", db, "-taxidlist", taxids,
         "-outfmt", "%a\t%T\t%t\t%s"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        errors="replace")

    n = skipped = capped = total = 0
    used = {}
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
            if budget and taxid not in priority:
                spent = used.get(taxid, 0)
                if spent >= budget:
                    capped += 1
                    continue
                used[taxid] = spent + len(seq)
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
        print(f"{skipped:,} ignorées (< {min_length} pb)")
    if capped:
        print(f"{capped:,} écartées par le plafond")
    print(f"{out} : {size/1e9:.2f} Go")
    print(f"E-value attendue pour 24 pb : ~{1.06e-6 * (total/32e6):.1e}")
    print(f"\nÀ rapatrier, puis sur la machine de démonstration :")
    print(f"    python3 prepare_blastdb.py --fasta {os.path.basename(out)}")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="Extrait de nt les séquences des taxons imagés")
    ap.add_argument("--db", required=True, help="chemin de la banque nt")
    ap.add_argument("--taxids", default="taxids.txt")
    ap.add_argument("--priority", default="priority_taxids.txt",
                    help="taxons jamais plafonnés (ceux de la démonstration)")
    ap.add_argument("--out", default="biodetective_subset.fasta")
    ap.add_argument("--budget", type=int, default=5_000_000,
                    help="bases par espèce ; 0 = aucun plafond")
    ap.add_argument("--min-length", type=int, default=0,
                    help="ignorer les séquences plus courtes")
    ap.add_argument("--survey", action="store_true",
                    help="mesurer les volumes sans rien extraire")
    args = ap.parse_args()

    if not os.path.exists(args.taxids):
        print(f"{args.taxids} introuvable.", file=sys.stderr)
        return 1
    n = sum(1 for l in open(args.taxids) if l.strip())

    priority = set()
    if os.path.exists(args.priority):
        priority = {l.strip() for l in open(args.priority) if l.strip()}
    else:
        print(f"({args.priority} absent : aucun taxon exempté de plafond)")
    print(f"{n:,} taxons demandés, {len(priority)} prioritaires")

    if args.survey:
        return survey(args.db, args.taxids, priority)
    return extract(args.db, args.taxids, args.out, args.min_length,
                   args.budget, priority)


if __name__ == "__main__":
    sys.exit(main())
