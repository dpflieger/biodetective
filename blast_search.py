#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
blast_search.py — Identification par BLAST contre la banque BioDetective.

Remplace la table de correspondance exacte : on interroge une vraie banque
construite à partir de nt, restreinte aux espèces qui possèdent une image.
On récupère donc un classement par E-value, un pourcentage d'identité réel,
et un alignement à montrer.

    python3 blast_search.py ATCGATCGATCGATCGATCGATCG

La banque doit avoir été construite au préalable :
    makeblastdb -in biodetective_subset.fasta -dbtype nucl -out blastdb/biodetective
"""

import asyncio
import os
import shutil
import subprocess
import sys

BLAST_DB = os.environ.get("BIODETECTIVE_BLASTDB", "blastdb/biodetective")


def find_blastn():
    """Localise blastn sans exiger que conda soit activé.

    Sur la machine de démonstration, python3 est celui du système
    (/usr/bin/python3) alors que blastn vient de conda : il n'est donc pas
    dans le PATH. Chercher au bon endroit évite un « blastn introuvable »
    au pire moment.
    """
    forced = os.environ.get("BIODETECTIVE_BLASTN")
    if forced:
        return forced
    found = shutil.which("blastn")
    if found:
        return found
    candidates = [
        os.path.join(os.path.dirname(sys.executable), "blastn"),
        os.path.expanduser("~/miniforge3/bin/blastn"),
        os.path.expanduser("~/miniconda3/bin/blastn"),
        os.path.expanduser("~/anaconda3/bin/blastn"),
        os.path.expanduser("~/mambaforge/bin/blastn"),
    ]
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return None


BLASTN = find_blastn()

# blastn-short est réglé pour les requêtes courtes : c'est le cas, le
# Brickopore produit 24 bases.
TASK = "blastn-short"

# Seuil volontairement large. Contre une banque de ~25 Mpb, une requête de
# 24 pb parfaite tombe vers 1e-6 ; on laisse passer bien au-delà pour montrer
# aussi les correspondances partielles, qui font partie de la démonstration.
EVALUE = 1.0
MAX_HITS = 50

# Au-delà, on considère que l'organisme n'est pas identifié. Une brique fausse
# sur 24 donne ~96 % : il faut rester tolérant.
MIN_IDENTITY = 85.0
MIN_COVERAGE = 80.0

_FIELDS = ("sseqid pident length mismatch gapopen qstart qend "
           "sstart send evalue bitscore qseq sseq qlen")


def blast_available():
    return BLASTN is not None


def db_available():
    return os.path.exists(BLAST_DB + ".nin") or os.path.exists(BLAST_DB + ".nal")


def _parse(stdout):
    """Transforme la sortie tabulée en liste de hits ordonnés."""
    hits = []
    for line in stdout.splitlines():
        p = line.rstrip("\n").split("\t")
        if len(p) < 14:
            continue
        sseqid = p[0]
        # L'identifiant est « taxid|accession » : le taxid vient de la banque,
        # pas d'un fichier de taxonomie annexe.
        taxid, _, accession = sseqid.partition("|")
        try:
            taxid = int(taxid)
        except ValueError:
            continue
        qlen = int(p[13]) or 1
        length = int(p[2])
        hits.append({
            "taxonomy_id": taxid,
            "accession": accession or sseqid,
            "percent_identity": round(float(p[1]), 1),
            "alignment_length": length,
            "mismatches": int(p[3]),
            "gaps": int(p[4]),
            "query_start": int(p[5]), "query_end": int(p[6]),
            "subject_start": int(p[7]), "subject_end": int(p[8]),
            "evalue": float(p[9]),
            "bitscore": float(p[10]),
            "query_seq": p[11],
            "subject_seq": p[12],
            "coverage": round(100.0 * length / qlen, 1),
        })
    # BLAST trie déjà par bitscore, mais on ne s'en remet pas à l'ordre reçu.
    hits.sort(key=lambda h: (h["evalue"], -h["bitscore"]))
    return hits


def midline(qseq, sseq):
    """Ligne centrale d'un alignement BLAST : « | » là où ça correspond."""
    return "".join("|" if a == b and a != "-" else " "
                   for a, b in zip(qseq, sseq))


def is_confident(hit):
    """Un hit assez bon pour annoncer un organisme à un enfant."""
    return (hit["percent_identity"] >= MIN_IDENTITY
            and hit["coverage"] >= MIN_COVERAGE)


def _cmd(query_path):
    return [
        BLASTN, "-task", TASK, "-query", query_path, "-db", BLAST_DB,
        "-evalue", str(EVALUE), "-max_target_seqs", str(MAX_HITS),
        "-num_threads", "4", "-dust", "no", "-soft_masking", "false",
        "-outfmt", "6 " + _FIELDS,
    ]


async def search(sequence, workdir="/tmp"):
    """Lance blastn en tâche de fond et retourne la liste des hits."""
    if not blast_available():
        raise RuntimeError(
            "blastn est introuvable. Installer BLAST+ "
            "(conda install -c bioconda blast) ou indiquer le chemin dans "
            "BIODETECTIVE_BLASTN.")
    if not db_available():
        raise RuntimeError(
            f"Banque BLAST absente ({BLAST_DB}). "
            f"Lancer makeblastdb, voir CLAUDE.md.")

    path = os.path.join(workdir, f"bd_query_{os.getpid()}_{id(sequence)}.fa")
    with open(path, "w") as fh:
        fh.write(f">requete\n{sequence}\n")
    try:
        proc = await asyncio.create_subprocess_exec(
            *_cmd(path), stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE)
        out, err = await proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(
                f"blastn a échoué : {err.decode('utf-8', 'replace')[:300]}")
        return _parse(out.decode("utf-8", "replace"))
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


def search_sync(sequence, workdir="/tmp"):
    """Version bloquante, pour les scripts et les tests."""
    path = os.path.join(workdir, f"bd_query_{os.getpid()}.fa")
    with open(path, "w") as fh:
        fh.write(f">requete\n{sequence}\n")
    try:
        r = subprocess.run(_cmd(path), capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"blastn a échoué : {r.stderr[:300]}")
        return _parse(r.stdout)
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


def main():
    if len(sys.argv) != 2:
        print(__doc__.strip())
        return 2
    from identify import normalize, SequenceError
    try:
        seq = normalize(sys.argv[1])
    except SequenceError as exc:
        print(f"Séquence invalide : {exc}")
        return 1

    hits = search_sync(seq)
    print(f"Requête : {seq} ({len(seq)} pb)")
    print(f"{len(hits)} hits\n")
    for h in hits[:5]:
        flag = "✓" if is_confident(h) else " "
        print(f"{flag} taxid {h['taxonomy_id']:<9} {h['accession']:<14} "
              f"id {h['percent_identity']:5.1f} %  cov {h['coverage']:5.1f} %  "
              f"E {h['evalue']:.2g}  score {h['bitscore']}")
    if hits:
        h = hits[0]
        print(f"\nAlignement du meilleur hit :")
        print(f"  requete  {h['query_start']:>4} {h['query_seq']} {h['query_end']}")
        print(f"           {'':>4} {midline(h['query_seq'], h['subject_seq'])}")
        print(f"  sujet    {h['subject_start']:>4} {h['subject_seq']} {h['subject_end']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
