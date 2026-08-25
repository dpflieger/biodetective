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
from datetime import datetime
import re
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

# Chaque analyse laisse une trace vérifiable : la commande exacte, le tableau
# que l'application a réellement lu, et le rapport détaillé de blastn. Mettre
# BIODETECTIVE_BLAST_RESULTS à vide désactive l'archivage.
RESULTS_DIR = os.environ.get("BIODETECTIVE_BLAST_RESULTS", "blast_results")

# staxids vient de la banque elle-même (nt + taxdb). Sur notre banque locale
# construite à partir d'un FASTA « taxid|accession », il est vide et l'on
# retombe sur l'identifiant. Les deux montages fonctionnent donc.
_FIELDS = ("sseqid pident length mismatch gapopen qstart qend "
           "sstart send evalue bitscore qseq sseq qlen slen staxids stitle")


def blast_available():
    return BLASTN is not None


def db_available():
    return os.path.exists(BLAST_DB + ".nin") or os.path.exists(BLAST_DB + ".nal")


# Repère le type de molécule dans la description NCBI, pour dire à un enfant
# « chromosome 1 » ou « ADN du chloroplaste » plutôt qu'un numéro d'accession.
LOCUS_PATTERNS = [
    (re.compile(r"\bchromosome\s+(\w+)", re.I), "chromosome {}"),
    (re.compile(r"\blinkage group\s+(\w+)", re.I), "groupe de liaison {}"),
    (re.compile(r"\bmitochondri", re.I), "génome mitochondrial"),
    (re.compile(r"\bchloroplast|\bplastid", re.I), "génome chloroplastique"),
    (re.compile(r"\bplasmid\s+(\S+)", re.I), "plasmide {}"),
    (re.compile(r"\bscaffold\s+(\S+)", re.I), "scaffold {}"),
    (re.compile(r"\bcontig\s+(\S+)", re.I), "contig {}"),
    (re.compile(r"\bcomplete genome", re.I), "génome complet"),
]


def describe_locus(title):
    """Nomme la molécule touchée. None si la description ne dit rien."""
    if not title:
        return None
    for rx, label in LOCUS_PATTERNS:
        m = rx.search(title)
        if m:
            return label.format(*m.groups()) if m.groups() else label
    return None


def _write_archive(archive_id, sequence, cmd, table, report):
    """Dépose le résultat brut dans RESULTS_DIR. Ne lève jamais.

    L'archivage ne doit sous aucun prétexte faire échouer une analyse en
    cours devant un enfant : toute erreur d'écriture est avalée.
    """
    if not RESULTS_DIR:
        return None
    try:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        path = os.path.join(RESULTS_DIR, f"{stamp}_{archive_id[:8]}.txt")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(f"# BioDetective — analyse du "
                     f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            fh.write(f"# job      : {archive_id}\n")
            fh.write(f"# requete  : {sequence} ({len(sequence)} pb)\n")
            fh.write(f"# banque   : {BLAST_DB}\n")
            fh.write(f"# commande : {' '.join(cmd)}\n")
            fh.write("\n===== Tableau lu par l'application (outfmt 6) =====\n")
            fh.write("# " + "\t".join(_FIELDS.split()) + "\n")
            fh.write(table or "(aucun hit)\n")
            if report:
                fh.write("\n===== Rapport blastn detaille (outfmt 0) =====\n")
                fh.write(report)
        return path
    except OSError:
        return None


INDEX_COLUMNS = ("date", "job", "sequence", "resultat", "taxid",
                 "identite", "couverture", "evalue", "score", "fichier")


def append_index(row):
    """Ajoute une ligne au récapitulatif TSV, pour parcourir la journée
    sans ouvrir les fichiers un par un."""
    if not RESULTS_DIR:
        return
    try:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        path = os.path.join(RESULTS_DIR, "index.tsv")
        new = not os.path.exists(path)
        with open(path, "a", encoding="utf-8") as fh:
            if new:
                fh.write("\t".join(INDEX_COLUMNS) + "\n")
            fh.write("\t".join(str(row.get(c, "")) for c in INDEX_COLUMNS)
                     + "\n")
    except OSError:
        pass


def append_verdict(path, text):
    """Ajoute la conclusion de l'application au fichier archivé."""
    if not path:
        return
    try:
        with open(path, "a", encoding="utf-8") as fh:
            fh.write("\n===== Verdict de l'application =====\n")
            fh.write(text.rstrip() + "\n")
    except OSError:
        pass


# Étiquettes de banque du format historique du NCBI : « ref|NC_003070.9| ».
DB_TAGS = {"ref", "gb", "emb", "dbj", "tpg", "tpe", "tpd", "gnl", "lcl",
           "pir", "prf", "sp", "pdb", "pat", "bbs", "gi"}


def split_seqid(sseqid):
    """(taxid éventuel, accession) à partir d'un identifiant BLAST.

    Trois formes coexistent selon la banque :
        3702|NC_003070.9   notre FASTA maison, le taxid est en tête
        ref|NC_003070.9|   nt construit avec -parse_seqids
        NC_003070.9        identifiant nu
    """
    parts = [x for x in sseqid.split("|") if x]
    if not parts:
        return None, sseqid
    if parts[0].isdigit() and len(parts) > 1:
        return int(parts[0]), parts[1]          # taxid|accession
    if parts[0].lower() in DB_TAGS and len(parts) > 1:
        # gi|123|ref|NC_… : la dernière étiquette connue précède l'accession
        for i in range(len(parts) - 1, 0, -1):
            if parts[i - 1].lower() in DB_TAGS:
                return None, parts[i]
        return None, parts[1]
    return None, parts[0]


def _parse(stdout):
    """Transforme la sortie tabulée en liste de hits ordonnés."""
    hits = []
    for line in stdout.splitlines():
        p = line.rstrip("\n").split("\t")
        if len(p) < 14:
            continue
        sseqid = p[0]
        # Deux origines possibles pour le taxid : la colonne staxids quand la
        # banque en porte (nt avec taxdb), sinon l'identifiant « taxid|acc »
        # de notre FASTA maison.
        staxids = p[15] if len(p) > 15 else ""
        embedded, accession = split_seqid(sseqid)

        taxid = None
        # Sans taxdb, blastn écrit « 0 » et non une colonne vide : le prendre
        # pour argent comptant donnerait taxid 0 sur toute la ligne.
        first = staxids.split(";")[0].strip()
        if first and first not in ("0", "N/A"):
            try:
                taxid = int(first)
            except ValueError:
                taxid = None
        if taxid is None:
            taxid = embedded
        if taxid is None:
            continue
        qlen = int(p[13]) or 1
        length = int(p[2])
        slen = int(p[14]) if len(p) > 14 and p[14].isdigit() else None
        # stitle peut contenir des tabulations : on reprend tout le reste.
        stitle = "\t".join(p[16:]).strip() if len(p) > 16 else ""
        # blastdbcmd recopie parfois l'identifiant en tête du titre.
        if stitle.startswith(sseqid):
            stitle = stitle[len(sseqid):].strip()
        sstart, send = int(p[7]), int(p[8])
        hits.append({
            "taxonomy_id": taxid,
            "accession": accession or sseqid,
            "percent_identity": round(float(p[1]), 1),
            "alignment_length": length,
            "mismatches": int(p[3]),
            "gaps": int(p[4]),
            "query_start": int(p[5]), "query_end": int(p[6]),
            "subject_start": sstart,
            "subject_end": send,
            "subject_length": slen,
            "subject_title": stitle,
            # sstart > send : l'alignement est sur le brin complémentaire.
            "strand": "moins" if sstart > send else "plus",
            "locus": describe_locus(stitle),
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


def _cmd(query_path, outfmt="6 " + _FIELDS):
    return [
        BLASTN, "-task", TASK, "-query", query_path, "-db", BLAST_DB,
        "-evalue", str(EVALUE), "-max_target_seqs", str(MAX_HITS),
        "-num_threads", "4", "-dust", "no", "-soft_masking", "false",
        "-outfmt", outfmt,
    ]


async def search(sequence, workdir="/tmp", archive_id=None):
    """Lance blastn en tâche de fond et retourne (hits, chemin d'archive)."""
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
    async def run(cmd):
        proc = await asyncio.create_subprocess_exec(
            *cmd, stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE)
        out, err = await proc.communicate()
        return proc.returncode, out.decode("utf-8", "replace"), err.decode(
            "utf-8", "replace")

    try:
        cmd = _cmd(path)
        archiving = bool(archive_id and RESULTS_DIR)

        # Le rapport lisible demande un second appel : blastn n'émet qu'un
        # format à la fois. Les deux sont indépendants, donc lancés
        # ensemble — en série ils coûtaient 584 ms au lieu de 315.
        if archiving:
            (rc, table, err), (rc2, report, _e2) = await asyncio.gather(
                run(cmd), run(_cmd(path, outfmt="0")))
        else:
            rc, table, err = await run(cmd)
            rc2, report = 1, None

        if rc != 0:
            raise RuntimeError(f"blastn a échoué : {err[:300]}")

        archive = _write_archive(archive_id, sequence, cmd, table,
                                 report if rc2 == 0 else None) \
            if archiving else None
        return _parse(table), archive
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
        return _parse(r.stdout)   # sans archivage : usage scripts et tests
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
