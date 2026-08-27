#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simulate.py — Analyse simulée : un organisme au hasard, une vraie statistique.

Pourquoi
--------
La banque BLAST utilisable est encore en construction (voir CLAUDE.md). En
attendant, la démonstration doit tourner de bout en bout : un enfant assemble
ses briques, colle sa séquence, et voit toujours un organisme. Ce module
remplace donc `blast_search.search()` en gardant exactement la même signature
et la même forme de résultat — le reste de l'application ne sait pas la
différence.

Ce qui est inventé, ce qui ne l'est pas
---------------------------------------
INVENTÉ  : le choix de l'organisme, le nombre de mésappariements, la position
           du hit sur le chromosome.
RÉEL     : l'organisme, sa photo, sa lignée, son génome, ses chromosomes et
           leurs numéros d'accession — tout vient de `biodetective.db` et de
           `genome_stats.json`. Les cousins sont ses vrais cousins.
           Le score et la E-value sont calculés par les **vraies** formules de
           Karlin-Altschul, avec les paramètres que blastn imprime lui-même
           pour `blastn-short` : les chiffres affichés sont ceux que BLAST
           donnerait pour cet alignement contre cette banque.

Les brins préparés de `sequences.json` continuent de rendre LEUR organisme :
ce sont les exemples à donner aux enfants, ils ne doivent pas tomber au hasard.

    python3 simulate.py GCAGAATAAGTGCATTGAACTTAA
"""

import json
import math
import os
import random
import sqlite3
import sys
from datetime import datetime

import blast_search

# --------------------------------------------------------------------------
# Statistique de Karlin-Altschul
# --------------------------------------------------------------------------
# Relevés dans nos propres rapports blastn (outfmt 0, tâche blastn-short,
# reward 1 / penalty -3). Ce ne sont pas des valeurs inventées : elles sont
# reproductibles avec `blastn -task blastn-short`.
LAMBDA, K, H = 1.37, 0.711, 1.31
REWARD, PENALTY = 1, -3

# Taille de la banque annoncée. Celle du sous-ensemble matérialisé qui est sur
# le disque — la E-value est proportionnelle à cette taille, donc l'annoncer
# juste est ce qui rend le chiffre crédible.
DB_SEQUENCES = int(os.environ.get("BIODETECTIVE_SIM_DBSEQ", 44772))
DB_LETTERS = int(os.environ.get("BIODETECTIVE_SIM_DBLEN", 32020836))
DB_NAME = "BioDetective"


def length_adjustment(qlen, db_letters=DB_LETTERS, db_seqs=DB_SEQUENCES):
    """Correction de bord de BLAST, par point fixe — comme le fait BLAST.

    Un alignement ne peut pas commencer trop près de la fin d'une séquence :
    BLAST retranche donc une longueur `ell` à la requête et à chaque séquence
    de la banque avant de calculer l'espace de recherche.

    Vérifié : pour une requête de 24 pb contre 44 772 séquences /
    32 020 836 lettres, on retrouve exactement les 313 940 280 que blastn
    imprime sous « Effective search space used ».
    """
    ell = 0.0
    for _ in range(20):
        m = max(1.0, qlen - ell)
        n = max(1.0, db_letters - db_seqs * ell)
        new = math.log(K * m * n) / H
        if abs(new - ell) < 1e-6:
            break
        ell = new
    ell = math.floor(ell)
    m = max(1.0, qlen - ell)
    n = max(1.0, db_letters - db_seqs * ell)
    return ell, m, n


def score_and_evalue(qlen, matches, mismatches):
    """(bitscore, E-value) pour un alignement sans trou, formules exactes."""
    raw = matches * REWARD + mismatches * PENALTY
    bits = (LAMBDA * raw - math.log(K)) / math.log(2)
    _ell, m, n = length_adjustment(qlen)
    evalue = K * m * n * math.exp(-LAMBDA * raw)
    return round(bits, 1), evalue


# --------------------------------------------------------------------------
# Fabrication de l'alignement
# --------------------------------------------------------------------------

COMPLEMENT = str.maketrans("ATCG", "TAGC")

# Nombre de mésappariements du meilleur hit. Sur 24 briques : 100 %, 95,8 %,
# 91,7 %, 87,5 %. Le 100 % reste le cas le plus fréquent — c'est celui qui
# fait plaisir — mais pas le seul, sinon le chiffre ne veut plus rien dire.
# La valeur tirée est ensuite bornée par k_max() : au-delà, la E-value
# dépasserait le seuil au-dessus duquel blastn ne rapporte plus rien, et
# afficher « E = 16 » ne tromperait personne.
TOP_MISMATCHES = [(0, 55), (1, 35), (2, 10)]


def _weighted(pairs):
    values, weights = zip(*pairs)
    return random.choices(values, weights=weights, k=1)[0]


def mutate(seq, k):
    """Recopie `seq` en changeant k bases au hasard."""
    if k <= 0:
        return seq
    out = list(seq)
    for i in random.sample(range(len(seq)), min(k, len(seq))):
        out[i] = random.choice([b for b in "ATCG" if b != out[i]])
    return "".join(out)


def k_max(qlen):
    """Le plus grand nombre de mésappariements qui reste présentable.

    Deux bornes, la plus basse gagne :
      * le seuil de confiance de l'application (MIN_IDENTITY) ;
      * le seuil de E-value de blastn — un alignement moins bon ne serait
        tout simplement pas rapporté, et l'afficher sonnerait faux.

    Sur 24 pb la seconde est de loin la plus stricte : elle autorise 2
    mésappariements (E = 0,067), pas 3 (E = 16).
    """
    conf = int(qlen * (100.0 - blast_search.MIN_IDENTITY) / 100.0)
    best = 0
    for k in range(0, conf + 1):
        _bits, e = score_and_evalue(qlen, qlen - k, k)
        if e > blast_search.EVALUE:
            break
        best = k
    return best


def make_hit(sequence, organism, replicons, mismatches):
    """Un hit complet, de la même forme que ceux de blast_search._parse()."""
    qlen = len(sequence)
    matches = qlen - mismatches
    subject = mutate(sequence, mismatches)
    bits, evalue = score_and_evalue(qlen, matches, mismatches)

    accession, locus, title, slen = pick_replicon(organism, replicons)
    # Position tirée au hasard sur le réplicon, avec de la marge aux deux
    # bouts pour que l'idéogramme ne colle jamais à une extrémité. Sans
    # longueur connue, une fenêtre modeste : les coordonnées n'apparaissent
    # alors que dans la gouttière de l'alignement.
    span = slen if slen else 20_000
    margin = max(qlen, int(span * 0.02))
    start = random.randint(margin, max(margin + 1, span - margin))
    plus = random.random() < 0.5
    sstart, send = (start, start + qlen - 1) if plus else (
        start + qlen - 1, start)
    if not plus:
        subject = subject.translate(COMPLEMENT)[::-1]
        query_seq = sequence.translate(COMPLEMENT)[::-1]
    else:
        query_seq = sequence

    return {
        "taxonomy_id": organism["taxonomy_id"],
        "accession": accession,
        "percent_identity": round(100.0 * matches / qlen, 1),
        "alignment_length": qlen,
        "mismatches": mismatches,
        "gaps": 0,
        "query_start": 1, "query_end": qlen,
        "subject_start": sstart, "subject_end": send,
        "subject_length": slen,
        "subject_title": title,
        "strand": "moins" if sstart > send else "plus",
        "locus": locus,
        "evalue": evalue,
        "bitscore": bits,
        "query_seq": query_seq,
        "subject_seq": subject,
        "coverage": 100.0,
    }


def load_replicon_map(path="genome_stats.json"):
    """taxid -> {nom de chromosome: accession}.

    C'est l'inverse de l'index de `api.load_replicons()`, qui va de
    l'accession vers le nom : ici on part du chromosome pour en tirer une
    accession réelle.
    """
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for taxid, info in json.load(fh).items():
            if info.get("replicons"):
                out[int(taxid)] = info["replicons"]
    return out


def pick_replicon(organism, replicons):
    """(accession, étiquette, titre, longueur) du support touché.

    Les accessions sont RÉELLES : elles viennent de l'assemblage de référence
    du NCBI via genome_stats.json. La longueur d'un chromosome, elle, est
    estimée — les rapports NCBI donnent la taille du génome et le nombre de
    chromosomes, jamais la longueur de chacun.
    """
    reps = replicons.get(organism["taxonomy_id"]) or {}
    size_mb = organism.get("genome_size_mb")
    total = max(100_000, int((size_mb or 0) * 1e6))
    name = organism["scientific_name"]

    if not reps and not size_mb:
        # Ni chromosomes, ni taille de génome : c'est le cas de 3 des 32
        # brins préparés (amanite, dionée, salamandre). On n'a alors rien de
        # vrai à dire sur le support, alors on n'en dit rien — pas de locus,
        # pas d'idéogramme. L'alignement, lui, s'affiche quand même : c'est
        # de toute façon le morceau que les enfants regardent.
        return "séquence de référence", None, f"{name} reference sequence", None

    if reps:
        # De préférence un autosome : la longueur affichée est une moyenne
        # (génome ÷ nombre de chromosomes), franchement fausse pour un
        # chromosome sexuel, souvent bien plus court.
        count = organism.get("chromosome_count") or len(reps)
        slen = max(50_000, total // max(1, count))
        keys = sorted(n for n in reps if n.isdigit()) or sorted(reps)
        chrom = random.choice(keys)
        return (reps[chrom], f"chromosome {chrom}",
                f"{name} chromosome {chrom}, complete sequence", slen)

    # Aucune carte des chromosomes pour cette espèce — c'est le cas de 9 des
    # 32 brins préparés. On se garde bien de baptiser « chromosome 1 » un
    # génome entier : on situe le hit dans l'assemblage, dont la taille est
    # réelle. Le dessin y gagne même en éloquence : « ta séquence est là,
    # dans 250 Mb d'ADN ».
    acc = organism.get("assembly_accession") or "—"
    return (acc, "génome complet", f"{name} genome assembly {acc}", total)


# --------------------------------------------------------------------------
# Choix de l'organisme et de sa parenté
# --------------------------------------------------------------------------

# Ordre de préférence des rangs pour aller chercher des cousins : du plus
# proche au plus lointain. Ce sont les colonnes de la table organisms.
COUSIN_RANKS = ["genus", "family", "order_name", "class_name", "phylum"]

# Mésappariements ajoutés à mesure qu'on s'éloigne. Le dégradé est le propos
# de l'écran « ses cousins ». Il reste faible parce qu'il l'est réellement :
# sur le brin « Lion », les vrais cousins sortaient à 100 % et 95,5 % — une
# région conservée l'est chez toute la famille. Un dégradé plus spectaculaire
# serait plus joli et faux.
RANK_PENALTY = {"genus": 1, "family": 1, "order_name": 2,
                "class_name": 2, "phylum": 2}

_SHOWABLE = ("image_path IS NOT NULL AND image_path != '' "
             "AND (common_name_fr IS NOT NULL OR common_name_en IS NOT NULL)")


class Simulator:
    """Tire un organisme, lui fabrique un alignement et une parenté."""

    def __init__(self, db_path, replicons, identifier=None,
                 results_dir=None):
        self.db_path = db_path
        self.replicons = replicons
        self.identifier = identifier
        self.results_dir = (results_dir if results_dir is not None
                            else blast_search.RESULTS_DIR)
        self.pool = self._load_pool()

    # -- pool ------------------------------------------------------------

    def _connect(self):
        conn = sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True)
        conn.row_factory = sqlite3.Row
        return conn

    def _load_pool(self):
        """Les taxids tirables au sort.

        On ne tire pas dans les 25 545 organismes : la plupart n'ont ni nom
        commun, ni génome connu, et donneraient un écran vide au moment le
        plus attendu. On restreint à ceux qui ont une photo, un nom que l'on
        peut lire à voix haute, et un assemblage de référence dont on connaît
        les chromosomes — soit tout ce qu'il faut pour remplir la fiche,
        l'idéogramme et la localisation du hit.
        """
        with self._connect() as conn:
            rows = conn.execute(
                f"SELECT taxonomy_id FROM organisms WHERE {_SHOWABLE} "
                "AND genome_size_mb IS NOT NULL "
                "AND chromosome_count IS NOT NULL "
                "AND chromosome_count > 0").fetchall()
            pool = [r["taxonomy_id"] for r in rows
                    if r["taxonomy_id"] in self.replicons]
            if len(pool) < 50:      # genome_stats.json absent : on élargit
                rows = conn.execute(
                    f"SELECT taxonomy_id FROM organisms WHERE {_SHOWABLE}"
                ).fetchall()
                pool = [r["taxonomy_id"] for r in rows]
        return pool

    def __len__(self):
        return len(self.pool)

    # -- organismes ------------------------------------------------------

    def _fiche(self, taxid):
        with self._connect() as conn:
            row = conn.execute(
                "SELECT * FROM organisms WHERE taxonomy_id = ?",
                (taxid,)).fetchone()
        return dict(row) if row else None

    def _cousins(self, ref, wanted=6):
        """De vrais parents de `ref`, du plus proche au plus lointain."""
        found, seen = [], {ref["taxonomy_id"]}
        with self._connect() as conn:
            for rank in COUSIN_RANKS:
                if len(found) >= wanted:
                    break
                value = ref.get(rank)
                if not value:
                    continue
                marks = ",".join("?" * len(seen))
                rows = conn.execute(
                    f"SELECT * FROM organisms WHERE {rank} = ? "
                    f"AND taxonomy_id NOT IN ({marks}) AND {_SHOWABLE} "
                    "ORDER BY RANDOM() LIMIT 2",
                    [value, *seen]).fetchall()
                for r in rows:
                    o = dict(r)
                    seen.add(o["taxonomy_id"])
                    found.append((rank, o))
                    if len(found) >= wanted:
                        break
        return found

    def choose(self, sequence):
        """(fiche, raison) de l'organisme retenu pour cette séquence.

        Un brin préparé rend toujours son propre organisme : ce sont les
        exemples imprimés que l'on donne aux enfants, et ils doivent tomber
        juste à tous les coups.
        """
        if self.identifier:
            known = self.identifier.identify(sequence)
            if known:
                fiche = self._fiche(known)
                if fiche:
                    return fiche, "brin préparé"
        # Le tirage dépend de la séquence : recoller deux fois la même
        # séquence redonne le même organisme. Un enfant qui recommence pour
        # montrer à ses parents ne doit pas obtenir une autre espèce.
        rng = random.Random(sequence)
        for _ in range(10):
            fiche = self._fiche(rng.choice(self.pool))
            if fiche:
                return fiche, "tirage au sort"
        return None, "pool vide"

    # -- recherche -------------------------------------------------------

    async def search(self, sequence, workdir="/tmp", archive_id=None):
        """Même signature que blast_search.search() : (hits, archive)."""
        return self.search_sync(sequence, archive_id=archive_id)

    def search_sync(self, sequence, archive_id=None):
        organism, reason = self.choose(sequence)
        if not organism:
            return [], None

        qlen = len(sequence)
        cap = k_max(qlen)
        # On laisse toujours un cran libre au-dessus du meilleur hit, pour
        # que les cousins puissent être strictement moins bons que lui.
        top = min(_weighted(TOP_MISMATCHES), max(0, cap - 1))
        # Un brin préparé sort à 100 %. C'est la séquence qu'on tend à
        # l'enfant sur un carton : elle doit tomber juste, et joliment.
        if reason == "brin préparé":
            top = 0

        best = make_hit(sequence, organism, self.replicons, top)
        cousins = [make_hit(sequence, cousin, self.replicons,
                            min(top + RANK_PENALTY[rank], cap))
                   for rank, cousin in self._cousins(organism)]
        cousins.sort(key=lambda h: (h["evalue"], -h["bitscore"]))
        # L'organisme retenu passe devant, quoi qu'il arrive : l'application
        # prend le premier hit sûr de la liste, et ce doit être celui-là.
        hits = [best] + cousins

        archive = self._archive(archive_id, sequence, hits, reason)
        return hits, archive

    # -- trace -----------------------------------------------------------

    def _archive(self, archive_id, sequence, hits, reason):
        """Trace vérifiable, en tête de laquelle il est écrit que c'est
        une simulation. C'est la garantie qu'on ne confondra jamais, plus
        tard, une journée simulée avec une journée de vrai BLAST."""
        if not archive_id or not self.results_dir:
            return None
        try:
            os.makedirs(self.results_dir, exist_ok=True)
            stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            path = os.path.join(self.results_dir,
                                f"{stamp}_{archive_id[:8]}_SIMULATION.txt")
            ell, m, n = length_adjustment(len(sequence))
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("# ===================================\n")
                fh.write("# ANALYSE SIMULEE — AUCUN BLAST REEL\n")
                fh.write("# ===================================\n")
                fh.write(f"# date     : "
                         f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                fh.write(f"# job      : {archive_id}\n")
                fh.write(f"# requete  : {sequence} ({len(sequence)} pb)\n")
                fh.write(f"# organisme: {reason}\n")
                fh.write(f"# banque   : {DB_NAME} — {DB_SEQUENCES} sequences, "
                         f"{DB_LETTERS} lettres (annoncee, non interrogee)\n")
                fh.write(f"# Karlin-Altschul : lambda {LAMBDA}, K {K}, H {H}, "
                         f"reward {REWARD}, penalty {PENALTY}\n")
                fh.write(f"# espace de recherche effectif : {int(m * n)} "
                         f"(ajustement {int(ell)})\n")
                fh.write("\n===== Hits fabriques =====\n")
                fh.write("# taxid\taccession\tpident\tlength\tmismatch\t"
                         "sstart\tsend\tevalue\tbitscore\ttitre\n")
                for h in hits:
                    fh.write("\t".join(str(x) for x in (
                        h["taxonomy_id"], h["accession"],
                        h["percent_identity"], h["alignment_length"],
                        h["mismatches"], h["subject_start"], h["subject_end"],
                        f"{h['evalue']:.3g}", h["bitscore"],
                        h["subject_title"])) + "\n")
                fh.write("\n===== Alignement du meilleur hit =====\n")
                best = hits[0]
                fh.write(f"Query  {best['query_start']:>6}  "
                         f"{best['query_seq']}  {best['query_end']}\n")
                fh.write(f"       {'':>6}  "
                         f"{blast_search.midline(best['query_seq'], best['subject_seq'])}\n")
                fh.write(f"Sbjct  {best['subject_start']:>6}  "
                         f"{best['subject_seq']}  {best['subject_end']}\n")
            return path
        except OSError:
            return None


# --------------------------------------------------------------------------

def main():
    if len(sys.argv) != 2:
        print(__doc__.strip())
        return 2
    from identify import Identifier, normalize

    seq = normalize(sys.argv[1])
    sim = Simulator("biodetective.db", load_replicon_map(), Identifier())
    print(f"Pool de {len(sim)} organismes tirables")
    hits, archive = sim.search_sync(seq, archive_id="cli-test")
    for h in hits:
        print(f"  taxid {h['taxonomy_id']:>8}  {h['accession']:<16} "
              f"{h['percent_identity']:>5} %  E {h['evalue']:.2g}  "
              f"score {h['bitscore']}  {h['subject_title'][:60]}")
    if archive:
        print(f"Trace : {archive}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
