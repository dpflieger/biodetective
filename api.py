#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
api.py — Backend FastAPI de BioDetective.

    python3 api.py
    → http://localhost:8000        état de l'API
    → http://localhost:8000/docs   documentation interactive

Pas de BLAST : l'identification passe par identify.py (voir CLAUDE.md).
Les routes s'appellent /analyze et non /blast pour ne pas mentir sur ce
qu'elles font.
"""

import asyncio
import json
import logging
import os
import random
import sqlite3
import sys
import time
import uuid
from contextlib import asynccontextmanager
from typing import Optional

from fastapi import APIRouter, FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

import blast_search
import remote_blast
import simulate
from identify import Identifier, SequenceError, normalize

# --------------------------------------------------------------------------
# Réglages — les trois premiers sont ceux qu'on retouche le jour de la démo
# --------------------------------------------------------------------------

# MODE SIMULATION. Tant que la banque BLAST définitive n'est pas construite,
# chaque séquence rend un organisme tiré au sort (voir simulate.py) : la
# démonstration tourne de bout en bout et la durée de recherche est
# entièrement sous notre contrôle, puisque plus rien ne cherche vraiment.
#
#     BIODETECTIVE_SIMULATE=0 python3 api.py      # vrai blastn
SIMULATE = os.environ.get("BIODETECTIVE_SIMULATE", "1").lower() not in (
    "0", "false", "non", "no", "")

# La recherche en table est instantanée. Cette attente EST la mise en scène :
# elle laisse le temps aux images de défiler et à l'enfant de s'installer.
#
# La durée est tirée au hasard pour que deux enfants qui se suivent ne voient
# pas exactement la même chose, et pour qu'une « analyse difficile » arrive de
# temps en temps.
#
# ⚠️ Le tirage est INDÉPENDANT du résultat, et il a lieu avant même la
# recherche. Si les analyses longues aboutissaient plus souvent à une séquence
# inconnue, l'opérateur — puis les enfants — apprendraient à lire la réponse
# avant l'écran de résultat, et toute la mise en scène tomberait.
#
# (poids, (durée mini, durée maxi), étiquette de journal)
#
# Plafond : 10 s. Au-delà, l'enfant décroche et la file d'attente s'allonge —
# et en mode simulation l'attente n'a plus aucune justification technique, ce
# n'est plus que de la mise en scène.
ANALYSIS_TIERS = [
    (35, (4.5, 6.0), "rapide"),
    (35, (6.0, 7.5), "normale"),
    (22, (7.5, 9.0), "approfondie"),
    (8, (9.0, 10.0), "très longue"),
]

# Moyenne ~6,8 s. Fixer une durée unique se fait sans toucher au code, utile
# si la file d'attente s'allonge en pleine journée :
#     BIODETECTIVE_FIXED_DELAY=3 python3 api.py
_fixed = os.environ.get("BIODETECTIVE_FIXED_DELAY")
ANALYSIS_FIXED_SECONDS = float(_fixed) if _fixed else None

# Taille du pool d'images tiré au démarrage pour l'animation de recherche.
IMAGE_POOL_SIZE = 400

# Au-delà, les jobs terminés les plus anciens sont oubliés. Sur une journée
# entière de démonstration, sans ça le dictionnaire grossit indéfiniment.
MAX_JOBS = 200

DB_PATH = os.environ.get("BIODETECTIVE_DB", "biodetective.db")
GENOME_STATS = "genome_stats.json"
TREE_STATS = "tree_stats.json"

# Historique des analyses. Persisté sur disque : une journée de démonstration
# est longue et un redémarrage du backend ne doit pas l'effacer.
HISTORY_FILE = os.environ.get("BIODETECTIVE_HISTORY", "history.json")
HISTORY_MAX = 40
# Frontend compilé par « npm run build ». Absent = mode développement,
# où react-scripts sert l'interface sur le port 3000.
BUILD_DIR = os.environ.get("BIODETECTIVE_BUILD", "build")
HOST = "0.0.0.0"
PORT = 8000

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("biodetective.log", encoding="utf-8"),
    ],
)
log = logging.getLogger("api")


# --------------------------------------------------------------------------
# Accès base
# --------------------------------------------------------------------------

def connect():
    """Connexion SQLite en lecture seule.

    Une connexion par requête : SQLite n'aime pas le partage entre threads,
    et le coût d'ouverture est négligeable sur un fichier local de 11 Mo.
    """
    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def organism_dict(row):
    """Ligne SQLite -> dict prêt pour le frontend.

    Les champs descriptifs vides sont retournés à None plutôt qu'en chaîne
    vide, pour que React puisse masquer la ligne d'un simple test.
    """
    if row is None:
        return None
    d = {k: (row[k] if row[k] not in ("",) else None) for k in row.keys()}
    return {
        "taxonomy_id": d["taxonomy_id"],
        "scientific_name": d["scientific_name"],
        "common_name_fr": d["common_name_fr"],
        "common_name_en": d["common_name_en"],
        # Nom à afficher en grand : le français s'il existe, l'anglais sinon,
        # et à défaut le nom scientifique. Le frontend n'a pas à arbitrer.
        "display_name": (d["common_name_fr"] or d["common_name_en"]
                         or d["scientific_name"]),
        "organism_type": d["organism_type"],
        "kingdom": d["kingdom"],
        "phylum": d["phylum"],
        "class_name": d["class_name"],
        "order_name": d["order_name"],
        "family": d["family"],
        "genus": d["genus"],
        "species": d["species"],
        "genome_size_mb": d["genome_size_mb"],
        "chromosome_count": d["chromosome_count"],
        "assembly_accession": d["assembly_accession"],
        "size_info": d["size_info"],
        "habitat": d["habitat"],
        "danger_level": d["danger_level"],
        "discovery_info": d["discovery_info"],
        "fun_facts": d["fun_facts"],
        "ecological_role": d["ecological_role"],
        "human_utility": d["human_utility"],
        "image": {
            "url": d["image_path"],
            "source": d["image_source"],
            "license": d["image_license"],
            "attribution": d["image_attribution"],
        } if d["image_path"] else None,
    }


# --------------------------------------------------------------------------
# État de l'application
# --------------------------------------------------------------------------

class State:
    identifier: Optional[Identifier] = None
    image_pool: list = []
    jobs: dict = {}
    organism_count: int = 0
    history: list = []
    # taxid -> {accession sans version: nom du chromosome}
    replicons: dict = {}
    simulator: Optional["simulate.Simulator"] = None


state = State()


def organisms_by_taxid(taxids):
    """Fiches de plusieurs organismes en une requête."""
    if not taxids:
        return {}
    marks = ",".join("?" * len(taxids))
    with connect() as conn:
        rows = conn.execute(
            f"SELECT * FROM organisms WHERE taxonomy_id IN ({marks})",
            list(taxids)).fetchall()
    return {r["taxonomy_id"]: organism_dict(r) for r in rows}


def strip_version(acc):
    """NC_003070.9 -> NC_003070. Les versions divergent entre sources."""
    return acc.rsplit(".", 1)[0] if "." in acc else acc


def load_replicons(path):
    """Index accession -> nom de chromosome, par espèce.

    Permet de nommer « chromosome 1 » le chromosome touché de façon fiable,
    à partir de son accession, plutôt qu'en analysant un texte libre.
    """
    if not os.path.exists(path):
        log.warning("%s absent : pas d'idéogramme", path)
        return {}
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    out = {}
    for taxid, v in data.items():
        reps = v.get("replicons") or {}
        if reps:
            out[int(taxid)] = {strip_version(a): n for n, a in reps.items()}
    log.info("%s : réplicons connus pour %d espèces", path, len(out))
    return out


def load_history():
    if not os.path.exists(HISTORY_FILE):
        return []
    try:
        with open(HISTORY_FILE, encoding="utf-8") as fh:
            data = json.load(fh)
        return data.get("entries", [])[:HISTORY_MAX]
    except (OSError, ValueError) as exc:
        # Un historique corrompu ne doit pas empêcher l'application de
        # démarrer devant une file d'enfants.
        log.warning("Historique illisible (%s), on repart de zéro", exc)
        return []


def save_history():
    """Écriture atomique : un arrêt brutal ne laisse pas un fichier tronqué."""
    tmp = HISTORY_FILE + ".tmp"
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump({"entries": state.history}, fh, ensure_ascii=False)
        os.replace(tmp, HISTORY_FILE)
    except OSError as exc:
        log.warning("Historique non enregistré : %s", exc)


def record_history(job):
    """Ajoute une analyse terminée en tête de l'historique."""
    o = job.get("organism")
    entry = {
        "at": time.time(),
        "sequence": job["sequence"],
        "matched": bool(job.get("matched")),
        "analysis_time": job.get("analysis_time"),
        "organism": o,
        "blast": job.get("blast"),
        "blast_archive": job.get("blast_archive"),
        "simulated": job.get("simulated", False),
        "remote": job.get("remote"),
    }
    state.history.insert(0, entry)
    del state.history[HISTORY_MAX:]
    save_history()


def load_image_pool():
    """Tire une fois pour toutes un pool d'images pour l'animation.

    ORDER BY RANDOM() sur 25 545 lignes à chaque analyse serait un gâchis :
    on tire IMAGE_POOL_SIZE images au démarrage et on pioche dedans ensuite.
    """
    with connect() as conn:
        rows = conn.execute(
            "SELECT taxonomy_id, scientific_name, common_name_fr, "
            "       common_name_en, image_path, organism_type "
            "FROM organisms WHERE image_path IS NOT NULL AND image_path != '' "
            "ORDER BY RANDOM() LIMIT ?", (IMAGE_POOL_SIZE,)).fetchall()
    return [{
        "taxonomy_id": r["taxonomy_id"],
        "scientific_name": r["scientific_name"],
        "display_name": (r["common_name_fr"] or r["common_name_en"]
                         or r["scientific_name"]),
        "organism_type": r["organism_type"],
        "url": r["image_path"],
    } for r in rows]


@asynccontextmanager
async def lifespan(app: FastAPI):
    if not os.path.exists(DB_PATH):
        log.error("%s introuvable. Lancer d'abord : python3 data_pipeline.py",
                  DB_PATH)
        raise RuntimeError(f"{DB_PATH} introuvable")

    with connect() as conn:
        state.organism_count = conn.execute(
            "SELECT count(*) FROM organisms").fetchone()[0]

    state.identifier = Identifier()
    state.image_pool = load_image_pool()
    state.history = load_history()
    state.replicons = load_replicons(GENOME_STATS)

    if SIMULATE:
        state.simulator = simulate.Simulator(
            DB_PATH, simulate.load_replicon_map(GENOME_STATS),
            state.identifier)
        log.warning("=" * 62)
        log.warning("MODE SIMULATION — aucun blastn n'est exécuté.")
        log.warning("Chaque séquence rend un organisme tiré au sort parmi "
                    "%d.", len(state.simulator))
        log.warning("Les %d brins préparés rendent bien LEUR organisme.",
                    len(state.identifier))
        log.warning("Repasser au vrai BLAST : BIODETECTIVE_SIMULATE=0")
        log.warning("=" * 62)
    elif not blast_search.blast_available():
        log.error("blastn introuvable : aucune analyse ne pourra aboutir.")
    elif not blast_search.db_available():
        log.error("Banque BLAST absente (%s). Voir CLAUDE.md.",
                  blast_search.BLAST_DB)
    else:
        log.info("BLAST   : %s", blast_search.BLASTN)
        log.info("Banque  : %s", blast_search.BLAST_DB)

    log.info("Base    : %s (%d organismes)", DB_PATH, state.organism_count)
    log.info("Démo    : %d séquences connues", len(state.identifier))
    log.info("Images  : pool de %d en mémoire", len(state.image_pool))
    log.info("Histori.: %d analyses reprises de %s",
             len(state.history), HISTORY_FILE)
    if ANALYSIS_FIXED_SECONDS is not None:
        log.info("Durée   : fixée à %.1f s", ANALYSIS_FIXED_SECONDS)
    else:
        lo = min(t[1][0] for t in ANALYSIS_TIERS)
        hi = max(t[1][1] for t in ANALYSIS_TIERS)
        avg = sum(w * (a + b) / 2 for w, (a, b), _ in ANALYSIS_TIERS) / sum(
            t[0] for t in ANALYSIS_TIERS)
        log.info("Durée   : %.1f à %.1f s, moyenne %.1f s", lo, hi, avg)
    log.info("Prêt sur http://localhost:%d  (docs : /docs)", PORT)
    yield
    log.info("Arrêt de l'API")


app = FastAPI(
    title="BioDetective",
    description="Identification d'organismes pour la Fête de la Science (IBMP)",
    version="1.0.0",
    lifespan=lifespan,
)

# Le proxy de react-scripts suffit en développement, mais si le frontend est
# un jour servi depuis un autre port ou une autre machine de la salle,
# l'absence de CORS se manifeste par un écran blanc difficile à diagnostiquer.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Toutes les routes de l'API vivent sous /api pour laisser la racine au
# frontend compilé. Sans ce préfixe, GET / ne pourrait pas être à la fois
# l'état de l'API et la page d'accueil.
api = APIRouter(prefix="/api")


# --------------------------------------------------------------------------
# Modèles
# --------------------------------------------------------------------------

class AnalyzeRequest(BaseModel):
    sequence: str = Field(..., description="Séquence ADN issue du Brickopore")


# --------------------------------------------------------------------------
# Routes de consultation
# --------------------------------------------------------------------------

@api.get("/")
def root():
    """État de l'API — sert de test de vie au frontend."""
    return {
        "app": "BioDetective",
        "status": "ok",
        "organisms": state.organism_count,
        "known_sequences": len(state.identifier) if state.identifier else 0,
        "analysis_fixed_seconds": ANALYSIS_FIXED_SECONDS,
        "simulated": SIMULATE,
        "jobs_in_memory": len(state.jobs),
    }


@api.get("/stats")
def stats():
    """Statistiques de la base, pour un écran d'accueil ou un contrôle rapide."""
    with connect() as conn:
        by_type = conn.execute(
            "SELECT organism_type, count(*) n FROM organisms "
            "GROUP BY organism_type ORDER BY n DESC").fetchall()
        named_fr = conn.execute(
            "SELECT count(*) FROM organisms "
            "WHERE common_name_fr IS NOT NULL").fetchone()[0]
        named_en = conn.execute(
            "SELECT count(*) FROM organisms "
            "WHERE common_name_en IS NOT NULL").fetchone()[0]
    return {
        "organisms": state.organism_count,
        "with_french_name": named_fr,
        "with_english_name": named_en,
        "known_sequences": len(state.identifier) if state.identifier else 0,
        "simulated": SIMULATE,
        "simulation_pool": len(state.simulator) if state.simulator else 0,
        "by_type": [{"type": r["organism_type"], "count": r["n"]}
                    for r in by_type],
    }


@api.get("/organism/{taxonomy_id}")
def organism(taxonomy_id: int):
    """Fiche complète d'un organisme."""
    with connect() as conn:
        row = conn.execute("SELECT * FROM organisms WHERE taxonomy_id=?",
                           (taxonomy_id,)).fetchone()
    if row is None:
        raise HTTPException(404, f"Aucun organisme avec le taxid {taxonomy_id}")
    return organism_dict(row)


@api.get("/random-images")
def random_images(count: int = 8):
    """Images aléatoires pour l'animation de l'écran de recherche.

    Puisées dans le pool chargé au démarrage, jamais en base.
    """
    count = max(1, min(count, 50))
    if not state.image_pool:
        return {"images": []}
    k = min(count, len(state.image_pool))
    return {"images": random.sample(state.image_pool, k)}


# --------------------------------------------------------------------------
# Analyse
# --------------------------------------------------------------------------

def draw_duration():
    """Tire la durée d'une analyse. Sans lien avec la séquence ni le résultat."""
    if ANALYSIS_FIXED_SECONDS is not None:
        return float(ANALYSIS_FIXED_SECONDS), "fixe"
    weights = [t[0] for t in ANALYSIS_TIERS]
    low, high = 0, 0
    tier = random.choices(ANALYSIS_TIERS, weights=weights, k=1)[0]
    _w, (low, high), label = tier
    return random.uniform(low, high), label


# Du plus précis au plus large. Le premier rang où deux organismes
# coïncident dit à quel point ils sont cousins.
RANKS = [
    ("genus", "genre"),
    ("family", "famille"),
    ("order_name", "ordre"),
    ("class_name", "classe"),
    ("phylum", "embranchement"),
    ("kingdom", "règne"),
]


def shared_rank(a, b):
    """(clé, étiquette, valeur) du rang commun le plus précis. None sinon."""
    for key, label in RANKS:
        va, vb = a.get(key), b.get(key)
        if va and vb and va == vb:
            return key, label, va
    return None


def build_relatives(hits, fiches, best):
    """Les autres hits, situés par rapport à l'organisme identifié.

    C'est la lecture honnête d'une liste de hits BLAST : ce ne sont pas des
    « moins bonnes réponses », ce sont les parents de l'organisme trouvé.
    Sur un brin d'Arabidopsis, les suivants sont moutarde, chou et colza —
    toute la famille des Brassicacées.
    """
    ref = fiches.get(best["taxonomy_id"])
    if not ref:
        return []
    out, seen = [], {best["taxonomy_id"]}
    for h in hits:
        taxid = h["taxonomy_id"]
        if taxid in seen:
            continue
        o = fiches.get(taxid)
        if not o:
            continue
        sr = shared_rank(ref, o)
        if not sr:
            continue
        seen.add(taxid)
        key, label, value = sr
        out.append({
            "taxonomy_id": taxid,
            "display_name": o["display_name"],
            "scientific_name": o["scientific_name"],
            "organism_type": o["organism_type"],
            "image": (o.get("image") or {}).get("url"),
            "rank_key": key,
            "rank_label": label,
            "rank_value": value,
            "percent_identity": h["percent_identity"],
            "evalue": h["evalue"],
        })
        if len(out) >= 6:
            break
    # Du plus proche au plus lointain : l'ordre des rangs, pas celui de BLAST.
    order = {k: i for i, (k, _l) in enumerate(RANKS)}
    out.sort(key=lambda r: (order[r["rank_key"]], r["evalue"]))
    return out


def build_ideogram(hit, organism):
    """De quoi dessiner le chromosome touché et y placer le hit.

    On ne dessine que le chromosome atteint, à sa longueur réelle : c'est la
    seule dont BLAST nous donne la taille exacte (slen). Représenter les
    autres supposerait des longueurs qu'on n'a pas, et un caryotype inventé
    vaudrait moins que pas de dessin du tout.
    """
    length = hit.get("subject_length")
    if not length or length < 1000:
        return None      # un marqueur de 800 pb ne se dessine pas

    taxid = hit["taxonomy_id"]
    name = state.replicons.get(taxid, {}).get(
        strip_version(hit.get("accession", "")))
    # À défaut du nom officiel, l'étiquette tirée du titre.
    label = (f"Chromosome {name}" if name
             else (hit.get("locus") or "Séquence de référence").capitalize())

    start, end = hit["subject_start"], hit["subject_end"]
    pos = min(start, end)
    return {
        "label": label,
        "accession": hit.get("accession"),
        "length": length,
        "position": pos,
        "fraction": round(pos / length, 6),
        "strand": hit.get("strand"),
        "chromosome_count": organism.get("chromosome_count"),
        "is_chromosome": bool(name),
    }


def purge_jobs():
    """Oublie les jobs terminés les plus anciens au-delà de MAX_JOBS."""
    if len(state.jobs) <= MAX_JOBS:
        return
    done = sorted(
        (j for j in state.jobs.values() if j["status"] != "running"),
        key=lambda j: j["created_at"])
    for job in done[: len(state.jobs) - MAX_JOBS]:
        state.jobs.pop(job["job_id"], None)


async def run_analysis(job_id: str, sequence: str):
    """Tâche de fond : attend, puis résout la séquence.

    L'attente précède volontairement la recherche pour que la durée perçue
    soit constante quel que soit le résultat — sinon une correspondance
    trouvée reviendrait plus vite qu'une inconnue, ce qui se remarque.
    """
    job = state.jobs[job_id]
    try:
        # L'attente et la recherche courent ENSEMBLE. Auparavant la seconde
        # s'ajoutait à la première ; désormais elle se cache derrière, tant
        # qu'elle reste plus courte — ce qui laisse de la marge pour une
        # banque plus grosse sans rallonger la démonstration.
        #
        # La durée perçue reste donc celle du tirage, indépendante du
        # résultat : c'est la propriété à préserver.
        engine = state.simulator if SIMULATE else blast_search
        (hits, archive), _ = await asyncio.gather(
            engine.search(sequence, archive_id=job_id),
            asyncio.sleep(job["planned_duration"]),
        )
        job["blast_archive"] = archive
        job["simulated"] = SIMULATE

        # Une fiche par taxon touché, en une seule requête SQL.
        fiches = organisms_by_taxid([h["taxonomy_id"] for h in hits[:30]])

        ranked = []
        for h in hits[:5]:
            o = fiches.get(h["taxonomy_id"])
            ranked.append({
                "taxonomy_id": h["taxonomy_id"],
                "accession": h["accession"],
                "scientific_name": o["scientific_name"] if o else None,
                "display_name": o["display_name"] if o else None,
                "organism_type": o["organism_type"] if o else None,
                "has_image": bool(o and o["image"]),
                "percent_identity": h["percent_identity"],
                "coverage": h["coverage"],
                "evalue": h["evalue"],
                "bitscore": h["bitscore"],
                "confident": blast_search.is_confident(h),
            })

        # Le meilleur hit sûr ET dont on possède la fiche : inutile
        # d'annoncer un organisme qu'on ne saurait pas illustrer.
        best = next((h for h in hits
                     if blast_search.is_confident(h)
                     and fiches.get(h["taxonomy_id"])), None)

        job["blast"] = {
            "hits": ranked,
            "total_hits": len(hits),
            "alignment": None,
            "relatives": [],
            # Le frontend s'en sert pour masquer le rapport brut : un rapport
            # simulé n'a rien à faire sous les yeux d'un visiteur.
            "simulated": SIMULATE,
        }
        if best:
            job["blast"]["alignment"] = {
                "query_seq": best["query_seq"],
                "midline": blast_search.midline(best["query_seq"],
                                                best["subject_seq"]),
                "subject_seq": best["subject_seq"],
                "query_start": best["query_start"],
                "query_end": best["query_end"],
                "subject_start": best["subject_start"],
                "subject_end": best["subject_end"],
                "subject_length": best.get("subject_length"),
                "subject_title": best.get("subject_title"),
                "locus": best.get("locus"),
                "strand": best.get("strand"),
                "ideogram": build_ideogram(best, fiches[best["taxonomy_id"]]),
                "accession": best["accession"],
                "percent_identity": best["percent_identity"],
                "evalue": best["evalue"],
                "bitscore": best["bitscore"],
                "mismatches": best["mismatches"],
                "gaps": best["gaps"],
            }
            job["blast"]["relatives"] = build_relatives(hits, fiches, best)
            job.update(status="completed", matched=True,
                       organism=fiches[best["taxonomy_id"]])
            blast_search.append_verdict(archive,
                f"Organisme retenu : "
                f"{fiches[best['taxonomy_id']]['scientific_name']} "
                f"(taxid {best['taxonomy_id']}, {best['accession']})\n"
                f"Identite {best['percent_identity']} %, "
                f"couverture {best['coverage']} %, E {best['evalue']:.3g}, "
                f"score {best['bitscore']}")
            blast_search.append_index({
                "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                "job": job_id[:8], "sequence": sequence,
                "resultat": fiches[best["taxonomy_id"]]["scientific_name"],
                "taxid": best["taxonomy_id"],
                "identite": best["percent_identity"],
                "couverture": best["coverage"],
                "evalue": f"{best['evalue']:.3g}",
                "score": best["bitscore"],
                "fichier": os.path.basename(archive) if archive else "",
            })
            log.info("[%s]%s %s -> %s  id %.1f%%  E %.2g", job_id[:8],
                     " SIM" if SIMULATE else "", sequence,
                     fiches[best["taxonomy_id"]]["scientific_name"],
                     best["percent_identity"], best["evalue"])
        else:
            job.update(status="completed", matched=False, organism=None)
            blast_search.append_verdict(archive,
                f"Aucun hit juge sur (seuils : identite "
                f"{blast_search.MIN_IDENTITY} %, couverture "
                f"{blast_search.MIN_COVERAGE} %). {len(hits)} hits bruts.")
            blast_search.append_index({
                "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                "job": job_id[:8], "sequence": sequence,
                "resultat": "SEQUENCE INCONNUE",
                "fichier": os.path.basename(archive) if archive else "",
            })
            log.info("[%s] %s -> aucun hit sûr (%d hits bruts)",
                     job_id[:8], sequence, len(hits))
    except Exception as exc:                      # noqa: BLE001
        # Une exception non rattrapée laisserait le job en "running" et le
        # frontend tournerait indéfiniment. Mieux vaut une erreur affichée.
        job.update(status="error", error_message="Erreur interne pendant l'analyse.")
        log.exception("[%s] échec de l'analyse : %s", job_id[:8], exc)
    finally:
        job["analysis_time"] = round(time.time() - job["started_at"], 2)
        # Même une séquence inconnue mérite sa ligne : c'est le cas le plus
        # fréquent, et l'opérateur veut savoir combien d'enfants sont passés.
        if job["status"] != "running":
            record_history(job)


def job_payload(job):
    return {
        "job_id": job["job_id"],
        "status": job["status"],
        "matched": job.get("matched"),
        "organism": job.get("organism"),
        "sequence": job["sequence"],
        "planned_duration": round(job.get("planned_duration", 0), 1),
        "blast": job.get("blast"),
        "blast_archive": job.get("blast_archive"),
        "simulated": job.get("simulated", False),
        "remote": job.get("remote"),
        "analysis_time": job.get("analysis_time"),
        "error_message": job.get("error_message"),
    }


@api.post("/analyze")
async def analyze(req: AnalyzeRequest):
    """Lance une analyse et retourne immédiatement un job_id à interroger."""
    try:
        sequence = normalize(req.sequence)
    except SequenceError as exc:
        # 400 : la séquence est inutilisable. À ne pas confondre avec une
        # séquence valide mais inconnue, qui répond 200 avec matched=false.
        raise HTTPException(400, str(exc))

    duration, tier = draw_duration()
    job_id = str(uuid.uuid4())
    state.jobs[job_id] = {
        "job_id": job_id,
        "status": "running",
        "sequence": sequence,
        "planned_duration": duration,
        "created_at": time.time(),
        "started_at": time.time(),
    }
    log.info("[%s] analyse %s de %.1f s", job_id[:8], tier, duration)
    purge_jobs()
    asyncio.create_task(run_analysis(job_id, sequence))
    return job_payload(state.jobs[job_id])


@api.get("/analyze/{job_id}")
def analyze_result(job_id: str):
    """Résultat d'une analyse. Interrogé toutes les 500 ms par le frontend."""
    job = state.jobs.get(job_id)
    if job is None:
        raise HTTPException(404, "Analyse inconnue ou déjà oubliée.")
    return job_payload(job)


@api.delete("/analyze/{job_id}")
def analyze_delete(job_id: str):
    if state.jobs.pop(job_id, None) is None:
        raise HTTPException(404, "Analyse inconnue.")
    return {"deleted": job_id}


async def run_remote(job_id, sequence):
    """BLAST distant, en tâche de fond. Hors du parcours de démonstration."""
    job = state.jobs[job_id]
    try:
        # urllib est bloquant : sans thread, l'API entière se figerait
        # pendant la trentaine de secondes que dure la recherche.
        hits, rid, elapsed = await asyncio.to_thread(
            remote_blast.search, sequence, "core_nt", None, False)
        fiches = organisms_by_taxid(
            [h["taxonomy_id"] for h in hits[:10] if h.get("taxonomy_id")])
        ranked = []
        for h in hits[:10]:
            o = fiches.get(h.get("taxonomy_id"))
            ranked.append({
                "accession": h["accession"],
                "taxonomy_id": h.get("taxonomy_id"),
                "scientific_name": o["scientific_name"] if o else None,
                "display_name": o["display_name"] if o else None,
                "subject_title": h.get("subject_title"),
                "percent_identity": h["percent_identity"],
                "coverage": h.get("coverage"),
                "evalue": h["evalue"],
                "bitscore": h["bitscore"],
                "known_here": o is not None,
            })
        job.update(status="completed", matched=bool(hits), organism=None,
                   remote={"rid": rid, "elapsed": round(elapsed, 1),
                           "database": "core_nt", "hits": ranked,
                           "total_hits": len(hits)})
        log.info("[%s] distant : %d hits en %.0f s (RID %s)",
                 job_id[:8], len(hits), elapsed, rid)
    except Exception as exc:                          # noqa: BLE001
        job.update(status="error",
                   error_message=f"BLAST distant impossible : {exc}")
        log.warning("[%s] distant en échec : %s", job_id[:8], exc)
    finally:
        job["analysis_time"] = round(time.time() - job["started_at"], 2)


@api.post("/remote-blast")
async def remote_blast_start(req: AnalyzeRequest):
    """Soumet la séquence au BLAST public du NCBI.

    ⚠️ Compter une trentaine de secondes, parfois des minutes : sans commune
    mesure avec l'écran de recherche. Outil de vérification, pas de démo.
    Le résultat se relit avec GET /api/analyze/{job_id}.
    """
    try:
        sequence = normalize(req.sequence)
    except SequenceError as exc:
        raise HTTPException(400, str(exc))

    job_id = str(uuid.uuid4())
    state.jobs[job_id] = {
        "job_id": job_id, "status": "running", "sequence": sequence,
        "planned_duration": 0, "remote": None,
        "created_at": time.time(), "started_at": time.time(),
    }
    purge_jobs()
    asyncio.create_task(run_remote(job_id, sequence))
    log.info("[%s] BLAST distant lancé : %s", job_id[:8], sequence)
    return job_payload(state.jobs[job_id])


@api.get("/blast-report/{job_id}")
def blast_report(job_id: str):
    """Le rapport blastn brut d'une analyse, tel qu'archivé.

    Ce que voit un bioinformaticien quand il demande « et ça donne quoi
    vraiment ? » : la commande, le tableau, les alignements, les paramètres
    de Karlin-Altschul.
    """
    job = state.jobs.get(job_id)
    path = job.get("blast_archive") if job else None
    if not path:
        raise HTTPException(404, "Aucun rapport pour cette analyse.")

    # Le chemin vient de nos propres écritures, mais on vérifie tout de même
    # qu'il reste dans le dossier d'archives : une route qui rend un fichier
    # ne doit jamais pouvoir en rendre un autre.
    root = os.path.realpath(blast_search.RESULTS_DIR or ".")
    full = os.path.realpath(path)
    if not full.startswith(root + os.sep) or not os.path.isfile(full):
        raise HTTPException(404, "Rapport introuvable.")

    with open(full, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    return {"job_id": job_id, "file": os.path.basename(full),
            "bytes": len(text), "report": text}


# L'arbre du vivant tel qu'on le montre au public. Volontairement grossier :
# trois domaines, quelques grands groupes sous les eucaryotes, les virus à
# part. Chaque noeud dit à quoi le reconnaître, en une phrase.
#
# « Procaryote » n'y figure pas comme branche, et c'est délibéré : le mot
# désigne une cellule sans noyau, pas un groupe de parenté. Bactéries et
# archées sont deux domaines distincts, aussi éloignés l'un de l'autre que
# de nous. Le terme est expliqué à part, comme une description.
TREE = {
    "domaines": [
        {
            "id": "bacteria", "nom": "Bactéries", "latin": "Bacteria",
            "match": {"superkingdom": ["Bacteria"]}, "stat": "bacteria",
            "phrase": "Une seule cellule, sans noyau. Les êtres vivants les "
                      "plus nombreux de la planète — il y en a plus dans ta "
                      "bouche que d'humains sur Terre.",
        },
        {
            "id": "archaea", "nom": "Archées", "latin": "Archaea",
            "match": {"superkingdom": ["Archaea"]}, "stat": "archaea",
            "phrase": "Elles ressemblent à des bactéries mais forment une "
                      "branche à part. Beaucoup vivent là où rien d'autre ne "
                      "tient : sources brûlantes, lacs salés.",
        },
        {
            "id": "eukaryota", "nom": "Eucaryotes", "latin": "Eukaryota",
            "match": {"superkingdom": ["Eukaryota"]}, "stat": "eukaryota",
            "phrase": "Leurs cellules rangent l'ADN dans un noyau. Nous en "
                      "faisons partie, avec les animaux, les plantes et les "
                      "champignons.",
            "enfants": [
                {"id": "metazoa", "nom": "Animaux", "latin": "Metazoa",
                 "match": {"kingdom": ["Metazoa"]}, "stat": "metazoa",
                 "phrase": "Ils mangent d'autres êtres vivants et se "
                           "déplacent, au moins un moment de leur vie."},
                {"id": "plantae", "nom": "Plantes", "latin": "Viridiplantae",
                 "match": {"kingdom": ["Viridiplantae"]},
                 "stat": "viridiplantae",
                 "phrase": "Elles fabriquent leur nourriture avec la lumière "
                           "du soleil."},
                {"id": "fungi", "nom": "Champignons", "latin": "Fungi",
                 "match": {"kingdom": ["Fungi"]}, "stat": "fungi",
                 "phrase": "Ni plantes ni animaux : ils digèrent leur "
                           "nourriture autour d'eux, puis l'absorbent."},
                {"id": "protistes", "nom": "Algues et protistes",
                 "latin": "autres eucaryotes",
                 "match": {"kingdom_not": ["Metazoa", "Viridiplantae",
                                           "Fungi"]},
                 "stat": "autres_eucaryotes",
                 "phrase": "Tout le reste : algues rouges et brunes, "
                           "diatomées, amibes. Souvent minuscules, souvent "
                           "oubliés."},
            ],
        },
    ],
    "a_part": {
        "id": "virus", "nom": "Virus", "latin": "Virus", "stat": "viruses",
        "phrase": "Ni tout à fait vivants, ni tout à fait inertes. Ils n'ont "
                  "pas de cellule et doivent en emprunter une pour se "
                  "reproduire.",
    },
    "note": "« Procaryote » veut dire « cellule sans noyau » : cela décrit "
            "les bactéries et les archées, mais ce n'est pas une branche de "
            "l'arbre. Ces deux domaines sont aussi éloignés l'un de l'autre "
            "qu'ils le sont de nous.",
}


def _tree_where(match):
    """Traduit la règle d'un noeud en clause SQL."""
    if "superkingdom" in match:
        vals = match["superkingdom"]
        return (f"superkingdom IN ({','.join('?' * len(vals))})", list(vals))
    if "kingdom" in match:
        vals = match["kingdom"]
        return (f"kingdom IN ({','.join('?' * len(vals))})", list(vals))
    if "kingdom_not" in match:
        vals = match["kingdom_not"]
        return ("superkingdom = 'Eukaryota' AND "
                f"(kingdom IS NULL OR kingdom NOT IN "
                f"({','.join('?' * len(vals))}))", list(vals))
    return ("1=0", [])


def load_tree_stats():
    """Espèces décrites par grand groupe, produites par build_tree_stats.py."""
    if not os.path.exists(TREE_STATS):
        log.warning("%s absent : l'arbre n'affichera que notre collection",
                    TREE_STATS)
        return {}
    with open(TREE_STATS, encoding="utf-8") as fh:
        return json.load(fh)


def _tree_node(conn, node, stats):
    where, args = _tree_where(node["match"])
    n = conn.execute(
        f"SELECT count(*) FROM organisms WHERE {where}", args).fetchone()[0]
    # Une vignette pour incarner le groupe : sans image, un noeud n'est
    # qu'un mot latin de plus.
    row = conn.execute(
        f"SELECT scientific_name, common_name_fr, common_name_en, image_path "
        f"FROM organisms WHERE ({where}) AND image_path IS NOT NULL "
        f"ORDER BY RANDOM() LIMIT 1", args).fetchone()
    out = {k: node[k] for k in ("id", "nom", "latin", "phrase")}
    # « décrites » : le nombre d'espèces connues de la science. « en photo » :
    # ce que nous savons illustrer. L'écart est énorme, et il est parlant.
    out["described"] = stats.get(node.get("stat", ""))
    out["count"] = n
    out["exemple"] = {
        "nom": (row["common_name_fr"] or row["common_name_en"]
                or row["scientific_name"]),
        "scientific_name": row["scientific_name"],
        "image": row["image_path"],
    } if row else None
    return out


@api.get("/tree")
def tree():
    """L'arbre du vivant, avec les effectifs réels de notre collection."""
    stats = load_tree_stats()
    with connect() as conn:
        domaines = []
        for d in TREE["domaines"]:
            node = _tree_node(conn, d, stats)
            node["enfants"] = [_tree_node(conn, c, stats)
                               for c in d.get("enfants", [])]
            domaines.append(node)

        # Les virus se répartissent sur plusieurs realms (-viria) : on les
        # compte par différence plutôt que d'énumérer une liste qui bougera.
        viruses = conn.execute(
            "SELECT count(*) FROM organisms WHERE organism_type "
            "IN ('Virus', 'Viroïde')").fetchone()[0]
        vrow = conn.execute(
            "SELECT scientific_name, common_name_fr, common_name_en, "
            "image_path FROM organisms WHERE organism_type IN "
            "('Virus','Viroïde') AND image_path IS NOT NULL "
            "ORDER BY RANDOM() LIMIT 1").fetchone()
        total = conn.execute("SELECT count(*) FROM organisms").fetchone()[0]

    a_part = dict(TREE["a_part"])
    a_part["described"] = stats.get("viruses")
    a_part["count"] = viruses
    a_part["exemple"] = {
        "nom": (vrow["common_name_fr"] or vrow["common_name_en"]
                or vrow["scientific_name"]),
        "scientific_name": vrow["scientific_name"],
        "image": vrow["image_path"],
    } if vrow else None

    return {"total": total,
            "total_described": stats.get("total_species"),
            "domaines": domaines, "a_part": a_part, "note": TREE["note"]}


@api.get("/history")
def history(limit: int = 12, matched_only: bool = False):
    """Dernières analyses, la plus récente en tête."""
    limit = max(1, min(limit, HISTORY_MAX))
    # L'index est celui de la liste complète : le frontend filtre parfois sur
    # les seules réussites, et un clic doit retrouver la bonne analyse.
    rows = [(i, h) for i, h in enumerate(state.history)
            if h["matched"] or not matched_only]
    out = []
    for i, h in rows[:limit]:
        o = h.get("organism") or {}
        out.append({
            "index": i,
            "at": h["at"],
            "sequence": h["sequence"],
            "matched": h["matched"],
            "analysis_time": h.get("analysis_time"),
            "taxonomy_id": o.get("taxonomy_id"),
            "display_name": o.get("display_name"),
            "scientific_name": o.get("scientific_name"),
            "organism_type": o.get("organism_type"),
            "image": (o.get("image") or {}).get("url"),
            "percent_identity": ((h.get("blast") or {}).get("alignment")
                                 or {}).get("percent_identity"),
        })
    found = sum(1 for h in state.history if h["matched"])
    return {
        "entries": out,
        "total": len(state.history),
        "found": found,
        "unknown": len(state.history) - found,
    }


@api.get("/history/{index}")
def history_entry(index: int):
    """Détail complet d'une analyse passée, pour la réafficher à l'écran."""
    if not 0 <= index < len(state.history):
        raise HTTPException(404, "Analyse absente de l'historique.")
    h = state.history[index]
    return {
        "status": "completed",
        "matched": h["matched"],
        "organism": h.get("organism"),
        "blast": h.get("blast"),
        "sequence": h["sequence"],
        "analysis_time": h.get("analysis_time"),
        "from_history": True,
    }


@api.delete("/history")
def history_clear():
    """Vide l'historique — entre deux groupes, ou en fin de journée."""
    n = len(state.history)
    state.history = []
    save_history()
    log.info("Historique vidé (%d analyses)", n)
    return {"cleared": n}


@api.post("/reload")
def reload_sequences():
    """Recharge sequences.json sans redémarrer l'API.

    Permet de corriger ou d'ajouter une séquence en pleine journée de démo
    sans couper le service.
    """
    try:
        state.identifier = Identifier()
    except (OSError, ValueError) as exc:
        raise HTTPException(500, f"Rechargement impossible : {exc}")
    log.info("sequences.json rechargé : %d séquences", len(state.identifier))
    return {"known_sequences": len(state.identifier)}


app.include_router(api)


# Le montage statique vient APRÈS include_router : FastAPI teste les routes
# dans l'ordre de déclaration, donc /api/... est résolu avant d'atteindre le
# catch-all du frontend.
if os.path.isdir(BUILD_DIR):
    @app.get("/", include_in_schema=False)
    def index():
        return FileResponse(os.path.join(BUILD_DIR, "index.html"))

    app.mount("/", StaticFiles(directory=BUILD_DIR, html=True), name="frontend")
    log.info("Frontend servi depuis %s/", BUILD_DIR)
else:
    @app.get("/", include_in_schema=False)
    def index_missing():
        return {
            "message": f"Frontend non compilé ({BUILD_DIR}/ absent).",
            "action": "Lancer « npm run build », ou « npm start » pour le "
                      "mode développement sur le port 3000.",
            "api": "/api",
            "docs": "/docs",
        }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host=HOST, port=PORT, log_level="warning")
