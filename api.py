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
from identify import Identifier, SequenceError, normalize

# --------------------------------------------------------------------------
# Réglages — les trois premiers sont ceux qu'on retouche le jour de la démo
# --------------------------------------------------------------------------

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
ANALYSIS_TIERS = [
    (55, (2.5, 4.0), "rapide"),
    (30, (4.5, 7.0), "normale"),
    (12, (7.5, 10.0), "approfondie"),
    (3, (10.5, 13.0), "très longue"),
]

# Moyenne ~4,9 s. Fixer une durée unique se fait sans toucher au code, utile
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

    if not blast_search.blast_available():
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
        await asyncio.sleep(job["planned_duration"])

        hits = await blast_search.search(sequence)

        # Une fiche par taxon touché, en une seule requête SQL.
        fiches = organisms_by_taxid([h["taxonomy_id"] for h in hits[:10]])

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
                "accession": best["accession"],
                "percent_identity": best["percent_identity"],
                "evalue": best["evalue"],
                "bitscore": best["bitscore"],
                "mismatches": best["mismatches"],
                "gaps": best["gaps"],
            }
            job.update(status="completed", matched=True,
                       organism=fiches[best["taxonomy_id"]])
            log.info("[%s] %s -> %s  id %.1f%%  E %.2g", job_id[:8], sequence,
                     fiches[best["taxonomy_id"]]["scientific_name"],
                     best["percent_identity"], best["evalue"])
        else:
            job.update(status="completed", matched=False, organism=None)
            log.info("[%s] %s -> aucun hit sûr (%d hits bruts)",
                     job_id[:8], sequence, len(hits))
    except Exception as exc:                      # noqa: BLE001
        # Une exception non rattrapée laisserait le job en "running" et le
        # frontend tournerait indéfiniment. Mieux vaut une erreur affichée.
        job.update(status="error", error_message="Erreur interne pendant l'analyse.")
        log.exception("[%s] échec de l'analyse : %s", job_id[:8], exc)
    finally:
        job["analysis_time"] = round(time.time() - job["started_at"], 2)


def job_payload(job):
    return {
        "job_id": job["job_id"],
        "status": job["status"],
        "matched": job.get("matched"),
        "organism": job.get("organism"),
        "sequence": job["sequence"],
        "planned_duration": round(job.get("planned_duration", 0), 1),
        "blast": job.get("blast"),
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
