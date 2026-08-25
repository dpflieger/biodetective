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

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from identify import Identifier, SequenceError, normalize

# --------------------------------------------------------------------------
# Réglages — les trois premiers sont ceux qu'on retouche le jour de la démo
# --------------------------------------------------------------------------

# La recherche en table est instantanée. Ce délai EST la mise en scène :
# il laisse le temps aux images de défiler et à l'enfant de s'installer.
# Trop court, l'animation n'existe pas ; trop long, un enfant de 8 ans décroche.
ANALYSIS_DELAY_SECONDS = 4.0

# Taille du pool d'images tiré au démarrage pour l'animation de recherche.
IMAGE_POOL_SIZE = 400

# Au-delà, les jobs terminés les plus anciens sont oubliés. Sur une journée
# entière de démonstration, sans ça le dictionnaire grossit indéfiniment.
MAX_JOBS = 200

DB_PATH = os.environ.get("BIODETECTIVE_DB", "biodetective.db")
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

    log.info("Base    : %s (%d organismes)", DB_PATH, state.organism_count)
    log.info("Démo    : %d séquences connues", len(state.identifier))
    log.info("Images  : pool de %d en mémoire", len(state.image_pool))
    log.info("Délai   : %.1f s par analyse", ANALYSIS_DELAY_SECONDS)
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


# --------------------------------------------------------------------------
# Modèles
# --------------------------------------------------------------------------

class AnalyzeRequest(BaseModel):
    sequence: str = Field(..., description="Séquence ADN issue du Brickopore")


# --------------------------------------------------------------------------
# Routes de consultation
# --------------------------------------------------------------------------

@app.get("/")
def root():
    """État de l'API — sert de test de vie au frontend."""
    return {
        "app": "BioDetective",
        "status": "ok",
        "organisms": state.organism_count,
        "known_sequences": len(state.identifier) if state.identifier else 0,
        "analysis_delay_seconds": ANALYSIS_DELAY_SECONDS,
        "jobs_in_memory": len(state.jobs),
    }


@app.get("/stats")
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


@app.get("/organism/{taxonomy_id}")
def organism(taxonomy_id: int):
    """Fiche complète d'un organisme."""
    with connect() as conn:
        row = conn.execute("SELECT * FROM organisms WHERE taxonomy_id=?",
                           (taxonomy_id,)).fetchone()
    if row is None:
        raise HTTPException(404, f"Aucun organisme avec le taxid {taxonomy_id}")
    return organism_dict(row)


@app.get("/random-images")
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
        await asyncio.sleep(ANALYSIS_DELAY_SECONDS)

        taxid = state.identifier.identify(sequence)
        if taxid is None:
            job.update(status="completed", matched=False, organism=None)
            log.info("[%s] séquence inconnue : %s", job_id[:8], sequence)
        else:
            with connect() as conn:
                row = conn.execute(
                    "SELECT * FROM organisms WHERE taxonomy_id=?",
                    (taxid,)).fetchone()
            if row is None:
                # sequences.json pointe vers un taxid absent : c'est une vraie
                # panne de configuration, pas une séquence inconnue.
                # validate_sequences.py existe pour l'éviter.
                job.update(
                    status="error",
                    error_message=f"L'organisme {taxid} est introuvable en base.")
                log.error("[%s] taxid %s absent — lancer validate_sequences.py",
                          job_id[:8], taxid)
            else:
                job.update(status="completed", matched=True,
                           organism=organism_dict(row))
                log.info("[%s] %s -> %s", job_id[:8], sequence,
                         row["scientific_name"])
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
        "analysis_time": job.get("analysis_time"),
        "error_message": job.get("error_message"),
    }


@app.post("/analyze")
async def analyze(req: AnalyzeRequest):
    """Lance une analyse et retourne immédiatement un job_id à interroger."""
    try:
        sequence = normalize(req.sequence)
    except SequenceError as exc:
        # 400 : la séquence est inutilisable. À ne pas confondre avec une
        # séquence valide mais inconnue, qui répond 200 avec matched=false.
        raise HTTPException(400, str(exc))

    job_id = str(uuid.uuid4())
    state.jobs[job_id] = {
        "job_id": job_id,
        "status": "running",
        "sequence": sequence,
        "created_at": time.time(),
        "started_at": time.time(),
    }
    purge_jobs()
    asyncio.create_task(run_analysis(job_id, sequence))
    return job_payload(state.jobs[job_id])


@app.get("/analyze/{job_id}")
def analyze_result(job_id: str):
    """Résultat d'une analyse. Interrogé toutes les 500 ms par le frontend."""
    job = state.jobs.get(job_id)
    if job is None:
        raise HTTPException(404, "Analyse inconnue ou déjà oubliée.")
    return job_payload(job)


@app.delete("/analyze/{job_id}")
def analyze_delete(job_id: str):
    if state.jobs.pop(job_id, None) is None:
        raise HTTPException(404, "Analyse inconnue.")
    return {"deleted": job_id}


@app.post("/reload")
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


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host=HOST, port=PORT, log_level="warning")
