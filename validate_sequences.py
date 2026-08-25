#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_sequences.py — Vérifie sequences.json avant une démo.

À lancer après toute modification de sequences.json, et une dernière fois
le matin de la Fête de la Science. Un taxid absent ou une image morte ne
se voit sinon que devant l'enfant.

    python3 validate_sequences.py            # contrôles hors ligne
    python3 validate_sequences.py --images   # teste aussi chaque URL (réseau)
"""

import argparse
import sqlite3
import sys

import blast_search as bs
from identify import Identifier, SequenceError, normalize

DEFAULT_DB = "biodetective.db"
LEGO_MAX_RUN = 2      # pas plus de 2 briques identiques d'affilée
SLOW_IMAGE_BYTES = 1_000_000

OK, WARN, FAIL = "  ok  ", " WARN ", " FAIL "


class Report:
    def __init__(self):
        self.errors = 0
        self.warnings = 0

    def line(self, status, msg):
        if status is FAIL:
            self.errors += 1
        elif status is WARN:
            self.warnings += 1
        print(f"[{status}] {msg}")

    def section(self, title):
        print(f"\n--- {title} ---")


def longest_run(seq):
    best = run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def check_images(entries, rep):
    """Teste que chaque URL renvoie bien une image. Nécessite le réseau."""
    import urllib.error
    import urllib.request

    for e in entries:
        url = e["image_path"]
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "BioDetective/1.0"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                ctype = resp.headers.get("Content-Type", "")
                body = resp.read()
        except (urllib.error.URLError, OSError) as exc:
            rep.line(FAIL, f"{e['nom_fr']} : image injoignable ({exc})")
            continue

        if not ctype.startswith("image/"):
            rep.line(FAIL, f"{e['nom_fr']} : la réponse n'est pas une image "
                           f"(Content-Type: {ctype})")
        elif len(body) > SLOW_IMAGE_BYTES:
            rep.line(WARN, f"{e['nom_fr']} : image lourde "
                           f"({len(body) / 1e6:.1f} Mo), lente à afficher")
        else:
            rep.line(OK, f"{e['nom_fr']} : image {len(body) / 1e3:.0f} ko")


def main():
    ap = argparse.ArgumentParser(description="Valide sequences.json")
    ap.add_argument("--db", default=DEFAULT_DB)
    ap.add_argument("--images", action="store_true",
                    help="teste aussi chaque URL d'image (réseau requis)")
    args = ap.parse_args()

    rep = Report()

    # --- 1. Le fichier se charge, aucune séquence en double ---------------
    rep.section("Chargement de sequences.json")
    try:
        ident = Identifier()
    except (OSError, ValueError) as exc:
        print(f"[{FAIL}] {exc}")
        return 1
    rep.line(OK, f"{len(ident)} séquences chargées, aucun doublon")

    # --- 2. Les séquences sont utilisables par un enfant ------------------
    rep.section("Forme des séquences")
    lengths = set()
    for e in ident.entries:
        label = e.get("nom_fr", e["taxonomy_id"])
        try:
            seq = normalize(e["sequence"])
        except SequenceError as exc:
            rep.line(FAIL, f"{label} : {exc}")
            continue
        lengths.add(len(seq))
        run = longest_run(seq)
        if run > LEGO_MAX_RUN:
            rep.line(WARN, f"{label} : {run} briques identiques d'affilée "
                           f"({seq}) — risque d'erreur de comptage")
    if lengths:
        rep.line(OK, f"longueurs présentes : "
                     f"{', '.join(str(n) for n in sorted(lengths))} bases")

    # --- 3. Chaque taxid existe en base et possède une image --------------
    rep.section("Correspondance avec biodetective.db")
    try:
        conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    except sqlite3.OperationalError as exc:
        rep.line(FAIL, f"{args.db} illisible ({exc}). Lancer data_pipeline.py.")
        return 1
    conn.row_factory = sqlite3.Row

    resolved = []
    for e in ident.entries:
        taxid = int(e["taxonomy_id"])
        label = e.get("nom_fr", taxid)
        row = conn.execute(
            "SELECT scientific_name, common_name_fr, organism_type, image_path "
            "FROM organisms WHERE taxonomy_id=?", (taxid,)).fetchone()

        if row is None:
            rep.line(FAIL, f"{label} : taxid {taxid} absent de la base")
            continue
        if not row["image_path"]:
            rep.line(FAIL, f"{label} : taxid {taxid} sans image")
            continue
        if not row["common_name_fr"]:
            rep.line(WARN, f"{label} : pas de nom français en base "
                           f"— compléter common_names_fr.json")
        if e.get("nom_scientifique") and \
                e["nom_scientifique"] != row["scientific_name"]:
            rep.line(WARN, f"{label} : sequences.json dit "
                           f"'{e['nom_scientifique']}', la base dit "
                           f"'{row['scientific_name']}'")
        resolved.append({"nom_fr": label, "image_path": row["image_path"]})

    conn.close()
    if resolved:
        rep.line(OK, f"{len(resolved)}/{len(ident.entries)} organismes "
                     f"résolus avec une image")

    # --- 3 bis. Chaque brin désigne-t-il bien son organisme ? -------------
    rep.section("Vérification BLAST des brins")
    if not bs.blast_available():
        rep.line(FAIL, "blastn introuvable : aucune analyse ne fonctionnera")
    elif not bs.db_available():
        rep.line(FAIL, f"banque BLAST absente ({bs.BLAST_DB})")
    else:
        for e in ident.entries:
            label = e.get("nom_fr", e["taxonomy_id"])
            try:
                hits = bs.search_sync(normalize(e["sequence"]))
            except Exception as exc:                      # noqa: BLE001
                rep.line(FAIL, f"{label} : BLAST a échoué ({exc})")
                continue
            if not hits:
                rep.line(FAIL, f"{label} : aucun hit, ce brin ne donnera rien")
                continue
            top = hits[0]
            if top["taxonomy_id"] != int(e["taxonomy_id"]):
                rep.line(FAIL, f"{label} : BLAST place le taxid "
                               f"{top['taxonomy_id']} en tête, pas "
                               f"{e['taxonomy_id']}")
            elif not bs.is_confident(top):
                rep.line(WARN, f"{label} : hit peu sûr "
                               f"(id {top['percent_identity']} %, "
                               f"cov {top['coverage']} %)")
            else:
                rep.line(OK, f"{label} : id {top['percent_identity']} %, "
                             f"E {top['evalue']:.1g}")

    # --- 4. Les images répondent (optionnel, réseau) ----------------------
    if args.images:
        rep.section("Images distantes")
        check_images(resolved, rep)
    else:
        print("\n(images non testées — relancer avec --images avant la démo)")

    # --- Verdict ----------------------------------------------------------
    print()
    if rep.errors:
        print(f"ÉCHEC : {rep.errors} erreur(s), {rep.warnings} avertissement(s)")
        return 1
    print(f"OK : {len(ident)} séquences valides, "
          f"{rep.warnings} avertissement(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
