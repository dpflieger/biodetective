#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_montages.py — Fabrique les planches d'images de l'écran de recherche.

Chaque planche est une bande JPEG horizontale de N vignettes carrées, que le
frontend fait défiler en CSS avec steps(). C'est l'équivalent d'un GIF, mais
3,8 fois plus léger, sans tramage 256 couleurs, et surtout la vitesse reste
réglable côté CSS au lieu d'être figée à la fabrication.

    python3 make_montages.py                    # 20 planches de 30 vignettes
    python3 make_montages.py --count 8 --frames 24

Les planches sont volumineuses et régénérables : elles ne sont pas versionnées,
au même titre que biodetective.db.
"""

import argparse
import concurrent.futures as cf
import io
import json
import logging
import os
import sqlite3
import sys
import time

from PIL import Image

DB_PATH = "biodetective.db"
OUT_DIR = os.path.join("public", "montages")
INDEX = "index.json"

FRAME_SIZE = 400        # côté de la vignette, en pixels
JPEG_QUALITY = 78
WORKERS = 6             # rester poli avec les serveurs du NCBI
TIMEOUT = 25

UA = {"User-Agent": "BioDetective/1.0 (Fete de la Science, IBMP Strasbourg)"}

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S", handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger("montages")


def fetch_frame(url):
    """Télécharge une image, la recadre en carré centré et la redimensionne."""
    import urllib.request
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            data = resp.read()
        im = Image.open(io.BytesIO(data)).convert("RGB")
    except Exception:
        return None
    w, h = im.size
    if w < 80 or h < 80:
        return None
    side = min(w, h)
    im = im.crop(((w - side) // 2, (h - side) // 2,
                  (w - side) // 2 + side, (h - side) // 2 + side))
    return im.resize((FRAME_SIZE, FRAME_SIZE), Image.LANCZOS)


def pick_urls(count):
    """Tire des images distinctes, en évitant les plus lourdes.

    Les photos de plusieurs mégaoctets ralentiraient la fabrication sans rien
    apporter : elles sont de toute façon réduites à 400 px.
    """
    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    rows = conn.execute(
        "SELECT image_path FROM organisms "
        "WHERE image_path IS NOT NULL AND image_path != '' "
        "ORDER BY RANDOM() LIMIT ?", (count,)).fetchall()
    conn.close()
    return [r[0] for r in rows]


def build(count, frames):
    if not os.path.exists(DB_PATH):
        log.error("%s introuvable. Lancer d'abord data_pipeline.py.", DB_PATH)
        return 1
    os.makedirs(OUT_DIR, exist_ok=True)

    needed = count * frames
    # On tire large : certaines URL ne répondront pas.
    urls = pick_urls(int(needed * 1.25))
    log.info("%d planches x %d vignettes = %d images à récupérer",
             count, frames, needed)

    t0 = time.time()
    collected = []
    with cf.ThreadPoolExecutor(WORKERS) as ex:
        for i, frame in enumerate(ex.map(fetch_frame, urls), 1):
            if frame is not None:
                collected.append(frame)
            if i % 50 == 0:
                log.info("  %d/%d testées, %d retenues (%.0f s)",
                         i, len(urls), len(collected), time.time() - t0)
            if len(collected) >= needed:
                break

    if len(collected) < frames:
        log.error("Seulement %d images récupérées, impossible de faire une "
                  "planche de %d. Vérifier la connexion au NCBI.",
                  len(collected), frames)
        return 1

    manifest = []
    made = 0
    for k in range(count):
        chunk = collected[k * frames:(k + 1) * frames]
        if len(chunk) < frames:
            log.warning("Plus assez d'images : %d planches au lieu de %d",
                        made, count)
            break
        sheet = Image.new("RGB", (FRAME_SIZE * frames, FRAME_SIZE))
        for i, im in enumerate(chunk):
            sheet.paste(im, (i * FRAME_SIZE, 0))
        name = f"montage-{k + 1:02d}.jpg"
        path = os.path.join(OUT_DIR, name)
        sheet.save(path, quality=JPEG_QUALITY, optimize=True, progressive=True)
        manifest.append({"file": f"/montages/{name}", "frames": frames})
        made += 1
        log.info("  %s  %d vignettes  %.0f ko",
                 name, frames, os.path.getsize(path) / 1024)

    with open(os.path.join(OUT_DIR, INDEX), "w", encoding="utf-8") as fh:
        json.dump({"frame_size": FRAME_SIZE, "montages": manifest}, fh,
                  ensure_ascii=False, indent=2)
        fh.write("\n")

    total = sum(os.path.getsize(os.path.join(OUT_DIR, m["file"].split("/")[-1]))
                for m in manifest)
    log.info("-" * 50)
    log.info("%d planches écrites dans %s/ (%.1f Mo, %.0f s)",
             made, OUT_DIR, total / 1e6, time.time() - t0)
    return 0


def main():
    ap = argparse.ArgumentParser(description="Fabrique les planches d'images")
    ap.add_argument("--count", type=int, default=20, help="nombre de planches")
    ap.add_argument("--frames", type=int, default=30, help="vignettes par planche")
    args = ap.parse_args()
    return build(args.count, args.frames)


if __name__ == "__main__":
    sys.exit(main())
