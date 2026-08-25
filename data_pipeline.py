#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
data_pipeline.py — Construit biodetective.db à partir de new_taxdump/.

Point d'entrée : images.dmp. On ne garde que les taxons qui ont une image,
puis on va chercher leur nom et leur lignée dans names.dmp et rankedlineage.dmp.

Usage :
    python3 data_pipeline.py                 # construit biodetective.db
    python3 data_pipeline.py --db test.db    # vers un autre fichier
"""

import argparse
import json
import logging
import os
import sqlite3
import sys
import time
from collections import Counter

TAXDUMP_DIR = "new_taxdump"
DEFAULT_DB = "biodetective.db"
NAMES_FR = "common_names_fr.json"
GENOME_STATS = "genome_stats.json"

SEP = "\t|\t"       # séparateur de champs des .dmp NCBI
ROW_END = "\t|"     # fin de ligne des .dmp NCBI

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("biodetective.log", encoding="utf-8"),
    ],
)
log = logging.getLogger("pipeline")


# --------------------------------------------------------------------------
# Lecture des .dmp
# --------------------------------------------------------------------------

def read_dmp(path, expected_cols):
    """Génère les lignes d'un .dmp NCBI sous forme de listes de champs.

    Les fichiers font jusqu'à 400 Mo : on lit en streaming, jamais en mémoire.
    Les colonnes manquantes en fin de ligne sont comblées par des chaînes vides.
    """
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.endswith(ROW_END):
                line = line[: -len(ROW_END)]
            cols = line.split(SEP)
            if len(cols) < expected_cols:
                cols += [""] * (expected_cols - len(cols))
            yield [c.strip() for c in cols[:expected_cols]]


def first_field(line):
    """Extrait le taxid en tête de ligne sans découper toute la ligne.

    Sur rankedlineage.dmp (393 Mo, 2,7 M lignes) ce raccourci évite ~2,7 M
    découpages complets alors qu'on ne garde que 25 000 lignes.
    """
    return line.partition(SEP)[0]


# --------------------------------------------------------------------------
# Normalisation des sources d'images
# --------------------------------------------------------------------------

def normalize_source(raw):
    """images.dmp écrit la même source de plusieurs façons.

    On observe 'iNaturalist', 'iNaturalist.com', 'inaturalist.com',
    'inaturalist.org' pour une seule source, et 'Wikimedia Commons' avec un
    double espace. On normalise ici, pas à l'affichage.
    """
    s = " ".join(raw.split())
    low = s.lower()
    if "inaturalist" in low:
        return "iNaturalist"
    if "wikimedia" in low or "wikipedia" in low:
        return "Wikimedia Commons"
    if not s:
        return "Source inconnue"
    return s


# Priorité quand un taxon a plusieurs images : on préfère Wikimedia (photos
# généralement mieux cadrées et licences plus claires).
SOURCE_PRIORITY = {"Wikimedia Commons": 0, "iNaturalist": 1}


def load_images(path):
    """taxid -> (url, licence, attribution, source). Une seule image par taxon."""
    best = {}
    rows = 0
    for cols in read_dmp(path, 8):
        image_id, _key, url, license_, attribution, source, _props, taxid_list = cols
        if not url or not taxid_list:
            continue
        rows += 1
        source = normalize_source(source)
        try:
            rank = (SOURCE_PRIORITY.get(source, 9), int(image_id))
        except ValueError:
            rank = (SOURCE_PRIORITY.get(source, 9), 1 << 30)
        for taxid in taxid_list.split():
            current = best.get(taxid)
            if current is None or rank < current[0]:
                best[taxid] = (rank, (url, license_, attribution, source))
    log.info("images.dmp : %d lignes -> %d taxons avec image", rows, len(best))
    return {t: v for t, (_r, v) in best.items()}


def load_names(path, wanted):
    """taxid -> (nom scientifique, nom commun anglais) pour les taxids voulus."""
    sci, common = {}, {}
    for line in open(path, encoding="utf-8"):
        taxid = first_field(line)
        if taxid not in wanted:
            continue
        line = line.rstrip("\n")
        if line.endswith(ROW_END):
            line = line[: -len(ROW_END)]
        cols = line.split(SEP)
        if len(cols) < 4:
            continue
        name, name_class = cols[1].strip(), cols[3].strip()
        if name_class == "scientific name":
            sci[taxid] = name
        elif name_class == "genbank common name":
            common[taxid] = name           # le plus fiable
        elif name_class == "common name":
            common.setdefault(taxid, name)  # ne pas écraser un genbank
    log.info("names.dmp  : %d noms scientifiques, %d noms communs EN",
             len(sci), len(common))
    return sci, common


def load_lineages(path, wanted):
    """taxid -> dict des rangs, depuis rankedlineage.dmp.

    Colonnes : taxid | nom | species | genus | family | order | class |
               phylum | kingdom | superkingdom
    Le rang propre du taxon est vide (une espèce a sa colonne 'species' vide,
    son nom est dans la colonne 'nom') : on le rebouche avec le nom.
    """
    out = {}
    for line in open(path, encoding="utf-8"):
        taxid = first_field(line)
        if taxid not in wanted:
            continue
        line = line.rstrip("\n")
        if line.endswith(ROW_END):
            line = line[: -len(ROW_END)]
        cols = [c.strip() for c in line.split(SEP)]
        cols += [""] * (10 - len(cols))
        out[taxid] = {
            "name": cols[1], "species": cols[2], "genus": cols[3],
            "family": cols[4], "order": cols[5], "class": cols[6],
            "phylum": cols[7], "kingdom": cols[8], "superkingdom": cols[9],
        }
    log.info("rankedlineage.dmp : %d lignées", len(out))
    return out


def load_common_names_fr(path):
    """taxid -> nom français, depuis un fichier maintenu à la main.

    Le taxdump NCBI ne contient quasiment aucun nom commun français
    (les 'genbank common name' sont anglais). Les noms français vivent donc
    dans common_names_fr.json, à part, pour survivre à une reconstruction
    de la base. Fichier absent = base construite sans noms français.
    """
    if not os.path.exists(path):
        log.warning("%s absent : aucun nom français ne sera renseigné", path)
        return {}
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    log.info("%s : %d noms français", path, len(data))
    return {str(k): v for k, v in data.items()}


def load_genome_stats(path):
    """taxid -> taille du génome et nombre de chromosomes.

    Produit par build_genome_stats.py depuis les rapports NCBI. Le taxdump
    n'en dit rien : ces chiffres viennent des assemblages.
    """
    if not os.path.exists(path):
        log.warning("%s absent : ni taille de génome ni chromosomes", path)
        return {}
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    log.info("%s : %d génomes renseignés", path, len(data))
    return data


def load_ranks(path, wanted):
    """taxid -> rang ('species', 'genus', ...) depuis nodes.dmp."""
    out = {}
    for line in open(path, encoding="utf-8"):
        taxid = first_field(line)
        if taxid not in wanted:
            continue
        cols = line.split(SEP)
        if len(cols) >= 3:
            out[taxid] = cols[2].strip()
    log.info("nodes.dmp  : %d rangs", len(out))
    return out


# --------------------------------------------------------------------------
# organism_type : étiquette française parlante pour un enfant
# --------------------------------------------------------------------------
#
# NCBI fournit bien un « blast name » mais son vocabulaire compte ~250 valeurs
# (« hooded tickspiders », « micrognathozoans »...) : inutilisable tel quel
# devant un enfant de 8 ans. On dérive donc l'étiquette de la lignée.
#
# Règles évaluées dans l'ordre, première touchée gagne. Ajouter une règle =
# ajouter une ligne ; c'est fait pour être bidouillé la veille de la démo.

TYPE_RULES = [
    # (rang à tester, valeur attendue, étiquette française)

    # --- Vertébrés -------------------------------------------------------
    ("class", "Mammalia",          "Mammifère"),
    ("class", "Aves",              "Oiseau"),
    ("class", "Actinopteri",       "Poisson"),
    ("class", "Actinopterygii",    "Poisson"),
    ("class", "Cladistia",         "Poisson"),
    ("class", "Dipnoi",            "Poisson"),
    ("class", "Coelacanthimorpha", "Poisson"),
    ("class", "Myxini",            "Poisson"),
    ("class", "Hyperoartia",       "Poisson"),
    ("class", "Chondrichthyes",    "Requin & raie"),
    ("class", "Lepidosauria",      "Reptile"),
    ("class", "Testudines",        "Tortue"),
    ("class", "Crocodylia",        "Crocodile"),
    ("class", "Amphibia",          "Amphibien"),
    # Tunicier : Chordata mais pas vertébré, le repli par embranchement
    # les classerait à tort en « Vertébré ».
    ("class", "Ascidiacea",        "Animal marin"),
    ("class", "Thaliacea",         "Animal marin"),

    # --- Invertébrés -----------------------------------------------------
    ("class", "Insecta",           "Insecte"),
    ("class", "Arachnida",         "Araignée & scorpion"),
    ("class", "Malacostraca",      "Crustacé"),
    ("class", "Branchiopoda",      "Crustacé"),
    ("class", "Hexanauplia",       "Crustacé"),
    ("class", "Ostracoda",         "Crustacé"),
    ("class", "Chilopoda",         "Mille-pattes"),
    ("class", "Diplopoda",         "Mille-pattes"),
    ("class", "Gastropoda",        "Escargot & limace"),
    ("class", "Bivalvia",          "Coquillage"),
    ("class", "Cephalopoda",       "Poulpe & calmar"),
    ("class", "Anthozoa",          "Corail & anémone"),
    ("class", "Scyphozoa",         "Méduse"),
    ("class", "Hydrozoa",          "Méduse"),
    ("class", "Cubozoa",           "Méduse"),

    # --- Plantes ---------------------------------------------------------
    ("class", "Magnoliopsida",     "Plante à fleurs"),
    ("class", "Liliopsida",        "Plante à fleurs"),
    ("class", "Pinopsida",         "Conifère"),
    ("class", "Cycadopsida",       "Plante à graines"),
    ("class", "Gnetopsida",        "Plante à graines"),
    ("class", "Polypodiopsida",    "Fougère"),
    ("class", "Lycopodiopsida",    "Lycopode"),
    ("class", "Bryopsida",         "Mousse"),
    ("class", "Sphagnopsida",      "Mousse"),
    ("class", "Polytrichopsida",   "Mousse"),
    ("class", "Andreaeopsida",     "Mousse"),
    ("class", "Marchantiopsida",   "Hépatique"),
    ("class", "Jungermanniopsida", "Hépatique"),

    # --- Algues ----------------------------------------------------------
    ("class", "Ulvophyceae",       "Algue verte"),
    ("class", "Chlorophyceae",     "Algue verte"),
    ("class", "Trebouxiophyceae",  "Algue verte"),
    ("class", "Charophyceae",      "Algue verte"),
    ("class", "Florideophyceae",   "Algue rouge"),
    ("class", "Bangiophyceae",     "Algue rouge"),
    ("class", "Phaeophyceae",      "Algue brune"),
    ("class", "Bacillariophyceae", "Diatomée"),
    ("class", "Dinophyceae",       "Plancton"),

    # --- Champignons -----------------------------------------------------
    # Les Lecanoromycetes sont les lichens : 346 taxons imagés, ça mérite
    # sa propre étiquette plutôt que « Champignon ».
    ("class", "Lecanoromycetes",   "Lichen"),
    ("class", "Arthoniomycetes",   "Lichen"),
    ("class", "Lichinomycetes",    "Lichen"),
    ("class", "Candelariomycetes", "Lichen"),
    ("class", "Agaricomycetes",    "Champignon"),
    ("class", "Saccharomycetes",   "Levure"),

    # --- Règles par ordre ------------------------------------------------
    # Le NCBI n'attribue pas de rang « classe » aux tortues ni aux crocodiles :
    # Testudines et Crocodylia sont des ordres, la colonne classe est vide.
    ("order", "Testudines",        "Tortue"),
    ("order", "Crocodylia",        "Crocodile"),
    ("order", "Ceratodontiformes", "Poisson"),
    ("order", "Coelacanthiformes", "Poisson"),

    # --- Repli par embranchement ----------------------------------------
    ("phylum", "Arthropoda",       "Arthropode"),
    ("phylum", "Mollusca",         "Mollusque"),
    ("phylum", "Cnidaria",         "Cnidaire"),
    ("phylum", "Echinodermata",    "Étoile & oursin"),
    ("phylum", "Porifera",         "Éponge"),
    ("phylum", "Annelida",         "Ver"),
    ("phylum", "Nematoda",         "Ver"),
    ("phylum", "Platyhelminthes",  "Ver"),
    ("phylum", "Chordata",         "Vertébré"),
    ("phylum", "Rhodophyta",       "Algue rouge"),
    ("phylum", "Chlorophyta",      "Algue verte"),
    ("phylum", "Bacillariophyta",  "Diatomée"),
    ("phylum", "Streptophyta",     "Plante"),
    ("phylum", "Basidiomycota",    "Champignon"),
    ("phylum", "Ascomycota",       "Champignon"),
    ("phylum", "Oomycota",         "Moisissure"),
    ("phylum", "Ciliophora",       "Micro-organisme"),
    ("phylum", "Apicomplexa",      "Micro-organisme"),
    ("phylum", "Euglenozoa",       "Micro-organisme"),

    # --- Repli par règne -------------------------------------------------
    ("kingdom", "Metazoa",         "Animal"),
    ("kingdom", "Viridiplantae",   "Plante"),
    ("kingdom", "Fungi",           "Champignon"),

    # --- Repli par domaine ----------------------------------------------
    ("superkingdom", "Bacteria",   "Bactérie"),
    ("superkingdom", "Archaea",    "Archée"),
    ("superkingdom", "Eukaryota",  "Micro-organisme"),
]


# Suffixes imposés par l'ICTV à chaque rang viral. Se fier au seul
# « superkingdom » ne suffit pas : beaucoup de familles virales imagées
# (Baculoviridae, Fuselloviridae...) n'ont aucun realm renseigné.
VIRAL_SUFFIXES = (
    "viria",       # realm      (Riboviria)
    "virae",       # règne      (Orthornavirae)
    "viricota",    # phylum     (Negarnaviricota)
    "viricetes",   # classe     (Naldaviricetes)
    "virales",     # ordre      (Herpesvirales)
    "viridae",     # famille    (Baculoviridae)
    "viriformidae",
    "satellitidae",
)


def is_viral(lin, name):
    if lin.get("superkingdom") == "Viruses":
        return True
    candidates = [name] + [lin.get(r, "") for r in
                           ("superkingdom", "kingdom", "phylum", "class",
                            "order", "family")]
    return any(c.endswith(VIRAL_SUFFIXES) for c in candidates if c)


def classify(lin, name=""):
    """Étiquette française à partir de la lignée. Jamais vide."""
    # Les viroïdes ne sont pas des virus : pas de capside, juste un ARN nu.
    # Ils sont peu nombreux mais autant être exact.
    if name.endswith("viroidae") or lin.get("family", "").endswith("viroidae"):
        return "Viroïde"
    if is_viral(lin, name):
        return "Virus"

    for rank, value, label in TYPE_RULES:
        if lin.get(rank) == value:
            return label
    return "Organisme"


# --------------------------------------------------------------------------
# Écriture SQLite
# --------------------------------------------------------------------------

SCHEMA = """
PRAGMA encoding = 'UTF-8';

DROP TABLE IF EXISTS organisms;
CREATE TABLE organisms (
    taxonomy_id       INTEGER PRIMARY KEY,
    scientific_name   TEXT NOT NULL,
    common_name_fr    TEXT,
    common_name_en    TEXT,
    superkingdom      TEXT,
    kingdom           TEXT,
    phylum            TEXT,
    class_name        TEXT,
    order_name        TEXT,
    family            TEXT,
    genus             TEXT,
    species           TEXT,
    organism_type     TEXT,
    size_info         TEXT,
    habitat           TEXT,
    danger_level      TEXT,
    discovery_info    TEXT,
    fun_facts         TEXT,
    ecological_role   TEXT,
    human_utility     TEXT,
    image_path        TEXT,
    image_source      TEXT,
    image_license     TEXT,
    image_attribution TEXT,
    genome_size_mb    REAL,
    chromosome_count  INTEGER,
    assembly_accession TEXT,
    created_at        TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at        TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_scientific_name ON organisms(scientific_name);
CREATE INDEX idx_organism_type   ON organisms(organism_type);
CREATE INDEX idx_kingdom         ON organisms(kingdom);
CREATE INDEX idx_superkingdom    ON organisms(superkingdom);
"""

INSERT = """
INSERT OR REPLACE INTO organisms (
    taxonomy_id, scientific_name, common_name_fr, common_name_en,
    superkingdom, kingdom, phylum, class_name, order_name, family, genus, species,
    organism_type, image_path, image_source, image_license, image_attribution,
    genome_size_mb, chromosome_count, assembly_accession
) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
"""


def build(db_path, taxdump_dir):
    t0 = time.time()

    images_path = os.path.join(taxdump_dir, "images.dmp")
    names_path = os.path.join(taxdump_dir, "names.dmp")
    ranked_path = os.path.join(taxdump_dir, "rankedlineage.dmp")
    nodes_path = os.path.join(taxdump_dir, "nodes.dmp")

    for p in (images_path, names_path, ranked_path, nodes_path):
        if not os.path.exists(p):
            log.error("Fichier manquant : %s", p)
            log.error("Décompresser new_taxdump.tar.gz dans %s/", taxdump_dir)
            return 1

    images = load_images(images_path)
    wanted = set(images)
    sci, common_en = load_names(names_path, wanted)
    lineages = load_lineages(ranked_path, wanted)
    ranks = load_ranks(nodes_path, wanted)
    names_fr = load_common_names_fr(NAMES_FR)
    genomes = load_genome_stats(GENOME_STATS)

    rows, skipped = [], 0
    type_counts = Counter()

    for taxid in sorted(wanted, key=int):
        name = sci.get(taxid) or lineages.get(taxid, {}).get("name")
        if not name:
            skipped += 1
            continue

        lin = lineages.get(taxid, {})
        # rankedlineage laisse vide la colonne du rang propre du taxon :
        # une espèce a 'species' vide, un genre a 'genus' vide, etc.
        rank = ranks.get(taxid, "")
        if rank in ("species", "genus", "family", "order", "class",
                    "phylum", "kingdom") and not lin.get(rank):
            lin = dict(lin, **{rank: name})

        otype = classify(lin, name)
        type_counts[otype] += 1
        url, license_, attribution, source = images[taxid]

        rows.append((
            int(taxid), name, names_fr.get(taxid), common_en.get(taxid),
            lin.get("superkingdom") or None,
            lin.get("kingdom") or None, lin.get("phylum") or None,
            lin.get("class") or None, lin.get("order") or None,
            lin.get("family") or None, lin.get("genus") or None,
            lin.get("species") or None,
            otype, url, source, license_ or None, attribution or None,
            (genomes.get(taxid) or {}).get("genome_size_mb"),
            (genomes.get(taxid) or {}).get("chromosomes") or None,
            (genomes.get(taxid) or {}).get("assembly"),
        ))

    conn = sqlite3.connect(db_path)
    try:
        conn.executescript(SCHEMA)
        conn.executemany(INSERT, rows)
        conn.commit()
    finally:
        conn.close()

    log.info("--------------------------------------------------")
    log.info("%d organismes écrits dans %s (%.1f s)",
             len(rows), db_path, time.time() - t0)
    if skipped:
        log.warning("%d taxons ignorés (aucun nom trouvé)", skipped)
    log.info("Types les plus fréquents :")
    for label, n in type_counts.most_common(15):
        log.info("    %-22s %6d", label, n)
    with_genome = sum(1 for r in rows if r[17])
    log.info("Génomes renseignés : %d / %d", with_genome, len(rows))
    matched_fr = sum(1 for r in rows if r[2])
    log.info("Noms français renseignés : %d / %d", matched_fr, len(rows))
    unknown = type_counts.get("Organisme", 0)
    if unknown:
        log.warning("%d organismes non classés -> compléter TYPE_RULES", unknown)
    return 0


def main():
    ap = argparse.ArgumentParser(description="Construit biodetective.db")
    ap.add_argument("--db", default=DEFAULT_DB, help="fichier SQLite de sortie")
    ap.add_argument("--taxdump", default=TAXDUMP_DIR, help="dossier new_taxdump")
    args = ap.parse_args()
    return build(args.db, args.taxdump)


if __name__ == "__main__":
    sys.exit(main())
