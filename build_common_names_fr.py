#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_common_names_fr.py — Récolte les noms français sur Wikidata.

Le taxdump NCBI ne contient **aucun** nom commun français : ses
« genbank common name » sont anglais. Wikidata, elle, relie ses fiches aux
taxons NCBI par la propriété P685, et porte des noms vernaculaires par langue.

    python3 build_common_names_fr.py               # -> common_names_fr_wikidata.json
    python3 build_common_names_fr.py --limit 500   # essai rapide
    python3 build_common_names_fr.py --update-db   # écrit aussi dans la base

Puis, pour que la base les reprenne :

    python3 data_pipeline.py

Les noms écrits à la main dans `common_names_fr.json` **gagnent toujours** :
ce fichier-ci ne les touche pas, et data_pipeline.py les applique en dernier.

Quelle propriété, et pourquoi
-----------------------------
Trois sources possibles, mesurées sur les mêmes taxons :

    taxon        rdfs:label fr        P1843 fr
    Apis         « Apis mellifera »   « Abeille européenne »   <- P1843 gagne
    Arabidopsis  « Arabidopsis … »    « Arabette des dames »   <- P1843 gagne
    Zea mays     « maïs »             « Maïs »                 <- égalité
    E. coli      « Escherichia coli » (aucune)                 <- rien à prendre

`rdfs:label` est souvent le nom scientifique recopié : le prendre tel quel
remplirait `common_name_fr` de latin, ce qui est **pire que vide** — le nom
français passe avant l'anglais dans `display_name`, donc un faux nom français
masquerait un vrai nom anglais. On prend donc P1843 en premier, le label
seulement s'il diffère du nom scientifique.

Rendement mesuré sur 500 espèces tirées au sort : 487 connues de Wikidata,
84 avec un P1843 français, 35 de plus avec un label exploitable, soit
**24 %** — environ 6 100 noms sur les 25 545 organismes de la base.
"""

import argparse
import json
import os
import re
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

ENDPOINT = "https://query.wikidata.org/sparql"

# Wikidata exige un User-Agent identifiable et refuse les requêtes anonymes.
USER_AGENT = ("BioDetective/1.0 (https://github.com/dpflieger/Biodetective; "
              "IBMP Strasbourg; vulgarisation scientifique)")

DEFAULT_DB = "biodetective.db"
DEFAULT_OUT = "common_names_fr_wikidata.json"
MANUAL_FILE = "common_names_fr.json"

# 250 taxids par requête : au-delà, WDQS commence à rendre des 500 sur une
# clause VALUES trop grosse. Mesuré : 250 taxids en ~4 s.
BATCH = 250
PAUSE = 1.0          # entre deux requêtes, par correction
TIMEOUT = 120
RETRIES = 4

# P225 est le nom scientifique porté par la fiche Wikidata elle-même. On le
# demande pour pouvoir rejeter les candidats qui n'en sont qu'une recopie —
# y compris sous un SYNONYME, que la comparaison avec notre propre nom
# scientifique ne verrait pas passer (« Ephemerocybe congregata » chez nous,
# « Coprinellus congregatus » chez eux).
QUERY = """SELECT ?taxid ?vern ?label ?sciname WHERE {
  VALUES ?taxid { %s }
  ?item wdt:P685 ?taxid .
  OPTIONAL { ?item wdt:P1843 ?vern . FILTER(LANG(?vern) = "fr") }
  OPTIONAL { ?item rdfs:label ?label . FILTER(LANG(?label) = "fr") }
  OPTIONAL { ?item wdt:P225 ?sciname }
}"""

# Un binôme latin : « Apis mellifera », « Zea mays subsp. mays ».
BINOMIAL = re.compile(r"^[A-Z][a-z]+ [a-z][a-z-]+")


# --------------------------------------------------------------------------
# Tri des candidats
# --------------------------------------------------------------------------

def looks_scientific(name, scientific, genus, wd_sciname=None):
    """Ce « nom français » n'est-il que du latin déguisé ?

    C'est le filtre qui compte. Un nom scientifique recopié dans
    common_name_fr passerait devant un vrai nom anglais à l'affichage :
    mieux vaut laisser la case vide.

    Quatre pièges, tous rencontrés sur un essai de 1 000 taxons :
      * la recopie pure et simple du nom scientifique ;
      * un binôme commençant par notre genre ;
      * un **synonyme** sous un autre genre — d'où la comparaison avec le
        P225 de la fiche Wikidata, et pas seulement avec notre nom à nous ;
      * l'épithète seule (« Ni » pour *Trichoplusia ni*).
    """
    low = name.strip().lower()
    if len(low) < 3:
        return True
    if low == (scientific or "").lower():
        return True
    if wd_sciname and low == wd_sciname.strip().lower():
        return True
    parts = (scientific or "").split()
    if len(parts) >= 2 and low == parts[1].lower():
        return True                       # l'épithète prise pour un nom
    if genus and low == genus.lower():
        return True
    if genus and BINOMIAL.match(name) and low.startswith(genus.lower() + " "):
        return True
    return False


def best_vernacular(values):
    """Le meilleur des P1843. Wikidata en donne souvent plusieurs.

    Deux règles, dans l'ordre :
      * majuscule initiale d'abord — Wikidata double fréquemment une valeur
        en bas de casse (« Cigogne blanche » et « cigogne blanche ») ;
      * puis la plus courte, qui est le nom de base plutôt qu'une précision
        régionale (« Lion » plutôt que « Lion d'Afrique »).
    """
    if not values:
        return None
    upper = [v for v in values if v[:1].isupper()]
    return sorted(upper or values, key=lambda v: (len(v), v))[0]


def capitalize(name):
    """« lion » -> « Lion ». Les 33 noms écrits à la main le sont ainsi."""
    name = name.strip()
    return name[:1].upper() + name[1:] if name else name


def choose(entry, scientific, genus):
    """Le nom retenu pour un taxon, ou None s'il n'y a rien d'exploitable."""
    wd = entry.get("sciname")
    ok = lambda n: n and not looks_scientific(n, scientific, genus, wd)

    # On écarte d'abord les candidats latins, PUIS on choisit : sinon un
    # synonyme latin plus court éliminerait le vrai nom français.
    vern = best_vernacular([v for v in entry.get("vern", ()) if ok(v)])
    if vern:
        return capitalize(vern)
    if ok(entry.get("label")):
        return capitalize(entry["label"])
    return None


# --------------------------------------------------------------------------
# Interrogation
# --------------------------------------------------------------------------

def ask(taxids):
    """Un lot de taxids -> {taxid: {"vern": set, "label": str}}."""
    query = QUERY % " ".join('"%s"' % t for t in taxids)
    url = ENDPOINT + "?" + urllib.parse.urlencode(
        {"query": query, "format": "json"})
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as fh:
                data = json.load(fh)
            break
        except (urllib.error.URLError, OSError, ValueError) as exc:
            wait = 5 * (attempt + 1)
            print(f"    ! {exc} — nouvel essai dans {wait} s", file=sys.stderr)
            time.sleep(wait)
    else:
        raise RuntimeError("Wikidata ne répond pas après "
                           f"{RETRIES} tentatives.")

    out = {}
    for b in data["results"]["bindings"]:
        taxid = b["taxid"]["value"]
        e = out.setdefault(taxid, {"vern": set(), "label": None,
                                   "sciname": None})
        if "vern" in b:
            e["vern"].add(b["vern"]["value"])
        if "label" in b:
            e["label"] = b["label"]["value"]
        if "sciname" in b:
            e["sciname"] = b["sciname"]["value"]
    return out


# --------------------------------------------------------------------------
# Fichiers
# --------------------------------------------------------------------------

def load_previous(path):
    """Reprend une récolte interrompue plutôt que de tout redemander."""
    if not os.path.exists(path):
        return {}, set()
    with open(path, encoding="utf-8") as fh:
        doc = json.load(fh)
    return doc.get("names", {}), set(str(t) for t in doc.get("queried", []))


def save(path, names, queried):
    """Écriture atomique : une récolte de cinq minutes ne doit pas être
    tronquée par un Ctrl-C au mauvais moment."""
    doc = {
        "_comment": [
            "Noms français récoltés sur Wikidata (P1843, sinon rdfs:label fr).",
            "GÉNÉRÉ — ne pas éditer à la main : les corrections vont dans",
            "common_names_fr.json, qui gagne toujours sur ce fichier.",
            "Régénérer : python3 build_common_names_fr.py",
        ],
        "_source": "query.wikidata.org — P685 (NCBI taxid) -> P1843 / label fr",
        "_date": time.strftime("%Y-%m-%d"),
        "_queried_count": len(queried),
        "names": dict(sorted(names.items(), key=lambda kv: int(kv[0]))),
        "queried": sorted(int(t) for t in queried),
    }
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, ensure_ascii=False, indent=1)
    os.replace(tmp, path)


def update_db(db_path, names):
    """Écrit les noms directement dans la base, sans tout reconstruire.

    Raccourci commode quand new_taxdump/ n'est plus sur le disque.
    `data_pipeline.py` reste la voie normale, et la seule qui garantisse que
    la base entière est cohérente avec ses fichiers sources.
    """
    conn = sqlite3.connect(db_path)
    manual = {}
    if os.path.exists(MANUAL_FILE):
        with open(MANUAL_FILE, encoding="utf-8") as fh:
            manual = {str(k): v for k, v in json.load(fh).items()}
    n = 0
    for taxid, name in names.items():
        if taxid in manual:      # le fichier à la main a le dernier mot
            continue
        cur = conn.execute(
            "UPDATE organisms SET common_name_fr = ?, "
            "updated_at = CURRENT_TIMESTAMP WHERE taxonomy_id = ?",
            (name, int(taxid)))
        n += cur.rowcount
    conn.commit()
    conn.close()
    return n


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--db", default=DEFAULT_DB)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--batch", type=int, default=BATCH)
    ap.add_argument("--limit", type=int, default=0,
                    help="n'interroger que N taxons (essai)")
    ap.add_argument("--restart", action="store_true",
                    help="tout redemander, sans reprendre la récolte")
    ap.add_argument("--update-db", action="store_true",
                    help="écrire aussi dans la base, sans data_pipeline.py")
    args = ap.parse_args()

    if not os.path.exists(args.db):
        print(f"{args.db} introuvable. Lancer d'abord data_pipeline.py.")
        return 1

    conn = sqlite3.connect(args.db)
    conn.row_factory = sqlite3.Row
    rows = conn.execute("SELECT taxonomy_id, scientific_name, genus "
                        "FROM organisms").fetchall()
    conn.close()
    info = {str(r["taxonomy_id"]): (r["scientific_name"], r["genus"])
            for r in rows}

    names, queried = ({}, set()) if args.restart else load_previous(args.out)
    todo = [t for t in info if t not in queried]
    if args.limit:
        todo = todo[:args.limit]

    print(f"{len(info)} organismes en base")
    if queried:
        print(f"  {len(queried)} déjà interrogés, {len(names)} noms acquis")
    print(f"  {len(todo)} à interroger, par lots de {args.batch}")
    if not todo:
        print("Rien à faire.")
        return 0

    t0 = time.time()
    batches = (len(todo) + args.batch - 1) // args.batch
    try:
        for i in range(0, len(todo), args.batch):
            chunk = todo[i:i + args.batch]
            found = ask(chunk)
            new = 0
            for taxid in chunk:
                queried.add(taxid)
                entry = found.get(taxid)
                if not entry:
                    continue
                sci, genus = info[taxid]
                name = choose(entry, sci, genus)
                if name:
                    names[taxid] = name
                    new += 1
            done = i // args.batch + 1
            print(f"  lot {done}/{batches} — {new} noms — "
                  f"{len(names)} au total — {time.time() - t0:.0f} s")
            if done % 10 == 0:
                save(args.out, names, queried)   # filet de sécurité
            time.sleep(PAUSE)
    except KeyboardInterrupt:
        print("\nInterrompu — la récolte est sauvegardée, relancer pour "
              "reprendre où elle en est.")
    except RuntimeError as exc:
        print(f"\n{exc}\nLa récolte partielle est sauvegardée.")

    save(args.out, names, queried)
    print(f"\n{len(names)} noms français dans {args.out} "
          f"({100 * len(names) / max(1, len(queried)):.0f} % des "
          f"{len(queried)} taxons interrogés)")

    if args.update_db:
        n = update_db(args.db, names)
        print(f"{n} lignes mises à jour dans {args.db}")
    else:
        print("Pour que la base les reprenne : python3 data_pipeline.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
