#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
remote_blast.py — Interroger le BLAST du NCBI, à distance.

Utile pour vérifier ce que dit le NCBI sur une séquence, sans banque locale.
**Pas pour la démonstration** : une recherche distante prend des dizaines de
secondes à plusieurs minutes, là où l'écran de recherche dure 2,5 à 13 s.

    python3 remote_blast.py TCATTGTAAGATTGGAATAATTCAATTTCGACAT
    python3 remote_blast.py --database core_nt --taxid 3702 ATCG…

Le NCBI demande de s'identifier et de ne pas marteler ses serveurs :
    export BIODETECTIVE_NCBI_EMAIL="prenom.nom@exemple.fr"
L'adresse n'est envoyée que si cette variable est définie.
"""

import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "BioDetective"

# Le NCBI demande de ne pas réinterroger un même RID plus d'une fois par
# minute. On reste au-dessus de sa borne de politesse la plus stricte pour
# les soumissions, et raisonnable pour l'attente.
POLL_SECONDS = 20
MAX_WAIT_SECONDS = 600


def _email_params():
    email = os.environ.get("BIODETECTIVE_NCBI_EMAIL", "").strip()
    p = {"tool": TOOL}
    if email:
        p["email"] = email
    return p


def _get(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": f"{TOOL}/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def submit(sequence, database="core_nt", taxid=None, word_size=7,
           expect=1000.0, hitlist=50):
    """Dépose la requête et retourne (RID, durée estimée en secondes)."""
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "MEGABLAST": "off",          # mégablast ignore les requêtes courtes
        "DATABASE": database,
        "QUERY": sequence,
        "WORD_SIZE": str(word_size), # 7 : indispensable sous ~30 pb
        "EXPECT": str(expect),
        "HITLIST_SIZE": str(hitlist),
        "FILTER": "F",               # pas de masquage : nos brins sont courts
    }
    if taxid:
        params["TAXIDS"] = str(taxid)
    params.update(_email_params())

    body = _get(URL + "?" + urllib.parse.urlencode(params))
    rid = rtoe = None
    for line in body.splitlines():
        line = line.strip()
        if line.startswith("RID = "):
            rid = line[6:].strip()
        elif line.startswith("RTOE = "):
            try:
                rtoe = int(line[7:].strip())
            except ValueError:
                pass
    if not rid:
        snippet = " ".join(body.split())[:300]
        raise RuntimeError(f"Le NCBI n'a pas rendu de RID : {snippet}")
    return rid, rtoe or 10


def status(rid):
    """WAITING, READY, FAILED ou UNKNOWN."""
    params = {"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"}
    params.update(_email_params())
    body = _get(URL + "?" + urllib.parse.urlencode(params))
    st = "UNKNOWN"
    hits = None
    for line in body.splitlines():
        line = line.strip()
        if line.startswith("Status="):
            st = line[7:].strip()
        elif line.startswith("ThereAreHits="):
            hits = line[13:].strip() == "yes"
    return st, hits


def fetch_xml(rid):
    """Récupère le résultat en XML.

    FORMAT_TYPE=Tabular rend un tableau vide sur cette interface, même quand
    la recherche a trouvé : le XML est le seul format tabulable fiable, et il
    porte en prime les séquences alignées et la longueur du sujet.
    """
    params = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}
    params.update(_email_params())
    return _get(URL + "?" + urllib.parse.urlencode(params), timeout=300)


def parse_xml(text):
    """Transforme le XML BLAST en la même forme que les hits locaux."""
    import xml.etree.ElementTree as ET
    start = text.find("<?xml")
    if start > 0:
        text = text[start:]
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise RuntimeError(f"XML illisible du NCBI : {exc}")

    hits = []
    for hit in root.iter("Hit"):
        acc = (hit.findtext("Hit_accession") or "").strip()
        hit_def = (hit.findtext("Hit_def") or "").strip()
        try:
            slen = int(hit.findtext("Hit_len") or 0)
        except ValueError:
            slen = 0
        for hsp in hit.iter("Hsp"):
            def num(tag, cast=int, default=0):
                v = hsp.findtext(tag)
                try:
                    return cast(v)
                except (TypeError, ValueError):
                    return default
            align_len = num("Hsp_align-len") or 1
            ident = num("Hsp_identity")
            qs, qe = num("Hsp_query-from"), num("Hsp_query-to")
            hs, he = num("Hsp_hit-from"), num("Hsp_hit-to")
            qseq = hsp.findtext("Hsp_qseq") or ""
            sseq = hsp.findtext("Hsp_hseq") or ""
            hits.append({
                "accession": acc,
                "subject_title": hit_def,
                "subject_length": slen or None,
                "percent_identity": round(100.0 * ident / align_len, 1),
                "alignment_length": align_len,
                "mismatches": align_len - ident,
                "gaps": num("Hsp_gaps"),
                "query_start": qs, "query_end": qe,
                "subject_start": hs, "subject_end": he,
                "evalue": num("Hsp_evalue", float, 0.0),
                "bitscore": round(num("Hsp_bit-score", float, 0.0), 1),
                "query_seq": qseq,
                "subject_seq": sseq,
                "strand": "moins" if hs > he else "plus",
            })
    hits.sort(key=lambda h: (h["evalue"], -h["bitscore"]))
    return hits


def resolve_taxids(accessions):
    """accession -> (taxid, titre). Le Tabular du NCBI ne donne pas le taxid."""
    if not accessions:
        return {}
    params = {"db": "nuccore", "id": ",".join(accessions[:50]),
              "retmode": "json"}
    params.update(_email_params())
    try:
        data = json.loads(_get(EUTILS + "/esummary.fcgi?"
                               + urllib.parse.urlencode(params)))
    except Exception:
        return {}
    out = {}
    res = data.get("result", {})
    for uid in res.get("uids", []):
        rec = res.get(uid, {})
        val = (rec.get("taxid"), rec.get("title", ""))
        # Le XML de BLAST rend « AY052207 », esummary « AY052207.1 » :
        # on indexe sous les deux formes pour que la jonction se fasse.
        for key in (rec.get("accessionversion"), rec.get("caption")):
            if key:
                out[key] = val
                out[key.rsplit(".", 1)[0]] = val
    return out


def search(sequence, database="core_nt", taxid=None, verbose=True,
           max_wait=MAX_WAIT_SECONDS):
    """Chaîne complète. Retourne (hits, RID, secondes écoulées)."""
    t0 = time.time()
    rid, rtoe = submit(sequence, database=database, taxid=taxid)
    if verbose:
        print(f"  RID {rid} — le NCBI annonce ~{rtoe} s", flush=True)

    time.sleep(min(rtoe, 30))
    while True:
        st, has_hits = status(rid)
        waited = time.time() - t0
        if verbose:
            print(f"  {waited:6.0f} s  {st}", flush=True)
        if st == "READY":
            break
        if st == "FAILED":
            raise RuntimeError(f"Le NCBI a échoué sur le RID {rid}")
        if st == "UNKNOWN":
            raise RuntimeError(f"RID {rid} inconnu ou expiré")
        if waited > max_wait:
            raise TimeoutError(
                f"Toujours rien après {waited:.0f} s (RID {rid}, "
                f"consultable un moment sur le site du NCBI)")
        time.sleep(POLL_SECONDS)

    if has_hits is False:
        return [], rid, time.time() - t0
    hits = parse_xml(fetch_xml(rid))
    for h in hits:
        h["coverage"] = round(100.0 * h["alignment_length"] / len(sequence), 1)
    meta = resolve_taxids([h["accession"] for h in hits[:20]])
    for h in hits:
        tx, title = meta.get(h["accession"], (None, None))
        h["taxonomy_id"] = tx
        if title:
            h["subject_title"] = title
    return hits, rid, time.time() - t0


def main():
    ap = argparse.ArgumentParser(description="BLAST distant chez le NCBI")
    ap.add_argument("sequence")
    ap.add_argument("--database", default="core_nt",
                    help="core_nt, nt, refseq_rna, ref_euk_rep_genomes…")
    ap.add_argument("--taxid", help="restreindre à un taxon")
    ap.add_argument("--max-wait", type=int, default=MAX_WAIT_SECONDS)
    args = ap.parse_args()

    from identify import SequenceError, normalize
    try:
        seq = normalize(args.sequence)
    except SequenceError as exc:
        print(f"Séquence invalide : {exc}")
        return 1

    if not os.environ.get("BIODETECTIVE_NCBI_EMAIL"):
        print("(BIODETECTIVE_NCBI_EMAIL non définie : le NCBI demande une "
              "adresse de contact)")
    print(f"Requête : {seq} ({len(seq)} pb) contre {args.database}")
    try:
        hits, rid, el = search(seq, database=args.database, taxid=args.taxid,
                               max_wait=args.max_wait)
    except (RuntimeError, TimeoutError) as exc:
        print(f"Échec : {exc}")
        return 1

    print(f"\n{len(hits)} hits en {el:.0f} s (RID {rid})\n")
    for h in hits[:8]:
        print(f"  {h['accession']:<16} taxid {str(h.get('taxonomy_id')):<9} "
              f"id {h['percent_identity']:5.1f} %  E {h['evalue']:<9.2g} "
              f"score {h['bitscore']:<6} {h.get('subject_title','')[:44]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
