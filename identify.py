#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
identify.py — Identification d'un organisme à partir d'une séquence ADN.

Correspondance EXACTE contre sequences.json, après normalisation.
Pas de BLAST : voir CLAUDE.md pour le pourquoi.

Usage en ligne de commande (pratique pour tester sans lancer l'API) :
    python3 identify.py AACTGCAGCAAGGTAT
"""

import json
import os
import sys

SEQUENCES_FILE = "sequences.json"
VALID_BASES = frozenset("ATCG")

# Longueur maximale acceptée. Les brins du Brickopore font 16 bases ; on
# accepte large pour ne pas rejeter un enfant bavard, mais on borne pour
# éviter qu'un copier-coller d'un génome entier ne parte dans la base.
MAX_LENGTH = 1000


class SequenceError(ValueError):
    """Séquence inutilisable. Le message est destiné à être affiché tel quel."""


def normalize(raw):
    """Nettoie une séquence saisie à la main puis la valide.

    Un enfant colle rarement une chaîne propre : espaces, retours à la ligne,
    minuscules, parfois une numérotation de type FASTA. On absorbe tout ça,
    puis on refuse ce qui n'est franchement pas de l'ADN.

    Lève SequenceError avec un message affichable en l'état.
    """
    if raw is None:
        raise SequenceError("Aucune séquence n'a été saisie.")

    # Une ligne d'en-tête FASTA (>...) est ignorée plutôt que rejetée :
    # c'est un copier-coller courant et l'intention est claire.
    lines = [l for l in raw.splitlines() if not l.lstrip().startswith(">")]
    seq = "".join(lines)
    seq = "".join(seq.split()).upper()

    # Le U de l'ARN est toléré et converti : l'erreur est compréhensible
    # et la corriger vaut mieux que renvoyer l'enfant à son clavier.
    seq = seq.replace("U", "T")

    if not seq:
        raise SequenceError("Aucune séquence n'a été saisie.")
    if len(seq) > MAX_LENGTH:
        raise SequenceError(
            f"Cette séquence est bien trop longue ({len(seq)} bases). "
            f"Le Brickopore en produit 16."
        )

    bad = sorted(set(seq) - VALID_BASES)
    if bad:
        raise SequenceError(
            "Cette séquence contient des lettres qui ne sont pas de l'ADN : "
            + ", ".join(bad)
            + ". Seuls A, T, C et G sont acceptés."
        )
    return seq


class Identifier:
    """Table de correspondance séquence -> taxonomy_id, chargée une fois."""

    def __init__(self, path=SEQUENCES_FILE):
        self.path = path
        self.table = {}
        self.entries = []
        self.load()

    def load(self):
        if not os.path.exists(self.path):
            raise FileNotFoundError(
                f"{self.path} introuvable : aucune séquence ne pourra être "
                f"identifiée."
            )
        with open(self.path, encoding="utf-8") as fh:
            doc = json.load(fh)

        entries = doc.get("sequences", [])
        table, seen = {}, {}
        for entry in entries:
            seq = normalize(entry["sequence"])
            taxid = int(entry["taxonomy_id"])
            if seq in seen and seen[seq] != taxid:
                raise ValueError(
                    f"La séquence {seq} est associée à deux organismes "
                    f"différents ({seen[seq]} et {taxid}) dans {self.path}."
                )
            seen[seq] = taxid
            table[seq] = taxid
        self.table = table
        self.entries = entries

    def identify(self, raw):
        """Retourne le taxonomy_id, ou None si la séquence est inconnue.

        None n'est PAS une erreur : les enfants assemblent les briques
        librement, donc la plupart des séquences sont inconnues. C'est le
        cas nominal le plus fréquent.
        """
        return self.table.get(normalize(raw))

    def __len__(self):
        return len(self.table)


def main():
    if len(sys.argv) != 2:
        print(__doc__.strip())
        return 2
    ident = Identifier()
    try:
        seq = normalize(sys.argv[1])
    except SequenceError as exc:
        print(f"Séquence invalide : {exc}")
        return 1

    taxid = ident.identify(seq)
    print(f"Séquence normalisée : {seq} ({len(seq)} bases)")
    if taxid is None:
        print(f"Aucune correspondance parmi {len(ident)} séquences connues.")
        return 0

    match = next(e for e in ident.entries
                 if normalize(e["sequence"]) == seq)
    print(f"-> taxid {taxid} : {match['nom_fr']} "
          f"({match['nom_scientifique']}, {match['type']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
