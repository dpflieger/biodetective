# CLAUDE.md - Guide de développement BioDetective

## ⚠️ État du projet : APPLICATION COMPLÈTE ET FONCTIONNELLE

Ce fichier a d'abord été un cahier des charges. Tout y est maintenant réel et testé
dans un navigateur.

**Fait :**
- `data_pipeline.py` → `biodetective.db` (25 545 organismes, ~6 s de construction)
- `common_names_fr.json` → 33 noms français
- `sequences.json` → table de correspondance de la démo (33 organismes)
- `identify.py` → normalisation + recherche exacte
- `validate_sequences.py` → contrôle avant démo
- `api.py` → API FastAPI complète, 8 routes, testée de bout en bout
- Frontend React complet : accueil, recherche, résultat, séquence inconnue, erreur
- Polices Orbitron/Rajdhani servies en local (57 ko, aucune dépendance réseau)

**Reste à faire :** les fiches descriptives, le cache d'images local, et la
répétition générale en salle. Voir la roadmap.

Voir la [Roadmap](#-roadmap-ordre-de-construction) en fin de fichier pour l'ordre de construction.

---

## 🧬 Présentation du projet

**BioDetective** est une application de vulgarisation scientifique développée pour la **Fête de la Science** à l'IBMP (Institut de Biologie Moléculaire des Plantes, Strasbourg). Elle permet à des enfants et parents de découvrir la bioinformatique en identifiant des organismes à partir de séquences ADN construites avec des LEGOs colorés (dispositif Brickopore).

### Concept
Un séquenceur LEGO (Brickopore) génère des séquences ADN à partir de briques colorées. L'utilisateur colle cette séquence dans BioDetective, qui identifie l'organisme et affiche sa photo avec une fiche descriptive dans une interface futuriste façon série policière.

### ⚠️ Pas de BLAST réel — décision assumée

L'identification se fait par **table de correspondance exacte** (`sequences.json`), pas par BLAST.

**Pourquoi :** les séquences du Brickopore font 20 à 100 bases. Contre une base comme `nt`
(~10¹² lettres), un 20-mer parfait a une E-value attendue d'environ **35** — autrement dit
`blastn-short` ne trouverait rien d'exploitable en dessous de ~45-50 bases, et les courts motifs
conservés sont partagés par des milliers de taxons. Le « meilleur hit » serait souvent un
organisme arbitraire et sans photo. La table de correspondance donne un résultat **déterministe** :
on choisit à l'avance des espèces avec de belles photos et un fort intérêt pour les enfants.

Conséquences à ne jamais oublier en codant :
- **BLAST+ n'est pas une dépendance.** Ne pas réintroduire `blastn`, `makeblastdb` ou la base `nt`.
- Le pourcentage d'identité affiché vaut **toujours 100 %** (correspondance exacte).
- L'API garde une **exécution asynchrone avec polling** alors qu'une recherche en table est
  instantanée. C'est **volontaire** : l'attente est la mise en scène (images qui défilent,
  barre de progression). Voir [Durée d'analyse variable](#durée-danalyse-variable).

---

## 🏗️ Architecture

```
Biodetective/
├── data_pipeline.py         # [FAIT] new_taxdump/ → biodetective.db
├── common_names_fr.json     # [FAIT] taxid → nom français (maintenu à la main)
├── sequences.json           # [FAIT] séquence → taxid (cœur de la démo)
├── identify.py              # [FAIT] normalisation + recherche exacte
├── validate_sequences.py    # [FAIT] contrôle avant démo
├── make_montages.py         # [FAIT] fabrique les planches de l'écran de recherche
├── biodetective.db          # [GÉNÉRÉ] 25 545 organismes, 11 Mo
├── requirements.txt         # [FAIT]
├── api.py                   # [FAIT] API FastAPI, 8 routes
├── package.json             # [FAIT] config React (proxy → localhost:8000)
├── public/montages/         # [GÉNÉRÉ] 20 planches JPEG, 18,6 Mo, non versionnées
├── new_taxdump/             # données NCBI (non versionnées)
│   ├── names.dmp            # noms scientifiques et communs anglais
│   ├── nodes.dmp            # rang de chaque taxon
│   ├── images.dmp           # images — point d'entrée du pipeline
│   ├── rankedlineage.dmp    # lignée éclatée par rang
│   └── … (11 autres .dmp inutilisés)
├── src/                     # [FAIT] frontend React
│   ├── index.js
│   ├── App.js
│   ├── BioDetective.js      # composant unique, 5 écrans
│   └── biodetective.css
└── public/
    ├── index.html
    ├── favicon.svg
    └── fonts/               # Orbitron + Rajdhani en local + fonts.css
```

---

## 🔧 Stack technique

### Backend
- **Python 3.9** (version installée sur la machine) avec **FastAPI** et **Uvicorn**
- **SQLite** via le module `sqlite3` natif
- Identification par table de correspondance — **aucune dépendance bioinfo externe**
- **Encodage UTF-8 explicite** partout (noms d'organismes avec accents)

### Frontend
- **React 18** avec Create React App (`react-scripts 5.0.1`)
- **Node 22 / npm 10** (versions installées)
- **CSS externe** dans `biodetective.css` (pas de CSS-in-JS)
- **Hooks uniquement** : `useState`, `useEffect`, `useRef`
- **Pas de bibliothèque UI externe** (Tailwind, MUI, etc.)
- **Proxy** vers `http://localhost:8000` configuré dans `package.json`

### Base de données
- **SQLite** : `biodetective.db`
- **25 545 organismes** avec au moins une image (chiffre mesuré dans `images.dmp`)
- Images servies depuis des **URLs distantes NCBI** (voir [Images](#images))
- Taxonomie NCBI : Règne → Phylum → Classe → Ordre → Famille → Genre → Espèce

---

## 🚀 Démarrage

### Prérequis
```bash
# Python 3.9+
pip install -r requirements.txt

# Node.js 22
npm install
```

Pas de BLAST+, pas de conda bioconda, pas de base `nt` à télécharger.

### Construction des données (une seule fois)
```bash
python3 data_pipeline.py     # new_taxdump/*.dmp → biodetective.db
python3 make_montages.py     # → public/montages/, ~2 min, réseau requis
```
Les deux produisent des fichiers volumineux et régénérables, donc non versionnés.

### Lancement — le jour de la démo (une seule commande)
```bash
npm run build        # une fois, après toute modification du frontend
python3 api.py
# → tout sur http://localhost:8000
```
FastAPI sert lui-même le frontend compilé. **Un seul processus, un seul port,
aucun proxy.** C'est le mode à utiliser en salle : rien à expliquer à qui
redémarre la machine.

### Lancement — en développement (rechargement à chaud)
```bash
# Terminal 1
python3 api.py       # http://localhost:8000

# Terminal 2
npm start            # http://localhost:3000, proxy /api → 8000
```
Le port 3000 recharge à chaud mais impose le contrôle d'hôte de react-scripts
(voir `.env`). Le port 8000 n'a pas cette contrainte.

---

## 📡 API Backend

### Endpoints

**Toutes les routes sont préfixées par `/api`** : la racine `/` est occupée par
le frontend compilé.

| Méthode | Route | Description |
|---------|-------|-------------|
| GET | `/api/` | Statut de l'API |
| GET | `/api/stats` | Statistiques de la base de données |
| GET | `/api/random-images?count=8` | Images aléatoires pour l'animation de recherche |
| POST | `/api/analyze` | Soumettre une séquence ADN pour analyse |
| GET | `/api/analyze/{job_id}` | Récupérer le résultat d'un job |
| DELETE | `/api/analyze/{job_id}` | Supprimer un job |
| GET | `/api/organism/{taxonomy_id}` | Détails d'un organisme |
| POST | `/api/reload` | Recharger `sequences.json` sans redémarrer |
| GET | `/` | Frontend compilé (ou message d'aide si `build/` absent) |
| GET | `/docs` | Documentation interactive FastAPI |

Le montage statique est déclaré **après** `include_router` : FastAPI résout les
routes dans l'ordre de déclaration, donc `/api/...` est traité avant le
catch-all du frontend.

`POST /reload` permet de corriger ou d'ajouter une séquence en pleine journée de démo
sans couper le service.

> Les routes s'appellent `/analyze` et non `/blast` : il n'y a pas de BLAST derrière.

### Format requête
```json
POST /api/analyze
{
  "sequence": "ATCGATCGATCG..."
}
```

### Format réponse — correspondance trouvée
```json
{
  "job_id": "uuid-v4",
  "status": "completed",
  "matched": true,
  "organism": {
    "scientific_name": "Escherichia coli",
    "taxonomy_id": 562,
    "percent_identity": 100.0,
    "has_image": true,
    "image_data": {
      "url": "http://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/123",
      "source": "Wikimedia Commons",
      "license": "CC BY-SA 4.0",
      "attribution": "Nom du photographe"
    },
    "organism_type": "Bactérie",
    "kingdom": "Bacteria",
    "common_name_fr": null,
    "common_name_en": "E. coli"
  },
  "analysis_time": 3.2,
  "error_message": null
}
```

### Format réponse — aucune correspondance
```json
{
  "job_id": "uuid-v4",
  "status": "completed",
  "matched": false,
  "organism": null,
  "analysis_time": 3.2,
  "error_message": null
}
```

`matched: false` **n'est pas une erreur** : c'est le cas le plus fréquent (voir
[Écran « Séquence inconnue »](#écran--séquence-inconnue-)). `status: "error"` est réservé
aux vraies pannes (base illisible, séquence vide, caractères non-ADN).

### Statuts de job
`running` → `completed` (avec `matched` true ou false) ou `error`.

Une séquence **invalide** (lettres autres que ATCG, champ vide) ne crée pas de job :
`POST /api/analyze` répond directement **400** avec un message affichable tel quel.
Une séquence **valide mais inconnue** crée bien un job qui aboutit à
`200 / completed / matched:false`. Ne pas confondre les deux côté React.

### `display_name`
Chaque organisme renvoyé porte un champ `display_name` : nom français s'il existe,
sinon nom commun anglais, sinon nom scientifique. Le frontend affiche ce champ sans
avoir à arbitrer lui-même.

### Durée d'analyse variable

La recherche en table est instantanée, mais l'écran de recherche a besoin de
durer. La durée est **tirée au hasard** à chaque analyse, pour que deux enfants
qui se suivent ne voient pas la même chose et qu'une « analyse difficile »
arrive de temps en temps.

| Poids | Durée | Registre |
|-------|-------|----------|
| 55 % | 2,5 – 4 s | rapide |
| 30 % | 4,5 – 7 s | normale |
| 12 % | 7,5 – 10 s | approfondie |
| 3 % | 10,5 – 13 s | très longue |

Moyenne ≈ 4,9 s. Réglable dans `ANALYSIS_TIERS` en tête d'`api.py`.

#### ⚠️ La durée ne doit jamais dépendre du résultat

Le tirage a lieu **avant** la recherche et n'utilise ni la séquence ni le
résultat. Si les analyses longues aboutissaient plus souvent à une séquence
inconnue, l'opérateur — puis les enfants — apprendraient à lire la réponse avant
l'écran de résultat, et toute la mise en scène tomberait.

Vérifié sur 60 analyses : moyenne 5,19 s quand un organisme est trouvé,
5,33 s quand la séquence est inconnue. Écart 0,14 s, soit rien.
**Refaire cette mesure après toute modification d'`ANALYSIS_TIERS`.**

#### Messages par phases

Une attente de 10 s doit avoir l'air de chercher plus profond, pas d'être
bloquée. `SEARCH_PHASES` dans `BioDetective.js` change de registre au fil des
secondes : messages normaux, puis « Séquence complexe, analyse approfondie… »
à partir de 4,5 s, puis « Encore un instant, ça vient… » à partir de 8 s.

#### Accélérer en pleine journée

Si la file d'attente s'allonge, fixer une durée unique sans toucher au code :
```bash
BIODETECTIVE_FIXED_DELAY=3 python3 api.py
```

## 🔑 Table de correspondance (`sequences.json`)

C'est le fichier le plus important du projet : il décide de ce que voient les enfants.
**33 organismes** y sont définis, tous vérifiés (taxid existant, image joignable).

```json
{
  "sequences": [
    {
      "sequence": "AACTGCAGCAAGGTAT",
      "taxonomy_id": 3702,
      "nom_fr": "Arabette des dames",
      "nom_scientifique": "Arabidopsis thaliana",
      "type": "Plante à fleurs"
    }
  ]
}
```

Une liste d'objets plutôt qu'un dictionnaire plat : `nom_fr` permet de relire et corriger
le fichier sans avoir à résoudre les taxids de tête.

### Règles
- **Séquences de 16 bases**, jamais plus de **2 briques identiques d'affilée**. Cette
  contrainte n'est pas cosmétique : devant quatre briques identiques un enfant en compte
  trois ou cinq, et la séquence saisie ne correspond plus à rien.
- **Normalisation avant comparaison** (`identify.normalize`) : majuscules, suppression des
  espaces et retours à la ligne, en-tête FASTA `>` ignoré, `U` d'ARN converti en `T`.
- Seuls `A`, `T`, `C`, `G` sont acceptés ensuite ; sinon `SequenceError` dont le message est
  directement affichable (« Cette séquence contient des lettres qui ne sont pas de l'ADN : … »).
- Chaque `taxonomy_id` **doit** exister en base **et** avoir une image → `validate_sequences.py`.

### Composition de la démo
Plantes 8 (dont *Arabidopsis thaliana* et *Physcomitrium patens*, les modèles de l'IBMP),
mammifères 8, oiseaux 3 (dont la cigogne blanche, clin d'œil alsacien), insectes 4,
monde marin 5, plus amanite, salamandre, *E. coli*, levure et *Chlamydomonas*.

### Avant chaque démo
```bash
python3 validate_sequences.py --images
```
Vérifie les doublons, les runs de briques, les taxids, les images, et signale les images
lourdes. **5 images dépassent 1 Mo** (salamandre 5,6 Mo, tortue verte 2,7 Mo, séquoia 1,7 Mo,
amanite 1,5 Mo, pieuvre 1,2 Mo) : à précharger, sinon l'affichage traîne.

---

## 🗄️ Base de données SQLite

### Table `organisms`
```sql
CREATE TABLE organisms (
    taxonomy_id INTEGER PRIMARY KEY,    -- NCBI Taxonomy ID
    scientific_name TEXT NOT NULL,       -- Nom scientifique binomial
    common_name_fr TEXT,                 -- Nom commun français (33 renseignés, voir plus bas)
    common_name_en TEXT,                 -- Nom commun anglais (depuis names.dmp)
    kingdom TEXT,                        -- Règne
    phylum TEXT,                         -- Embranchement
    class_name TEXT,                     -- Classe
    order_name TEXT,                     -- Ordre
    family TEXT,                         -- Famille
    genus TEXT,                          -- Genre
    species TEXT,                        -- Espèce
    organism_type TEXT,                  -- Type simplifié (ex: "Bactérie", "Plante")
    size_info TEXT,                      -- VIDE
    habitat TEXT,                        -- VIDE
    danger_level TEXT,                   -- VIDE
    discovery_info TEXT,                 -- VIDE
    fun_facts TEXT,                      -- VIDE
    ecological_role TEXT,                -- VIDE
    human_utility TEXT,                  -- VIDE
    image_path TEXT,                     -- URL image (NCBI)
    image_source TEXT,                   -- Source normalisée (Wikimedia Commons, iNaturalist)
    image_license TEXT,                  -- Licence (CC BY-SA 4.0, etc.)
    image_attribution TEXT,              -- Auteur de la photo (souvent vide)
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_scientific_name ON organisms(scientific_name);
CREATE INDEX idx_organism_type ON organisms(organism_type);
CREATE INDEX idx_kingdom ON organisms(kingdom);
```

### ⚠️ Colonnes descriptives vides

`size_info`, `habitat`, `danger_level`, `discovery_info`, `fun_facts`, `ecological_role`,
`human_utility` **resteront NULL** dans un premier temps. Aucune de ces informations n'existe
dans le taxdump NCBI — il faudrait les écrire à la main ou interroger Wikidata/EOL.

Les colonnes sont créées dès maintenant pour ne pas avoir à migrer le schéma plus tard,
mais **le composant React doit masquer toute ligne dont la valeur est NULL ou vide**.
Une fiche avec sept champs « — » est pire que pas de fiche.

### `common_name_fr` vient d'un fichier à part
`names.dmp` ne contient quasiment que des noms communs anglais (`genbank common name`) :
**0 nom français** sur 25 545 organismes. Les noms français sont donc maintenus à la main
dans `common_names_fr.json` (taxid → nom), relu par `data_pipeline.py` à chaque construction.
Actuellement **33 renseignés** — les organismes de la démo. En ajouter = éditer ce fichier
et relancer le pipeline.

`common_name_en` est en revanche rempli pour **8 471 organismes** (33 %), donc l'anglais
peut servir de repli quand le français manque.

### Table `phylogenetic_tree` (non prioritaire)
```sql
CREATE TABLE phylogenetic_tree (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name_fr TEXT,
    name_en TEXT,
    parent_id INTEGER REFERENCES phylogenetic_tree(id),
    level INTEGER,      -- 1=Règne, 2=Phylum, etc.
    color_code TEXT,    -- Couleur hex pour affichage
    icon TEXT,          -- Emoji ou icône
    description_fr TEXT
);
```

---

## 🎨 Frontend React

### Flux de l'application
```
[Écran Accueil]
  → Saisie séquence ADN
  → Clic "LANCER L'ANALYSE"
      ↓
[Écran Recherche] (3-5 s)
  → Images aléatoires d'organismes qui défilent
  → Messages de progression
  → Polling /api/analyze/{job_id} toutes les 500 ms
      ↓
      ├── matched: true  → [Écran Résultat]
      ├── matched: false → [Écran Séquence inconnue]
      └── status: error  → [Écran Erreur]
```

### Écran Résultat
- Photo de l'organisme identifié
- « 100 % » en gros
- Fiche descriptive (type, règne, nom commun) — **lignes vides masquées**
- Crédit photo (source + licence + auteur) — obligatoire pour les licences CC
- Bouton « Nouvelle Analyse »

### Écran « Séquence inconnue »

**Ce sera l'écran le plus vu de la journée.** Les enfants assemblent les briques librement,
donc la plupart des séquences ne seront pas dans `sequences.json`.

Il doit être **valorisant, jamais une erreur** :
- Ton positif : « Séquence inconnue de nos bases ! Tu viens peut-être de découvrir
  une espèce nouvelle 🧬 » — surtout pas « Échec » ou « Introuvable ».
- Réutiliser l'esthétique futuriste, pas un bandeau rouge.
- Afficher la séquence saisie, joliment formatée par blocs de 10 bases.
- Bouton « Réessayer » bien visible.
- Idée à creuser : proposer une photo d'organisme au hasard comme « espèce voisine imaginaire ».

### États React principaux
```javascript
const [currentScreen, setCurrentScreen] = useState('home');  // 'home'|'search'|'result'|'unknown'|'error'
const [sequence, setSequence] = useState('');
const [searchImages, setSearchImages] = useState([]);
const [result, setResult] = useState(null);
const [progressMessage, setProgressMessage] = useState('');
const [isAPIConnected, setIsAPIConnected] = useState(false);
```

### Intervals/Timers à gérer (via useRef)
- `pollingInterval` : poll `/api/analyze/{job_id}` toutes les **500 ms** (l'analyse dure 3-5 s,
  un polling à 2 s ajouterait jusqu'à 2 s d'attente inutile)
- `imageRotationInterval` : rotation des images toutes les 300 ms
- `messageInterval` : rotation des messages toutes les 2 s
- **Important** : tous doivent être nettoyés par `clearAllIntervals()` à chaque sortie de
  l'écran de recherche — résultat, inconnu **et** erreur.

### Gestion d'erreur images
```jsx
// Toujours inclure un fallback onError sur les images
<img
  src={imageUrl}
  onError={(e) => { e.target.src = PLACEHOLDER_DATA_URI; }}
/>
```

---

### Deux pièges à ne pas réintroduire

**1. Le test de vie ne doit pas taper sur `/`.** La racine sert le frontend, dans
les deux modes. Un `fetch('/')` répond 200 avec du HTML même backend éteint : le
test réussirait toujours. On interroge `/api/stats` et on vérifie que la réponse
contient bien un `organisms` numérique.

**2. Le contrôle d'hôte de react-scripts répond HTTP 200.** Son corps vaut
« Invalid Host header », mais le statut est 200. Tester un accès en regardant le
code de retour ne prouve donc rien : il faut regarder le **contenu**. Ce piège ne
concerne que le port 3000 ; le port 8000 n'a aucun contrôle d'hôte.

**3. Pas de prop `key` sur l'image qui tourne.** Avec une `key` changeante, React
remonte un `<img>` neuf toutes les 300 ms ; un élément fraîchement monté n'a pas
encore peint, et le cadre de l'écran de recherche reste noir pendant toute
l'analyse. On réutilise le même noeud en ne changeant que `src`, et on ne fait
tourner que les images **déjà préchargées** (`ready`), jamais le pool brut.

---

### Planches de l'écran de recherche

L'attente doit ressembler à une recherche de série policière, pas à un diaporama.
Les images ne sont donc pas chargées une par une : `make_montages.py` fabrique
**20 planches** de 30 vignettes carrées (600 organismes distincts, 18,6 Mo), et le
frontend en tire une au hasard à chaque analyse.

Le défilement est **entièrement en CSS**, sans aucun timer JavaScript :

```css
.montage__strip { animation: montage-slide var(--durée) steps(30) infinite; }
@keyframes montage-slide { from { translateX(0) } to { translateX(-100%) } }
```

`steps(N)` sur une translation de -100 % tombe exactement sur chaque vignette,
puisque la bande fait N fois la largeur du cadre : le pas vaut `largeur/N`, soit
une vignette pile. Vérifié à la mesure — translations de 836, 1254, 1672 px pour
un cadre de 418 px, tous multiples exacts.

**Cadence : `FRAME_MS = 80` dans `BioDetective.js`** — 12 images par seconde.
C'est le réglage de l'effet. Le modifier ne demande **aucune refabrication**,
contrairement à un GIF dont la vitesse est figée à la création. C'est la raison
principale du choix de la planche plutôt que du GIF ; l'autre est le poids :
18,6 Mo contre 59,5 Mo pour le même contenu en GIF, et sans tramage 256 couleurs.

À 12 images/seconde, `prefers-reduced-motion` **fige la bande** sur sa première
vignette plutôt que de la ralentir. Vérifié.

Si `public/montages/` est absent, le frontend retombe sur une rotation image par
image via `/api/random-images` : dégradé mais fonctionnel.

---

## 🎨 Design System

### Palette de couleurs
```css
--color-primary: #00ff88;      /* Vert néon - couleur principale */
--color-secondary: #00ccff;    /* Bleu cyan - accents */
--color-accent: #ff6b00;       /* Orange - alertes, bouton retour */
--color-bg: #0c0c0c;           /* Fond très sombre */
--color-bg-card: rgba(0, 0, 0, 0.8);  /* Fond des cartes */
--color-text: #ffffff;          /* Texte secondaire */
```

### Briques d'ADN
Les séquences sont affichées en briques colorées (`DnaStrip`), pour que l'enfant
retrouve à l'écran ce qu'il a dans les mains. **Les couleurs de `BASE_COLORS` dans
`BioDetective.js` doivent être ajustées à celles des vraies briques du Brickopore.**

### Polices
```css
font-family: 'Orbitron', monospace;   /* Titres - futuriste */
font-family: 'Rajdhani', sans-serif;  /* Corps - lisible */
font-family: 'Courier New', monospace; /* Séquences ADN */
```

> Orbitron et Rajdhani viennent de Google Fonts. Elles doivent être **téléchargées dans
> `public/fonts/`** et servies en local : la salle n'a pas de connexion garantie et une
> police manquante casse toute l'identité visuelle.

### Animations clés
- `pulse` : effet pulsation sur les titres
- `fadeInOut` : images de recherche qui clignotent
- `progress` : barre de progression infinie pendant l'analyse
- `blink` : statut « en cours »
- `glitch` : transition vers l'écran résultat

---

## 🧬 Données NCBI

### Source
**new_taxdump.tar.gz** (déjà téléchargée et décompressée dans `new_taxdump/`).
Format général : séparateur `\t|\t`, fin de ligne `\t|\n`.

| Fichier | Usage |
|---------|-------|
| `names.dmp` | Noms scientifiques et communs anglais |
| `nodes.dmp` | Hiérarchie taxonomique (parent + rang) |
| `images.dmp` | Images — **le point d'entrée du pipeline** |
| `rankedlineage.dmp` | Lignée déjà éclatée par rang (évite de remonter `nodes.dmp` à la main) |

Les 11 autres `.dmp` (`citations`, `host`, `typematerial`, `taxidlineage`, …) ne servent pas.

### `images.dmp` — format réel
```
image_id | image_key | url | license | attribution | source | properties | taxid_list
```
- **25 570 lignes**, **25 545 taxonomy_ids uniques**
- `taxid_list` : taxids séparés par des espaces (une image peut couvrir plusieurs taxons)

### ⚠️ Les URLs ne sont pas des URLs Wikimedia

La colonne `url` contient des adresses **hébergées par le NCBI**, pas des liens directs :

```
http://www.ncbi.nlm.nih.gov/Taxonomy/taxi/images/3
```

Deux conséquences :
1. **C'est du `http://`, pas du `https://`.** Ça fonctionne parce que l'app est servie sur
   `http://localhost:3000` — aucun blocage de contenu mixte. **Si l'app passe un jour derrière
   HTTPS, toutes les images disparaissent.**
2. Le rendu dépend de la disponibilité du NCBI le jour J. « Wikimedia Commons » n'est que
   l'attribution de la colonne `source`.

### ⚠️ La colonne `source` doit être normalisée
Valeurs réellement présentes : `iNaturalist` (1 611), `iNaturalist.com` (1 648),
`inaturalist.com` (3), `inaturalist.org` (3) — quatre orthographes d'une même source.
Idem `Wikimedia Commons` (22 130) vs `Wikimedia  Commons` (double espace, 2) vs
`wikimedia.org` (1). Plus quelques marginaux (`expasy.org`, `CDC.gov`, `WoRMS`…).
**Normaliser dans `data_pipeline.py`**, pas à l'affichage.

### Pipeline de données (`data_pipeline.py`)
1. Charger `images.dmp` → ensemble des taxids ayant une image (+ url, licence, attribution, source normalisée)
2. Charger `names.dmp` (nom scientifique + nom commun anglais) et `rankedlineage.dmp` (lignée)
3. Pour chaque taxid avec image → construire un `Organism`
4. Déduire `organism_type` (« Bactérie », « Plante », « Champignon », « Animal »…) depuis le règne/phylum
5. `INSERT` dans `biodetective.db`

> `names.dmp` (280 Mo), `nodes.dmp` (282 Mo) et `rankedlineage.dmp` (393 Mo) sont volumineux :
> filtrer d'abord sur les 25 545 taxids d'`images.dmp`, ne pas tout charger en mémoire.

---

## ⚠️ Points d'attention

### Encodage
- SQLite : `conn.execute("PRAGMA encoding = 'UTF-8'")`
- Logging : `FileHandler('biodetective.log', encoding='utf-8')`
- Ouvrir tous les `.dmp` avec `encoding='utf-8'` explicite
- Les noms d'organismes contiennent des accents, ñ, ü, etc.

### Jobs
- Chaque job a un UUID unique stocké **en mémoire** dans un dict
- Les jobs ne sont pas persistés : si le backend redémarre, les jobs en cours sont perdus
- Prévoir une purge des jobs terminés (sinon fuite mémoire sur une journée entière de démo)

### Images
- URLs distantes → toujours gérer `onError` côté React
- L'endpoint `/random-images` fait `ORDER BY RANDOM()` : appelé à chaque analyse, donc
  **charger une fois au démarrage un pool de quelques centaines d'images en mémoire** et y piocher

### Le jour de la démo
- Vérifier la connexion internet **avant** l'ouverture au public (les images en dépendent)
- Lancer le script de validation de `sequences.json` (tous les taxids présents et avec image)
- Prévoir un jeu de séquences imprimées connues pour relancer la démo si un enfant bloque

---

## 🚧 Fonctionnalités à développer

- [ ] **Script de validation `sequences.json`** — vérifier que chaque taxid existe et a une image
- [ ] **Polices en local** — retirer la dépendance Google Fonts
- [ ] **Correspondance approximative** — si un enfant se trompe d'une brique, retrouver quand même
      l'organisme (distance d'édition sur `sequences.json`). Transformerait la majorité des
      « séquence inconnue » en résultats. **À arbitrer** : contredit le choix du lookup exact.
- [ ] **Fiches descriptives** — remplir `fun_facts` & co. à la main pour les organismes de `sequences.json`
- [ ] **Noms communs français** — enrichissement via Wikidata API
- [ ] **Images en local** — supprimer la dépendance au NCBI le jour J
- [ ] **Arbre phylogénétique interactif** — SVG/D3.js avec la table `phylogenetic_tree`
- [ ] **Historique des analyses** — garder les dernières séquences analysées
- [ ] **Son/musique** — ambiance sonore pendant la recherche

---

## 🗺️ Roadmap (ordre de construction)

1. ~~`git init` + `.gitignore`~~ ✅
2. ~~`data_pipeline.py` → `biodetective.db`~~ ✅
3. ~~`sequences.json` + `identify.py` + `validate_sequences.py`~~ ✅
4. ~~`api.py` : `/`, `/stats`, `/organism/{id}`~~ ✅
5. ~~`api.py` : `/analyze` + polling + `/random-images`~~ ✅
6. ~~Frontend : accueil → recherche → résultat~~ ✅
7. ~~Frontend : écran « séquence inconnue »~~ ✅
8. ~~Polices locales~~ ✅, ~~purge des jobs~~ ✅, ~~serveur unique~~ ✅, cache d'images **local** (reste à faire)
9. Répétition générale dans les conditions réelles de la salle

---

## 🎯 Contexte d'utilisation

- **Public** : enfants et parents lors de la Fête de la Science
- **Lieu** : IBMP Strasbourg, salle de présentation
- **Matériel** : ordinateur fixe avec écran large (pas de mobile)
- **Réseau** : internet disponible le jour J (confirmé) — **les images en dépendent**
- **Opérateur** : David (bioinformaticien), seul à utiliser l'interface
- **Séquences** : courtes (20-100 bases), assemblage **libre** des briques par les enfants
  → la plupart des séquences ne seront pas reconnues, c'est normal et prévu
- **Durée d'une démo** : 2-5 minutes par enfant

---

## 👤 Développeur

- **David** - Bioinformaticien, IBMP Strasbourg
- Environnement : Rocky Linux 9, HPC cluster « babibel », conda (env `base`)
- Expertise : Python, bioinformatique, génomique des plantes (Arabidopsis), Nanopore
- Niveau React/JS : intermédiaire (projet en apprentissage)
