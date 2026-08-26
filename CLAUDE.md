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

### Identification par BLAST réel

L'identification passe par un **vrai blastn** contre une banque construite à
partir de `nt`, restreinte aux espèces qui possèdent une image. On obtient donc
un classement par E-value, un pourcentage d'identité réel et un alignement à
montrer — bien plus parlant qu'un « 100 % » décrété.

**Pourquoi une banque restreinte, et pas `nt` :** les brins du Brickopore font
24 bases. Contre `nt` (~10¹² lettres) une requête aussi courte n'a aucune chance :
la E-value attendue dépasse 1. Contre une banque de ~25 Mpb, elle tombe vers
10⁻⁶. **La E-value est proportionnelle à la taille de la banque** — c'est la
contrainte qui gouverne tout le reste, et la raison pour laquelle
`extract_from_nt.py` ne prend qu'une séquence par espèce et rejette les génomes
complets.

Mesuré sur une base simulée de 16 Mpb :

| Longueur | Meilleur hit | E-value | Correct ? |
|----------|--------------|---------|-----------|
| 12 pb | une espèce au hasard | 1,6 | ✗ |
| 16 pb | la bonne | 0,010 | ✓ |
| **24 pb** | la bonne | **5,4 × 10⁻⁷** | ✓ |
| 30 pb | la bonne | 2,3 × 10⁻¹⁰ | ✓ |

12 briques est sous le plancher : une séquence sans rapport gagne. 24 est le
compromis retenu entre impact à l'écran et nombre de briques à assembler.

**Ce que BLAST apporte en plus :** une brique mal comptée sur 24 donne encore le
bon organisme, à 95,8 % d'identité. La table de correspondance exacte, elle,
échouait. C'était le principal risque de la démonstration.

**Séquences inventées :** sur 8 brins de 24 briques tirés au hasard, aucun ne
produit de hit jugé sûr. L'écran « séquence inconnue » reste donc le cas normal
pour un enfant qui improvise.

### Chaîne de construction de la banque

**Avec nt copié sur la machine** — la voie retenue, et de loin la plus simple :

```bash
blastdb_aliastool -db nt -taxidlist taxids.txt -dbtype nucl \
                  -out blastdb/biodetective -title "BioDetective"
python3 make_strips.py --length 24
python3 validate_sequences.py
```

`blastdb_aliastool` écrit un simple fichier `.nal` **de quelques centaines
d'octets** qui restreint nt à nos taxons. Aucune extraction, aucune copie de
séquences, et surtout **l'espace de recherche est réellement réduit** : blastn
annonce la taille du sous-ensemble, donc les E-value sont justes.

C'est ce dernier point qui départage les deux méthodes. `blastn -db nt
-taxidlist …` filtre bien les résultats mais **calcule les statistiques sur nt
entier** : vérifié, la taille annoncée reste inchangée, et un brin de 24 briques
retomberait vers E ≈ 0,1. L'alias, lui, donne 1,0 × 10⁻⁶.

⚠️ Copier **aussi les fichiers `taxdb.btd` et `taxdb.bti`** à côté de nt. Sans
eux, `-taxidlist` refuse de fonctionner et `staxids` ne rend que des zéros.

`extract_from_nt.py` et `prepare_blastdb.py` restent utiles si nt n'est pas
accessible localement, mais l'alias les rend inutiles dès qu'il l'est.

### La banque en service

Provisoire, issue de la première extraction plafonnée : **44 772 séquences,
21 489 espèces, 32,0 Mpb**. Elle fait tourner les 31 brins préparés, mais ne
contient qu'un marqueur par espèce — une séquence prise ailleurs dans un génome
n'y trouve rien. À remplacer par l'extraction sans plafond.

### Brins préparés : 32 sur 33

Un seul organisme résiste : le **papillon monarque**. BLAST y arrive, mais tous
ses fragments discriminants contiennent 3 briques identiques d'affilée, ce que la
règle LEGO interdit.

La **vigne**, longtemps impossible, est repassée dès que le chloroplaste de
*Vitis vinifera* a rejoint la banque : la séquence de référence d'origine était
partagée avec les autres vignes, le chloroplaste ne l'est pas. Illustration de la
règle générale — **plus la banque est riche, plus il y a de brins possibles**.

`make_strips.py` sait désormais tirer sa séquence de référence de la banque
BLAST elle-même, via `blastdbcmd`, quand aucun FASTA n'est présent sur le
disque. C'est indispensable avec un alias construit au-dessus de nt, où il n'y a
aucun FASTA à lire.

`taxids.txt` se régénère par
`sqlite3 biodetective.db "SELECT taxonomy_id FROM organisms;" > taxids.txt`.

⚠️ `blastn` vient de conda alors que `python3` est celui du système : il n'est
donc pas dans le PATH. `blast_search.find_blastn()` va le chercher dans les
emplacements conda usuels ; `BIODETECTIVE_BLASTN` force un chemin.

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
Le port 3000 recharge à chaud mais impose le contrôle d'hôte de react-scripts.

⚠️ **`.env.local` n'est pas versionné** et il est nécessaire pour joindre le
port 3000 par un nom d'hôte plutôt que par son IP. Le recréer après un clone :

```bash
echo 'DANGEROUSLY_DISABLE_HOST_CHECK=true' > .env.local
```

react-scripts n'autorise qu'un seul nom d'hôte, l'IP LAN qu'il détecte ; tout
autre nom reçoit « Invalid Host header ». Ce réglage lève une protection contre
le DNS rebinding : acceptable pour un serveur de développement sur le réseau de
l'institut, à ne pas reprendre ailleurs. Le port 8000 n'est pas concerné. Le port 8000 n'a pas cette contrainte.

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
| GET | `/api/history?limit=12&matched_only=` | Dernières analyses |
| GET | `/api/history/{index}` | Détail d'une analyse passée, pour la rejouer |
| DELETE | `/api/history` | Vider l'historique |
| POST | `/api/remote-blast` | BLAST distant chez le NCBI (lent, hors démo) |
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
- Crédit photo (source + licence + auteur) — **en surimpression, au survol
  seulement** : il alourdissait un écran qui doit d'abord montrer un organisme.
  Positionné hors du flux, pour que son apparition ne déplace rien. `focus-within`
  et `tabindex` le rendent atteignable au clavier, faute de souris.
  Les licences CC demandent l'attribution : la garder accessible, jamais la retirer
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

### Taille du génome et nombre de chromosomes

La fiche annonce la taille du génome et le nombre de chromosomes, quand on les
connaît. Ces chiffres ne sont **pas** dans le taxdump : ils viennent des
assemblages, via les rapports NCBI `GENOME_REPORTS`.

```bash
python3 build_genome_stats.py --download   # -> genome_stats.json
python3 data_pipeline.py                   # relit le JSON
```

Les rapports bruts (plus de 230 Mo) ne sont pas versionnés ; `genome_stats.json`,
compact, l'est. `data_pipeline.py` le relit à chaque construction, comme
`common_names_fr.json`.

**Couverture : 4 403 espèces sur 25 545 (17 %), mais 30 des 33 de la
démonstration.** La plupart des espèces photographiées n'ont jamais été
séquencées : les lignes vides sont simplement masquées.

⚠️ Les trois rapports **n'ont pas la même disposition de colonnes** : dans
`prokaryotes.txt` l'accession est en 19ᵉ colonne et les réplicons en 9ᵉ,
l'inverse d'`eukaryotes.txt`, et `viruses.txt` donne des **Kb** et non des Mb.
Les colonnes sont donc repérées **par leur intitulé**, jamais par leur rang.
Lire au rang donnait 11 Mb pour *Escherichia coli*, qui en fait 4,6.

Une espèce a souvent des dizaines d'assemblages : on retient celui qui porte
une vraie accession, puis des réplicons RefSeq (`NC_`/`NZ_`), c'est-à-dire le
génome de référence. `genome_stats.json` conserve aussi la correspondance
`chromosome → accession`, qui permettra de situer un hit sur le bon chromosome
et de le dessiner.

---

### Écran « L'arbre du vivant »

Accessible depuis l'accueil. Trois domaines, quelques grands groupes sous les
eucaryotes, les virus à part, avec une vignette par groupe tirée de la base.

Chaque groupe porte **deux chiffres** : les espèces **décrites par la science**,
et celles que nous avons **en photo**. Le premier vient de la taxonomie NCBI
entière, pas de notre collection — sinon l'arbre montrerait la forme de ce qui a
été photographié, pas celle du vivant.

| Groupe | Décrites | En photo |
|--------|---------:|---------:|
| Bactéries | 551 513 | 76 |
| Archées | 13 638 | 9 |
| Eucaryotes | 1 754 746 | 25 274 |
| — Animaux | 1 151 891 | 14 180 |
| — Plantes | 348 210 | 9 110 |
| — Champignons | 201 049 | 1 842 |
| — Algues et protistes | 53 596 | 142 |
| Virus | 80 011 | 186 |
| **Total** | **2 436 560** | **25 545** |

L'écart est le propos : 551 513 bactéries décrites, 76 photographiées. Ce que
l'on voit du vivant n'est pas ce qu'il est.

```bash
python3 build_tree_stats.py     # -> tree_stats.json, ~20 s
```

Ne comptent que les taxons de **rang `species`** : la taxonomie contient aussi
les genres, familles et sous-espèces, et les additionner gonflerait les chiffres
sans rien vouloir dire. Les virus sont reconnus par les suffixes ICTV, comme
dans `data_pipeline.py`.

#### « Procaryote » n'est pas une branche

Le mot décrit une **cellule sans noyau** ; il ne désigne pas un groupe de
parenté. Bactéries et archées sont deux domaines distincts, aussi éloignés
l'un de l'autre qu'ils le sont de nous — les regrouper sur l'arbre laisserait
croire à une parenté qui n'existe pas. L'écran l'explique dans un encadré
plutôt que de faire du mot une branche.

Les virus sont encadrés en pointillés, sans être rattachés à un domaine.

`superkingdom` a été ajouté à `organisms` pour cet écran : le rang existait
dans `rankedlineage.dmp` mais n'était pas repris. Relancer `data_pipeline.py`
après mise à jour.

---

### « Ses cousins » — la parenté déduite des hits

Une liste de hits BLAST n'est pas une liste de mauvaises réponses : ce sont les
**parents** de l'organisme trouvé. L'écran de résultat les présente donc comme
un arbre, groupés par le rang taxonomique qu'ils partagent avec lui.

Sur le brin « Lion » :

| Rang partagé | Cousins | Identité |
|---|---|---|
| genre *Panthera* | tigre, léopard, jaguar | 100 %, 95,5 %, 95,5 % |
| classe *Mammalia* | souris à poche soyeuse | 100 % |
| embranchement *Chordata* | deux rainettes | 95,5 % |

Le dégradé est la leçon : plus le cousin est proche, meilleur est le hit. C'est
le sens même d'un alignement, rendu visible sans un mot d'explication.

`shared_rank()` compare les deux lignées du plus précis au plus large et retient
le premier rang où elles coïncident. Aucun arbre externe : la taxonomie est déjà
dans `biodetective.db`.

Cet arbre a remplacé l'ancienne liste « autres correspondances », qui faisait
doublon dès lors que les mêmes organismes s'affichaient au-dessus avec leur
photo. Les chiffres — identité et E-value — sont passés sur les vignettes.

---

### Rapport BLAST brut, dans l'interface

`▸ Voir le rapport BLAST brut` déplie le fichier archivé : la commande, le
tableau, les alignements, les paramètres de Karlin-Altschul. C'est la réponse à
« et ça donne quoi, vraiment ? » — la question que pose un collègue
bioinformaticien devant le stand.

Servi par `GET /api/blast-report/{job_id}`, qui vérifie que le chemin demandé
reste bien dans `blast_results/` : une route qui rend un fichier ne doit jamais
pouvoir en rendre un autre. Masqué lors de la relecture d'une analyse passée,
dont le job n'existe plus en mémoire.

---

### Idéogramme

Le chromosome touché est dessiné à l'échelle, avec un marqueur à la position du
hit — un SVG, sans aucune bibliothèque.

> CHROMOSOME 1 — 1 DES 5 CHROMOSOMES · 30,4 Mb
> `1 ────────────┃──────────────── 30,4 Mb`  hit en 12 345 679

Le nom du chromosome vient de la correspondance `accession → nom` de
`genome_stats.json`, pas d'une analyse du titre : un hit sur `NC_003070.9`
devient « Chromosome 1 » de façon sûre. Sans correspondance, on retombe sur
l'étiquette tirée du titre (« Génome chloroplastique »).

**Seul le chromosome atteint est dessiné**, à sa longueur réelle donnée par
`slen`. Représenter les autres supposerait des longueurs dont nous ne disposons
pas : un caryotype complet exigerait la taille de chaque réplicon, que les
rapports NCBI ne donnent pas. Un caryotype inventé vaudrait moins que pas de
dessin du tout.

Rien n'est dessiné en dessous de 1 000 pb : un marqueur de 800 pb n'a pas de
géographie.

---

### Localisation du hit

L'écran de résultat annonce **où** la séquence a été trouvée : type de molécule,
position et brin.

> génome chloroplastique — position 90 001 – 90 034 sur 154 478 bases

La position vient directement de `sstart`/`send` : elle a toujours été dans la
sortie de BLAST. Ce qui manquait, c'était de savoir **dans quoi** : l'extraction
écrivait `>taxid|accession` et **jetait la description**. Le titre est désormais
conservé (`>taxid|accession description`), demandé à BLAST via `stitle`, et
`describe_locus()` en tire une étiquette lisible : « chromosome 1 », « génome
mitochondrial », « génome chloroplastique », « plasmide … ». Sans description
reconnaissable, la ligne n'est simplement pas affichée.

`sstart > send` signale un alignement sur le brin complémentaire, indiqué comme
tel. Vérifié sur le chloroplaste d'*Arabidopsis* : un fragment pris en 90 001
ressort en 90 001–90 034 sur le brin plus, et **une seconde fois en
148 648–148 615 sur le brin moins** — la répétition inversée du chloroplaste.
Le comportement est correct, et l'exemple est joli à montrer.

⚠️ Cela ne fonctionne que si la banque contient des enregistrements
génomiques. Avec l'ancienne extraction plafonnée à un marqueur par espèce, les
positions se rapportaient à un ADNc de 785 pb et n'apprenaient rien.

---

### BLAST distant chez le NCBI

`remote_blast.py` interroge le BLAST public du NCBI par son interface URL
(`Put` → `SearchInfo` → `Get`). En ligne de commande ou par
`POST /api/remote-blast`, dont le résultat se relit avec
`GET /api/analyze/{job_id}`.

```bash
export BIODETECTIVE_NCBI_EMAIL="prenom.nom@exemple.fr"   # le NCBI le demande
python3 remote_blast.py TCATTGTAAGATTGGAATAATTCAATTTCGACAT
```

**C'est un outil de vérification, jamais le moteur de la démonstration**, et la
mesure le montre nettement.

| | Banque locale | NCBI `core_nt` |
|---|---|---|
| Taille | 63 Mpb | 998 **G**pb |
| Durée | ~0,3 s | **33 à 53 s** |
| Brin « Lion », meilleur hit | *Panthera leo* | ***Panthera pardus*** (léopard) |
| E-value | 1,1 × 10⁻⁶ | 0,11 |

Deux raisons de garder la banque locale : la latence, sans commune mesure avec
un écran de recherche de 2,5 à 13 s, et surtout **la justesse**. Vingt-quatre
bases de COI sont partagées par tout le genre *Panthera* : contre `core_nt` le
léopard sort avant le lion. Restreindre la banque aux espèces que l'on sait
illustrer n'est pas qu'une économie, c'est ce qui rend la réponse juste.

Deux pièges rencontrés en écrivant le client :

- `ALIGNMENTS=0&DESCRIPTIONS=0` bornent le **nombre de lignes** rendues, pas
  leur verbosité : à zéro, le résultat revient vide alors que la recherche a
  trouvé.
- `FORMAT_TYPE=Tabular` rend un corps vide sur cette interface. Le **XML** est
  le format fiable, et il porte en prime les séquences alignées et la longueur
  du sujet. Les accessions y sont **sans version** (`AY052207`) là où esummary
  répond avec (`AY052207.1`) : indexer les deux formes, sinon aucun taxid ne
  se résout.

L'adresse de contact n'est transmise que si `BIODETECTIVE_NCBI_EMAIL` est
définie ; rien n'est codé en dur.

#### Pourquoi pas `blastn -remote`

BLAST+ sait faire la même chose tout seul, et ce serait plus simple à écrire.
Essayé : `blastn -task blastn-short -query … -db nt -remote` **n'a jamais rendu
la main en plus de 14 minutes** pour un brin de 24 bases, sans message d'erreur.

Ce n'est pas un blocage réseau — `strace` montre la résolution DNS puis une
connexion en HTTPS vers le NCBI, et le processus dort ensuite (`wchan =
hrtimer_nanosleep`, aucun socket ouvert) : il interroge à un rythme bien plus
lâche que nécessaire. Le même brin revient en **33 s** par l'interface URL.

`remote_blast.py` est donc plus long à écrire, mais vingt fois plus rapide et
il rend la main. À reconsidérer si une version future de BLAST+ resserre sa
cadence d'interrogation.

---

### Traces BLAST vérifiables

Chaque analyse dépose un fichier dans `blast_results/` :

```
blast_results/
├── index.tsv                       # récapitulatif, une ligne par analyse
├── 20260825-164422_d7b426f2.txt
└── …
```

Chaque fichier contient, dans l'ordre :

1. l'en-tête — date, job, requête, banque, **la commande blastn exacte** ;
2. le tableau `outfmt 6` que l'application a réellement lu ;
3. le rapport `outfmt 0` complet, avec alignements et paramètres de
   Karlin-Altschul ;
4. le verdict de l'application, avec les seuils appliqués.

C'est du blastn brut, pas notre restitution : on peut donc contrôler ce que
l'application a décidé, et rejouer la commande à l'identique.

`blastn` n'émet qu'un format à la fois, d'où deux appels. Ils sont **lancés en
parallèle** : 286 ms au lieu de 584 ms en série, soit 17 ms de plus qu'une
analyse sans archivage. Compter ~9 ko par analyse, soit 5 Mo pour 500.

`BIODETECTIVE_BLAST_RESULTS=""` désactive l'archivage. Une erreur d'écriture
n'interrompt jamais une analyse en cours devant un enfant.

---

### Historique des analyses

Les 40 dernières analyses sont conservées dans `history.json`, **écrit sur
disque** : une journée de démonstration est longue et un redémarrage du backend
ne doit pas l'effacer. L'écriture est atomique (fichier temporaire puis
`os.replace`), pour qu'un arrêt brutal ne laisse pas un fichier tronqué. Un
historique illisible est ignoré avec un avertissement plutôt que d'empêcher
l'application de démarrer devant une file d'enfants.

Les séquences inconnues **y figurent aussi** : c'est le cas le plus fréquent, et
l'opérateur veut savoir combien d'enfants sont passés, pas seulement combien ont
réussi. L'écran d'accueil n'affiche en revanche que les organismes identifiés
(`matched_only=true`), les vignettes étant là pour donner envie.

Un clic sur une vignette **réaffiche l'analyse sans relancer BLAST** (~140 ms) :
l'enfant qui revient avec ses parents veut revoir son organisme, pas refaire la
queue. Chaque entrée porte son `index` dans la liste complète, pour que le
filtrage sur les réussites ne décale pas les clics.

Pour repartir de zéro entre deux groupes :
```bash
curl -X DELETE http://localhost:8000/api/history
```

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
retrouve à l'écran ce qu'il a dans les mains. `BASE_COLORS` reprend les teintes
LEGO officielles des briques réellement utilisées par le Brickopore :

| Base | Couleur | Fond | Texte | Contraste |
|------|---------|------|-------|-----------|
| A | bleu | `#0055bf` | blanc | 6,9:1 |
| T | vert | `#4b9f4a` | noir | 5,9:1 |
| G | jaune | `#f2cd37` | noir | 12,6:1 |
| C | rouge | `#c91a09` | blanc | 5,8:1 |

La couleur du **texte change selon la base** : sur le bleu et le rouge, une
lettre noire tombe à 2,8 et 3,4 de contraste. Le vert et le jaune font
l'inverse. Recalculer si une teinte change.

T vert et C rouge sont la paire classiquement confondue par les daltoniens
(deutéranopie, protanopie). C'est déjà le cas des vraies briques, et les
changer trahirait le dispositif : la lettre inscrite sur chaque brique est le
canal de secours, elle doit rester lisible.

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
