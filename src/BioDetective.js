import React, { useCallback, useEffect, useRef, useState } from 'react';

// Couleurs des briques LEGO du Brickopore : A bleu, T vert, G jaune, C rouge.
// Teintes officielles LEGO, pour que l'enfant retrouve à l'écran exactement
// les briques qu'il a dans les mains.
//
// La couleur du texte varie parce que les quatre teintes ne se valent pas :
// sur le bleu et le rouge, une lettre noire tombe à 2,8 et 3,4 de contraste ;
// en blanc elle remonte à 6,9 et 5,8. Le vert et le jaune font l'inverse.
// Tous les couples retenus dépassent 4,5:1.
const BASE_COLORS = {
  A: { bg: '#0055bf', fg: '#ffffff' },  // Bright Blue
  T: { bg: '#4b9f4a', fg: '#0c0c0c' },  // Bright Green
  G: { bg: '#f2cd37', fg: '#0c0c0c' },  // Bright Yellow
  C: { bg: '#c91a09', fg: '#ffffff' },  // Bright Red
};

const UNKNOWN_BASE = { bg: '#555555', fg: '#ffffff' };

// L'analyse dure ~4 s côté serveur ; interroger toutes les 500 ms suffit à
// ce que l'attente perçue ne dépasse pas le délai voulu.
// Les routes de l'API vivent sous /api : la racine est occupée par le
// frontend lui-même, que ce soit react-scripts en développement ou FastAPI
// qui sert le build en production.
const API = '/api';

const POLL_MS = 500;
// Cadence du défilement, en millisecondes par vignette. C'est le réglage
// « effet série policière » : 300 ms faisait diaporama, 80 ms fait recherche.
// Modifiable ici sans refabriquer les planches, contrairement à un GIF.
const FRAME_MS = 80;

// Utilisée seulement en repli, quand les planches n'ont pas été fabriquées.
const IMAGE_ROTATE_MS = 120;
const MESSAGE_ROTATE_MS = 2000;

// En dessous, la saisie est un accident plutôt qu'une séquence.
const MIN_BASES = 4;

// Les analyses durent de 2,5 à 13 s (tirage côté serveur). Une attente longue
// doit avoir l'air de chercher plus profond, pas d'être bloquée : les messages
// changent donc de registre au fil des secondes.
const SEARCH_PHASES = [
  {
    after: 0,
    messages: [
      'Lecture de la séquence ADN…',
      'Comparaison avec les organismes connus…',
      'Analyse des correspondances…',
    ],
  },
  {
    after: 4500,
    messages: [
      'Aucune correspondance évidente…',
      'Élargissement à la base taxonomique complète…',
      'Recoupement des empreintes génétiques…',
      'Séquence complexe, analyse approfondie…',
    ],
  },
  {
    after: 8000,
    messages: [
      'Analyse approfondie en cours…',
      'Vérification des derniers candidats…',
      'Encore un instant, ça vient…',
    ],
  },
];

const SEARCH_MESSAGES = SEARCH_PHASES[0].messages;

// Affichés sur l'écran « séquence inconnue », qui sera de loin le plus vu :
// autant qu'on y apprenne quelque chose.
const DNA_FACTS = [
  "L'ADN de chacun de tes doigts est identique… mais tes empreintes, non !",
  'Déroulé, l’ADN d’une seule de tes cellules mesurerait 2 mètres.',
  'Tu partages environ 60 % de tes gènes avec une banane.',
  'Il y a 4 lettres dans l’ADN : A, T, C et G. Comme 4 couleurs de briques.',
  'Une bactérie recopie son ADN en 20 minutes seulement.',
  'Le plus grand génome connu appartient à une fougère, pas à un humain.',
];

// Image de repli si une URL distante ne répond pas : un point d'interrogation
// dessiné en SVG, donc jamais de requête réseau supplémentaire.
const PLACEHOLDER =
  'data:image/svg+xml;utf8,' +
  encodeURIComponent(
    `<svg xmlns="http://www.w3.org/2000/svg" width="400" height="400">
       <rect width="400" height="400" fill="#111"/>
       <text x="50%" y="50%" fill="#00ff88" font-size="120"
             text-anchor="middle" dominant-baseline="central"
             font-family="monospace">?</text>
     </svg>`
  );

/** Séquence ADN rendue en briques colorées, comme sur le Brickopore. */
function DnaStrip({ sequence, size = 'normal' }) {
  if (!sequence) return null;
  return (
    <div className={`dna-strip dna-strip--${size}`}>
      {sequence.split('').map((base, i) => {
        const c = BASE_COLORS[base] || UNKNOWN_BASE;
        return (
          <span
            key={i}
            className="dna-brick"
            style={{ background: c.bg, color: c.fg }}
            title={`Base ${i + 1} : ${base}`}
          >
            {base}
          </span>
        );
      })}
    </div>
  );
}

/** Sépare les milliers par une espace fine, comme il se doit en français. */
const fmt = (n) =>
  typeof n === 'number' ? n.toLocaleString('fr-FR') : n;

/** Le chromosome touché, dessiné à l'échelle, avec le hit repéré dessus.

   Seul le chromosome atteint est représenté : c'est le seul dont BLAST nous
   donne la longueur exacte. Dessiner les autres supposerait des tailles que
   nous n'avons pas, et un caryotype inventé vaudrait moins que rien. */
function Ideogram({ data }) {
  if (!data) return null;

  const W = 900;          // repère interne du SVG, mis à l'échelle en CSS
  const H = 54;
  const PAD = 26;
  const BAR = 26;
  const barW = W - 2 * PAD;
  const x = PAD + Math.min(1, Math.max(0, data.fraction)) * barW;
  const mb = (n) => `${(n / 1e6).toFixed(1)} Mb`;

  return (
    <div className="ideo">
      <div className="ideo__head">
        <span>
          {data.label}
          {data.chromosome_count && data.is_chromosome
            ? ` — 1 des ${data.chromosome_count} chromosomes`
            : ''}
        </span>
        <span className="ideo__len">{mb(data.length)}</span>
      </div>

      <svg
        className="ideo__svg"
        viewBox={`0 0 ${W} ${H}`}
        preserveAspectRatio="none"
        role="img"
        aria-label={`Position du hit sur ${data.label} : ${data.position} sur ${data.length} bases`}
      >
        <defs>
          <linearGradient id="chrom" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stopColor="#2b2b2b" />
            <stop offset="45%" stopColor="#4a4a4a" />
            <stop offset="100%" stopColor="#232323" />
          </linearGradient>
        </defs>

        <rect
          x={PAD} y={(H - BAR) / 2} width={barW} height={BAR}
          rx={BAR / 2} fill="url(#chrom)"
          stroke="rgba(0,255,136,0.35)" strokeWidth="1"
        />

        {/* Repères tous les dixièmes, pour donner l'échelle sans chiffrer. */}
        {Array.from({ length: 9 }, (_, i) => (
          <line
            key={i}
            x1={PAD + ((i + 1) / 10) * barW} x2={PAD + ((i + 1) / 10) * barW}
            y1={(H - BAR) / 2 + 4} y2={(H + BAR) / 2 - 4}
            stroke="rgba(255,255,255,0.10)" strokeWidth="1"
          />
        ))}

        {/* Le hit : trait plein sur toute la hauteur, plus une pointe. */}
        <line
          x1={x} x2={x} y1={(H - BAR) / 2 - 5} y2={(H + BAR) / 2 + 5}
          stroke="#00ff88" strokeWidth="3" strokeLinecap="round"
        />
        <polygon
          points={`${x - 6},${(H - BAR) / 2 - 6} ${x + 6},${(H - BAR) / 2 - 6} ${x},${(H - BAR) / 2 + 2}`}
          fill="#00ff88"
        />
      </svg>

      <div className="ideo__foot">
        <span>1</span>
        <span className="ideo__hit">
          hit en {fmt(data.position)}
          {data.strand === 'moins' ? ' (brin complémentaire)' : ''}
        </span>
        <span>{mb(data.length)}</span>
      </div>
    </div>
  );
}

/** Alignement BLAST, présenté comme dans un vrai rapport.

   Chaque base de la requête est colorée en vert si elle correspond, en
   orange sinon : un enfant repère instantanément sa brique mal placée. */
function Alignment({ data }) {
  if (!data) return null;
  const q = data.query_seq.split('');
  const m = data.midline.split('');
  const t = data.subject_seq.split('');
  return (
    <div className="align">
      <div className="align__head">
        <span>Alignement avec {data.accession}</span>
        <span>
          {data.mismatches > 0
            ? `${data.mismatches} différence${data.mismatches > 1 ? 's' : ''}`
            : 'correspondance parfaite'}
        </span>
      </div>
      {(data.locus || data.subject_length) && (
        <div className="locus">
          <span className="locus__where">
            {data.locus || 'Séquence de référence'}
          </span>
          <span className="locus__pos">
            position {fmt(data.subject_start)}
            {' – '}
            {fmt(data.subject_end)}
            {data.subject_length
              ? ` sur ${fmt(data.subject_length)} bases`
              : ''}
            {data.strand === 'moins' ? ' · brin complémentaire' : ''}
          </span>
        </div>
      )}
      <div className="align__body">
        <div className="align__row">
          <span className="align__label">Ta séquence</span>
          <span className="align__pos">{data.query_start}</span>
          <span className="align__seq">
            {q.map((c, i) => (
              <b key={i} className={m[i] === '|' ? 'ok' : 'ko'}>{c}</b>
            ))}
          </span>
          <span className="align__pos">{data.query_end}</span>
        </div>
        <div className="align__row align__row--mid">
          <span className="align__label" />
          <span className="align__pos" />
          <span className="align__seq">
            {m.map((c, i) => <b key={i}>{c === '|' ? '|' : '\u00a0'}</b>)}
          </span>
        </div>
        <div className="align__row">
          <span className="align__label">ADN de référence</span>
          <span className="align__pos">{data.subject_start}</span>
          <span className="align__seq">
            {t.map((c, i) => (
              <b key={i} className={m[i] === '|' ? 'ok' : 'ko'}>{c}</b>
            ))}
          </span>
          <span className="align__pos">{data.subject_end}</span>
        </div>
      </div>
    </div>
  );
}

/** Une ligne de la fiche. Ne rend rien si la valeur est vide. */
function Fact({ label, value }) {
  if (value === null || value === undefined || value === '') return null;
  return (
    <div className="fact">
      <span className="fact__label">{label}</span>
      <span className="fact__value">{value}</span>
    </div>
  );
}

export default function BioDetective() {
  const [screen, setScreen] = useState('home');
  const [input, setInput] = useState('');
  const [analysed, setAnalysed] = useState('');
  const [organism, setOrganism] = useState(null);
  const [errorMessage, setErrorMessage] = useState('');
  const [apiOk, setApiOk] = useState(null);
  const [pool, setPool] = useState([]);
  const [ready, setReady] = useState([]);
  const [montages, setMontages] = useState([]);
  const [montage, setMontage] = useState(null);
  const [imageIndex, setImageIndex] = useState(0);
  const [message, setMessage] = useState(SEARCH_MESSAGES[0]);
  const [fact, setFact] = useState(DNA_FACTS[0]);
  const [imageLoaded, setImageLoaded] = useState(false);
  const [elapsed, setElapsed] = useState(null);
  const [blast, setBlast] = useState(null);
  const [history, setHistory] = useState([]);
  const [histStats, setHistStats] = useState(null);
  const [fromHistory, setFromHistory] = useState(false);

  // Tous les intervalles vivent ici : un timer oublié continue de tourner
  // en fond et fait clignoter l'écran de résultat.
  const timers = useRef({ poll: null, image: null, message: null });
  const searchStart = useRef(0);

  const clearAllIntervals = useCallback(() => {
    Object.keys(timers.current).forEach((key) => {
      if (timers.current[key]) {
        clearInterval(timers.current[key]);
        timers.current[key] = null;
      }
    });
  }, []);

  // Vérification de l'API et préchargement du pool d'images au démarrage.
  useEffect(() => {
    let alive = true;
    (async () => {
      try {
        // On interroge /stats et non / : le serveur de développement de
        // react-scripts sert son propre index.html sur / et ne le relaie
        // pas au backend. Un fetch('/') répondrait donc 200 avec du HTML
        // même backend éteint — un test de vie qui réussit toujours.
        const r = await fetch(`${API}/stats`, { headers: { Accept: 'application/json' } });
        const data = await r.json();
        if (!r.ok || typeof data.organisms !== 'number') throw new Error('bad payload');
        if (alive) setApiOk(true);
      } catch {
        if (alive) setApiOk(false);
        return;
      }
      // Planches de l'écran de recherche, fabriquées par make_montages.py.
      try {
        const r = await fetch('/montages/index.json');
        const d = await r.json();
        if (alive && d.montages && d.montages.length) {
          setMontages(d.montages);
          setMontage(d.montages[Math.floor(Math.random() * d.montages.length)]);
        }
      } catch {
        /* pas de planches : on tombera sur le repli image par image */
      }
      try {
        const r = await fetch(`${API}/random-images?count=40`);
        const d = await r.json();
        if (alive) setPool(d.images || []);
      } catch {
        /* l'animation tournera sans images, ce n'est pas bloquant */
      }
    })();
    return () => {
      alive = false;
      clearAllIntervals();
    };
  }, [clearAllIntervals]);

  // Précharge le pool dans le cache du navigateur et ne garde que les images
  // effectivement arrivées. Sans ça, la rotation toutes les 300 ms remplace
  // chaque <img> avant que le NCBI n'ait fini de la servir, et le cadre de
  // l'écran de recherche reste noir pendant toute l'analyse.
  useEffect(() => {
    if (!pool.length) return undefined;
    let alive = true;
    pool.forEach((img) => {
      const el = new Image();
      el.onload = () => {
        if (!alive) return;
        setReady((r) => (r.some((x) => x.taxonomy_id === img.taxonomy_id)
          ? r
          : [...r, img]));
      };
      el.src = img.url;
    });
    return () => { alive = false; };
  }, [pool]);

  // Précharge la planche courante : elle pèse ~1 Mo, il ne faut pas la
  // télécharger au moment où l'enfant lance l'analyse.
  useEffect(() => {
    if (!montage) return;
    const el = new Image();
    el.src = montage.file;
  }, [montage]);

  const loadHistory = useCallback(async () => {
    try {
      const r = await fetch(`${API}/history?limit=10&matched_only=true`);
      const d = await r.json();
      setHistory(d.entries || []);
      setHistStats({ total: d.total, found: d.found, unknown: d.unknown });
    } catch {
      /* l'historique est un bonus, jamais un blocage */
    }
  }, []);

  useEffect(() => { loadHistory(); }, [loadHistory]);

  const reset = useCallback(() => {
    clearAllIntervals();
    loadHistory();
    setFromHistory(false);
    if (montages.length) {
      setMontage(montages[Math.floor(Math.random() * montages.length)]);
    }
    setInput('');
    setAnalysed('');
    setOrganism(null);
    setErrorMessage('');
    setImageLoaded(false);
    setScreen('home');
  }, [clearAllIntervals, montages, loadHistory]);

  // Réaffiche une analyse passée sans relancer BLAST : l'enfant qui revient
  // avec ses parents veut revoir son organisme, pas refaire la queue.
  const showPast = useCallback(async (index) => {
    try {
      const r = await fetch(`${API}/history/${index}`);
      if (!r.ok) return;
      const job = await r.json();
      clearAllIntervals();
      setAnalysed(job.sequence || '');
      setElapsed(job.analysis_time || null);
      setBlast(job.blast || null);
      setOrganism(job.organism);
      setImageLoaded(false);
      setFromHistory(true);
      setScreen(job.matched && job.organism ? 'result' : 'unknown');
    } catch {
      /* sans effet : on reste sur l'accueil */
    }
  }, [clearAllIntervals]);

  const fail = useCallback(
    (msg) => {
      clearAllIntervals();
      setErrorMessage(msg);
      setScreen('error');
    },
    [clearAllIntervals]
  );

  const finish = useCallback(
    (job) => {
      clearAllIntervals();
      setAnalysed(job.sequence || '');
      setElapsed(job.analysis_time || null);
      setBlast(job.blast || null);
      if (job.status === 'error') {
        fail(job.error_message || "L'analyse a échoué.");
      } else if (job.matched && job.organism) {
        setOrganism(job.organism);
        setImageLoaded(false);
        setScreen('result');
      } else {
        setFact(DNA_FACTS[Math.floor(Math.random() * DNA_FACTS.length)]);
        setScreen('unknown');
      }
    },
    [clearAllIntervals, fail]
  );

  const startAnalysis = useCallback(async () => {
    if (!input.trim()) return;
    clearAllIntervals();
    setOrganism(null);
    setErrorMessage('');
    setMessage(SEARCH_MESSAGES[0]);
    setScreen('search');

    // Animations de l'écran de recherche. Avec une planche, le défilement
    // est purement CSS et ne coûte aucun timer.
    if (!montage) {
      timers.current.image = setInterval(
        () => setImageIndex((i) => i + 1),
        IMAGE_ROTATE_MS
      );
    }
    searchStart.current = Date.now();
    let m = 0;
    timers.current.message = setInterval(() => {
      const elapsed = Date.now() - searchStart.current;
      // Dernière phase dont le seuil est franchi.
      const phase = SEARCH_PHASES.reduce(
        (acc, p) => (elapsed >= p.after ? p : acc),
        SEARCH_PHASES[0]
      );
      m = (m + 1) % phase.messages.length;
      setMessage(phase.messages[m]);
    }, MESSAGE_ROTATE_MS);

    let jobId;
    try {
      const r = await fetch(`${API}/analyze`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequence: input }),
      });
      const data = await r.json();
      if (!r.ok) {
        // 400 : séquence inutilisable. Le backend renvoie un message
        // déjà rédigé pour être lu par un enfant.
        fail(data.detail || 'Cette séquence ne peut pas être analysée.');
        return;
      }
      jobId = data.job_id;
    } catch {
      fail("Le serveur d'analyse ne répond pas. Vérifie qu'il est démarré.");
      return;
    }

    timers.current.poll = setInterval(async () => {
      try {
        const r = await fetch(`${API}/analyze/${jobId}`);
        if (!r.ok) {
          fail("L'analyse a été perdue. Relance-la.");
          return;
        }
        const job = await r.json();
        if (job.status !== 'running') finish(job);
      } catch {
        fail('La connexion au serveur a été interrompue.');
      }
    }, POLL_MS);
  }, [input, montage, clearAllIntervals, fail, finish]);

  const onKeyDown = (e) => {
    if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) startAnalysis();
  };

  const cleanInput = input.toUpperCase().replace(/[^ATCGU]/g, '').replace(/U/g, 'T');

  // ---------------------------------------------------------------- rendus

  const renderHome = () => (
    <div className="screen screen--home">
      <h1 className="title">BIO<span className="title__accent">DETECTIVE</span></h1>
      <p className="subtitle">Qui se cache derrière cet ADN&nbsp;?</p>

      <div className="panel">
        <label className="panel__label" htmlFor="seq">
          Séquence ADN du Brickopore
        </label>
        <textarea
          id="seq"
          className="seq-input"
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={onKeyDown}
          placeholder="Colle ou tape ta séquence : A T C G…"
          spellCheck="false"
          autoFocus
        />

        <div className="seq-meta">
          <span>{cleanInput.length} base{cleanInput.length > 1 ? 's' : ''}</span>
          <span className="seq-meta__hint">
            {cleanInput.length > 0 && cleanInput.length < MIN_BASES
              ? `Au moins ${MIN_BASES} bases`
              : 'Ctrl + Entrée pour lancer'}
          </span>
        </div>

        <DnaStrip sequence={cleanInput.slice(0, 60)} size="small" />

        <button
          className="btn btn--primary"
          onClick={startAnalysis}
          disabled={cleanInput.length < MIN_BASES || apiOk === false}
        >
          LANCER L'ANALYSE
        </button>
      </div>

      {history.length > 0 && (
        <div className="recent">
          <div className="recent__head">
            <span>Dernières découvertes</span>
            {histStats && (
              <span className="recent__count">
                {histStats.total} analyse{histStats.total > 1 ? 's' : ''}
                {' · '}
                {histStats.found} organisme{histStats.found > 1 ? 's' : ''} identifié
                {histStats.found > 1 ? 's' : ''}
              </span>
            )}
          </div>
          <div className="recent__row">
            {history.map((h) => (
              <button
                key={`${h.index}-${h.at}`}
                className="recent__item"
                onClick={() => showPast(h.index)}
                title={`${h.scientific_name} — ${h.percent_identity} %`}
              >
                <img
                  src={h.image || PLACEHOLDER}
                  alt=""
                  onError={(e) => { e.target.src = PLACEHOLDER; }}
                />
                <span className="recent__name">{h.display_name}</span>
              </button>
            ))}
          </div>
        </div>
      )}

      <div className={`api-status api-status--${apiOk === false ? 'ko' : apiOk ? 'ok' : 'wait'}`}>
        {apiOk === null && 'Connexion au serveur…'}
        {apiOk === true && 'Serveur connecté'}
        {apiOk === false && 'Serveur injoignable — lancer : python3 api.py'}
      </div>
    </div>
  );

  const renderSearch = () => {
    const img = ready.length ? ready[imageIndex % ready.length] : null;
    return (
      <div className="screen screen--search">
        <h2 className="search__title">ANALYSE EN COURS</h2>

        <div className="search__frame">
          {montage ? (
            /* Défilement 100 % CSS : une bande de N vignettes translatée par
               pas entiers. steps(N) sur une translation de -100 % tombe
               exactement sur chaque vignette, puisque la bande fait N fois
               la largeur du cadre. */
            <div className="montage">
              <img
                className="montage__strip"
                src={montage.file}
                alt=""
                style={{
                  animationDuration: `${montage.frames * FRAME_MS}ms`,
                  animationTimingFunction: `steps(${montage.frames})`,
                }}
              />
            </div>
          ) : img ? (
            /* Pas de prop key : elle remonterait un <img> neuf à chaque
               rotation, et un élément fraîchement monté n'a pas encore
               peint. On réutilise le même noeud en ne changeant que src,
               quasi instantané puisque l'image est déjà en cache. */
            <img
              className="search__image"
              src={img.url}
              alt=""
              onError={(e) => { e.target.src = PLACEHOLDER; }}
            />
          ) : (
            <div className="search__waiting">🧬</div>
          )}
          <div className="search__scan" />
        </div>

        <p className="search__message">{message}</p>
        <div className="progress"><div className="progress__bar" /></div>
        <DnaStrip sequence={cleanInput.slice(0, 40)} size="small" />
      </div>
    );
  };

  const renderResult = () => {
    if (!organism) return null;
    const img = organism.image;
    return (
      <div className="screen screen--result">
        <div className="result__grid">
          <div
            className="result__photo"
            tabIndex={0}
            title="Crédit photo au survol"
          >
            {!imageLoaded && <div className="result__loading">Chargement de la photo…</div>}
            <img
              src={img ? img.url : PLACEHOLDER}
              alt={organism.display_name}
              onLoad={() => setImageLoaded(true)}
              onError={(e) => { e.target.src = PLACEHOLDER; setImageLoaded(true); }}
              style={{ opacity: imageLoaded ? 1 : 0 }}
            />
            {img && (img.source || img.attribution) && (
              <p className="credit">
                Photo&nbsp;: {[img.attribution, img.source].filter(Boolean).join(' — ')}
                {img.license ? ` (${img.license.split(' (')[0]})` : ''}
              </p>
            )}
          </div>

          <div className="result__info">
            {fromHistory && (
              <div className="replay">Analyse déjà effectuée — historique</div>
            )}
            <div className="identity">
              IDENTIFIÉ À{' '}
              {blast && blast.alignment
                ? String(blast.alignment.percent_identity).replace('.', ',')
                : '100'}
              &nbsp;%
            </div>
            <h2 className="result__name">{organism.display_name}</h2>
            <p className="result__sci">{organism.scientific_name}</p>
            <div className="badge">{organism.organism_type}</div>

            <div className="facts">
              <Fact label="Règne" value={organism.kingdom} />
              <Fact label="Embranchement" value={organism.phylum} />
              <Fact label="Classe" value={organism.class_name} />
              <Fact label="Ordre" value={organism.order_name} />
              <Fact label="Famille" value={organism.family} />
              <Fact label="Genre" value={organism.genus} />
              <Fact
                label="Taille du génome"
                value={
                  organism.genome_size_mb
                    ? `${fmt(Math.round(organism.genome_size_mb))} millions de bases`
                    : null
                }
              />
              <Fact
                label="Chromosomes"
                value={organism.chromosome_count || null}
              />
              {/* Les champs ci-dessous sont vides pour l'instant : Fact ne
                  rend rien plutôt que d'afficher une ligne « — ». */}
              <Fact label="Taille" value={organism.size_info} />
              <Fact label="Habitat" value={organism.habitat} />
              <Fact label="Le savais-tu ?" value={organism.fun_facts} />
            </div>

            {blast && blast.alignment && (
              <p className="evalue">
                E-value <strong>{blast.alignment.evalue.toExponential(1)}</strong>
                <span className="evalue__hint">
                  {' '}— probabilité que ce soit un hasard
                </span>
              </p>
            )}

            <DnaStrip sequence={analysed} size="small" />
            {elapsed !== null && (
              <p className="result__timing">
                Analyse effectuée en{' '}
                {elapsed.toFixed(1).replace('.', ',')} secondes
              </p>
            )}
            <button className="btn btn--accent" onClick={reset}>
              NOUVELLE ANALYSE
            </button>
          </div>
        </div>

        <Ideogram data={blast && blast.alignment && blast.alignment.ideogram} />
        <Alignment data={blast && blast.alignment} />

        {blast && blast.hits && blast.hits.length > 1 && (
          <div className="others">
            <div className="others__title">
              Autres correspondances trouvées ({blast.total_hits} au total)
            </div>
            {blast.hits.slice(1, 5).map((h) => (
              <div className="others__row" key={h.accession}>
                <span className="others__name">
                  {h.display_name || h.scientific_name || h.accession}
                </span>
                <span className="others__num">{String(h.percent_identity).replace('.', ',')}&nbsp;%</span>
                <span className="others__num">E {h.evalue.toExponential(0)}</span>
              </div>
            ))}
          </div>
        )}
      </div>
    );
  };

  // L'écran le plus vu de la journée : les enfants assemblent les briques
  // librement, donc la plupart des séquences ne sont pas dans la table.
  // Il doit valoriser la découverte, jamais ressembler à une erreur.
  const renderUnknown = () => (
    <div className="screen screen--unknown">
      <div className="unknown__icon">🧬</div>
      <h2 className="unknown__title">SÉQUENCE INCONNUE</h2>
      <p className="unknown__lead">
        Cet ADN ne correspond à aucun organisme de nos bases.
        <br />
        Tu viens peut-être de découvrir une espèce nouvelle&nbsp;!
      </p>

      <div className="panel panel--tight">
        <div className="panel__label">Ta séquence</div>
        <DnaStrip sequence={analysed} />
        <p className="unknown__count">{analysed.length} bases analysées</p>
      </div>

      <p className="unknown__fact">
        <strong>Le savais-tu&nbsp;?</strong> {fact}
      </p>

      <button className="btn btn--primary" onClick={reset}>
        ESSAYER UNE AUTRE SÉQUENCE
      </button>
    </div>
  );

  const renderError = () => (
    <div className="screen screen--error">
      <div className="error__icon">⚠</div>
      <h2 className="error__title">PROBLÈME TECHNIQUE</h2>
      <p className="error__message">{errorMessage}</p>
      <button className="btn btn--accent" onClick={reset}>RETOUR</button>
    </div>
  );

  return (
    <div className="app">
      <div className="app__grid" />
      {screen === 'home' && renderHome()}
      {screen === 'search' && renderSearch()}
      {screen === 'result' && renderResult()}
      {screen === 'unknown' && renderUnknown()}
      {screen === 'error' && renderError()}
      <footer className="app__footer">
        Fête de la Science — IBMP Strasbourg · Brickopore
      </footer>
    </div>
  );
}
