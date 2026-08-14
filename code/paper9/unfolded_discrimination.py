#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# v-bump (revision, Aug 2026): path cascade for zeros_100.csv —
#   script dir, then CWD; no fixed mount path (review C6).
# v-bump (later revision): BASE-relative output defaults + CLI (--cache/--results/
#   --tables), cache version field 'unfolded-v2.1' with validation, --fresh (review N3).
"""
unfolded_discrimination.py — unfolded-scale discrimination battery (AnalysisLab_L1_L5)
AnalysisLab_L1_L5 · Juli 2026 · RH-frei · Pauca sed matura

EIN Skript, ALLE Zahlen der Unfolding-Diskriminationsbatterie. Deterministisch (kein RNG).
Cache defaults to <script-dir>/unfolded_cache.json; override with --cache; --fresh ignores an existing cache.  (--fresh oder Loeschen => voller Frischlauf mit Laufzeiten)

Phasen:
  P1  W_zeta : 31 Ordinaten via mp.zetazero (dps=40) + Abgleich zeros_100.csv
  P2  W_chi4 : L(s,chi4) = 4^{-s}[zeta(s,1/4)-zeta(s,3/4)]; Hardy-Z-Vorzeichenjagd
               + Argumentprinzip-Rechteck (Vollstaendigkeit)
  P3  W_DH   : Davenport-Heilbronn nach EULER-SIG-II §0b (xi aus FG); Ordinaten-only;
               Rechteckzaehlung (findet auch Off-Line-Paare — keine im Fenster erwartet,
               wird BEWIESEN, nicht angenommen)
  P4  Welten / Selbst-Entfaltung b~=N_W(b) / Observablen O1..O5
       (mpmath dps=40 primaer, numpy float64 Kreuzpfad) / vorregistriertes Gate §4
  P5  Explorativ (nicht gate-wirksam): eps~=0.10-Matrix, Robustheits-Vorzeichen,
       Tail-Check W_syn0 bis b~<=60

HARTER FEHLER (implementiert als Struktur): W_DH geht NUR ueber Ordinaten ein;
beta wird nirgends gespeichert oder verwendet. Kein RH/GUE/HP im Code.
"""
import json, os, sys, time
import numpy as np
import os as _os
def _z100():
    _sd = _os.path.dirname(_os.path.abspath(__file__))
    for _p in (_os.path.join(_sd, '..', '..', 'data', 'zeros_100.csv'),
               _os.path.join(_sd, 'data', 'zeros_100.csv'),
               _os.path.join('data', 'zeros_100.csv'),
               _os.path.join(_sd, 'zeros_100.csv'),
               'zeros_100.csv'):
        if _os.path.exists(_p): return _p
    raise SystemExit('zeros_100.csv not found (repo data/, script dir, CWD)')

from mpmath import mp, mpf, mpc

mp.dps = 40
import argparse as _ap
_p=_ap.ArgumentParser(); _BASE=_os.path.dirname(_os.path.abspath(__file__))
_p.add_argument('--cache',   default=_os.path.join(_BASE,'unfolded_cache.json'))
_p.add_argument('--results', default=_os.path.join(_BASE,'unfolded_results.json'))
_p.add_argument('--tables',  default=_os.path.join(_BASE,'unfolded_tables.md'))
_p.add_argument('--fresh', action='store_true', help='ignore existing cache')
_A,_ = _p.parse_known_args()
CACHE_VERSION = 'unfolded-v2.1'
CACHE_FILE = _A.cache
FRESH = _A.fresh
T_START = time.time()
TIMES = {}

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
LOGP = [mp.log(p) for p in PRIMES]
BTMAX = mpf(30)            # Fenster: b~ <= 30 (normativ, §0)
DELTA = mpf('0.25')        # Jitter-Amplitude (normativ, §1)
EPS_GATE = mpf('0.25')     # unfolding-convention default (derived)
EPS_ROB = mpf('0.10')      # Robustheit, nicht gate-tragend

# ---------------- glatte Zaehlfunktionen (RvM-Analoga) ----------------
# Konstanten: zeta 7/8 (normativ, §1); Grad-1 ungerade (a=1), kein Pol => +1/8
# (Herleitung: theta(t)/pi mit theta = (t/2)log(q t/(2 pi e)) + pi/8 + O(1/t));
# empirische Kalibrierung via Residuen r_k = N(gamma_k) - (k - 1/2) im Report.
def N_zeta(T):
    return T/(2*mp.pi)*mp.log(T/(2*mp.pi*mp.e)) + mpf(7)/8

def N_chi4(T):
    return T/(2*mp.pi)*mp.log(4*T/(2*mp.pi*mp.e)) + mpf(1)/8

def N_ded(T):
    return N_zeta(T) + N_chi4(T)

def N_dh(T):
    return T/(2*mp.pi)*mp.log(5*T/(2*mp.pi*mp.e)) + mpf(1)/8

def invN(NF, y, lo):
    """Inversion auf dem monoton steigenden Ast T > lo (Bisektion, 160 Schritte)."""
    lo = mpf(lo); hi = lo + 20
    while NF(hi) < y:
        hi = hi*2
    for _ in range(160):
        m = (lo + hi)/2
        if NF(m) < y:
            lo = m
        else:
            hi = m
    return (lo + hi)/2

# ---------------- L-Funktionen / Hardy-Z ----------------
# EULER-SIG-II §0b: xi^2 + 2*phi*xi - 1 = 0, xi = sqrt((5+sqrt5)/2) - (1+sqrt5)/2
XI = mp.sqrt((5 + mp.sqrt(5))/2) - (1 + mp.sqrt(5))/2

def L_chi4(s):
    return 4**(-s) * (mp.zeta(s, mpf(1)/4) - mp.zeta(s, mpf(3)/4))

def f_dh(s):
    return 5**(-s) * (mp.zeta(s, mpf(1)/5) + XI*mp.zeta(s, mpf(2)/5)
                      - XI*mp.zeta(s, mpf(3)/5) - mp.zeta(s, mpf(4)/5))

def theta_odd(t, q):
    """Hardy-Phase, ungerade Paritaet a=1, Leiter q, Wurzelzahl +1."""
    return t/2*mp.log(mpf(q)/mp.pi) + mp.im(mp.loggamma(mpf(3)/4 + mpc(0, 1)*t/2))

def Z_chi4(t):
    return mp.exp(mpc(0, 1)*theta_odd(t, 4)) * L_chi4(mpc(mpf(1)/2, t))

def Z_dh(t):
    return mp.exp(mpc(0, 1)*theta_odd(t, 5)) * f_dh(mpc(mpf(1)/2, t))

# ---------------- Nullstellenjagd (on-line) ----------------
def find_zeros(Zf, Thi, step=mpf('0.08')):
    """Vorzeichenwechsel von Re Z auf (0, Thi]; Verfeinerung dps=45.
    Rueckgabe: (Nullstellen aufsteigend, schlechteste Realitaet |Im Z|/|Z| auf Gitter)."""
    worst_real = mpf(0)
    brackets = []
    with mp.workdps(20):
        t = mpf('0.05')
        z = Zf(t); v = mp.re(z)
        nsteps = int(mp.ceil((Thi - t)/step)) + 1
        for i in range(1, nsteps + 1):
            t2 = mpf('0.05') + step*i
            z2 = Zf(t2); v2 = mp.re(z2)
            rr = abs(mp.im(z2))/max(abs(z2), mpf('1e-30'))
            if rr > worst_real:
                worst_real = rr
            if v*v2 < 0:
                brackets.append((t, t2))
            t, v = t2, v2
    zeros = []
    for (a, b) in brackets:
        with mp.workdps(20):
            aa, bb = mpf(a), mpf(b)
            fa = mp.re(Zf(aa))
            for _ in range(28):
                mid = (aa + bb)/2
                fm = mp.re(Zf(mid))
                if fa*fm <= 0:
                    bb = mid
                else:
                    aa, fa = mid, fm
        with mp.workdps(45):
            try:
                r = mp.findroot(lambda x: mp.re(Zf(x)), (aa, bb), solver='anderson')
            except Exception:
                fa2 = mp.re(Zf(aa))
                for _ in range(150):
                    mid = (aa + bb)/2
                    fm = mp.re(Zf(mid))
                    if fa2*fm <= 0:
                        bb = mid
                    else:
                        aa, fa2 = mid, fm
                r = (aa + bb)/2
            r = mpf(r)
        zeros.append(r)
    return sorted(zeros), worst_real

# ---------------- Argumentprinzip-Rechteck ----------------
def winding_rect(F, sig1, sig2, t1, t2, dps_c=18):
    """Umlaufzahl von F um Rechteck [sig1,sig2]x[t1,t2], CCW, adaptive Phasenverfolgung.
    Rueckgabe: (W, Anzahl Auswertungen, min|F| auf Kontur)."""
    stats = {'ev': 0, 'minabs': mpf('1e999')}
    sys.setrecursionlimit(300000)
    with mp.workdps(dps_c):
        def FF(z):
            stats['ev'] += 1
            v = F(z)
            a = abs(v)
            if a < stats['minabs']:
                stats['minabs'] = a
            return v
        total = [mpf(0)]
        def seg(z1, z2, f1, f2, depth):
            d = mp.arg(f2/f1)
            if (abs(d) > mpf('0.75') or min(abs(f1), abs(f2)) < mpf('1e-10')) and depth < 42:
                zm = safe((z1 + z2)/2)
                fm = FF(zm)
                seg(z1, zm, f1, fm, depth + 1)
                seg(zm, z2, fm, f2, depth + 1)
            else:
                total[0] += d
        edges = [(mpc(sig2, t1), mpc(sig2, t2), mpf('0.5')),
                 (mpc(sig2, t2), mpc(sig1, t2), mpf('0.05')),
                 (mpc(sig1, t2), mpc(sig1, t1), mpf('0.12')),
                 (mpc(sig1, t1), mpc(sig2, t1), mpf('0.05'))]
        def safe(z):
            # Einzelterme (Hurwitz) haben Pol bei s=1; Kombination regulaer, aber
            # numerisch inf-inf => Knoten minimal verschieben, falls exakt getroffen.
            if abs(z - 1) < mpf('1e-9'):
                return z + mpf('1e-6')
            return z
        for (za, zb, h) in edges:
            Lseg = abs(zb - za)
            nseg = max(4, int(Lseg/h) + 1)
            pts = [safe(za + (zb - za)*mpf(i)/nseg) for i in range(nseg + 1)]
            vals = [FF(p) for p in pts]
            for i in range(nseg):
                seg(pts[i], pts[i + 1], vals[i], vals[i + 1], 0)
        W = total[0]/(2*mp.pi)
    return W, stats['ev'], stats['minabs']

def real_axis_check(F, sig1, sig2, npts=301):
    """Vorzeichenkonstanz + min|F| auf dem reellen Segment [sig1,sig2] (Unterkante).
    Mittelpunkt-Gitter: trifft s=1 (Hurwitz-Polstelle der Einzelterme) nie exakt;
    die Kombinationen sind dort analytisch regulaer (Residuen heben sich)."""
    with mp.workdps(20):
        vals = []
        for i in range(npts):
            sg = sig1 + (sig2 - sig1)*(mpf(i) + mpf(1)/2)/npts
            v = mp.re(F(mpc(sg, 0)))
            assert mp.isfinite(v), 'nicht-finiter Wert bei sigma=%s' % mp.nstr(sg, 10)
            vals.append(v)
        min_abs = min(abs(v) for v in vals)
        signs = set((1 if v > 0 else -1) for v in vals)
    return min_abs, signs

def pick_Trect(zeros_sorted, Tneed):
    """Rechteckhoehe mittig in der Luecke oberhalb Tneed; zaehlt gefundene < T_rect."""
    below = [g for g in zeros_sorted if g <= Tneed]
    above = [g for g in zeros_sorted if g > Tneed]
    assert above, "Jagdbereich zu kurz — Thi erhoehen"
    Trect = (below[-1] + above[0])/2
    return Trect, len(below)

# ---------------- Cache ----------------
cache = {}
if os.path.exists(CACHE_FILE) and not FRESH:
    _c = json.load(open(CACHE_FILE))
    cache = _c if _c.get('_version')==CACHE_VERSION else {}
    cache['_version']=CACHE_VERSION
cache['_version']=CACHE_VERSION

def save_cache():
    json.dump(cache, open(CACHE_FILE, 'w'))

# =========================================================
# P1 — W_zeta
# =========================================================
tp = time.time()
if 'zeta' not in cache:
    zz = [mp.zetazero(k).imag for k in range(1, 32)]
    cache['zeta'] = [mp.nstr(g, 35) for g in zz]
    save_cache()
GZ = [mpf(s) for s in cache['zeta']]
TIMES['P1_zetazero_31'] = time.time() - tp

csv_g = {}
for line in open(_z100()):
    line = line.strip()
    if not line or line[0] not in '0123456789':
        continue
    parts = line.split(',')
    csv_g[int(parts[0])] = mpf(parts[1])
if len(csv_g) < 31:
    raise SystemExit(f'zeros_100.csv must contain at least 31 rows for the '
                     f'P1 cross-check; got {len(csv_g)}')
max_dev_csv = max(abs(GZ[k - 1] - csv_g[k]) for k in range(1, 32))

# =========================================================
# P2 — W_chi4
# =========================================================
tp = time.time()
Tmax_zeta = invN(N_zeta, BTMAX, 2*mp.pi)
Tmax_chi4 = invN(N_chi4, BTMAX, 2)
Tmax_ded = invN(N_ded, BTMAX, 4)
Tmax_dh = invN(N_dh, BTMAX, 2)
Thi_chi4 = max(Tmax_chi4, Tmax_ded) + mpf(2)
if 'chi4' not in cache:
    zc, worst = find_zeros(Z_chi4, Thi_chi4)
    cache['chi4'] = {'zeros': [mp.nstr(g, 30) for g in zc],
                     'worst_real': mp.nstr(worst, 5)}
    save_cache()
GC = [mpf(s) for s in cache['chi4']['zeros']]
TIMES['P2_chi4_zeros'] = time.time() - tp

tp = time.time()
if 'chi4_rect' not in cache:
    Trect_c, nfound_c = pick_Trect(GC, Tmax_chi4)
    mn_c, sg_c = real_axis_check(L_chi4, mpf('-0.5'), mpf('1.5'))
    Wc, ev_c, mina_c = winding_rect(L_chi4, mpf('-0.5'), mpf('1.5'), mpf(0), Trect_c)
    cache['chi4_rect'] = {'Trect': mp.nstr(Trect_c, 20), 'nfound': nfound_c,
                          'W': mp.nstr(Wc, 12), 'ev': ev_c,
                          'minabs': mp.nstr(mina_c, 5),
                          'bottom_minabs': mp.nstr(mn_c, 5),
                          'bottom_signs': sorted(sg_c)}
    save_cache()
TIMES['P2_chi4_rect'] = time.time() - tp

# =========================================================
# P3 — W_DH (Ordinaten-only)
# =========================================================
tp = time.time()
Thi_dh = Tmax_dh + mpf('2.5')
if 'dh' not in cache:
    zd, worstd = find_zeros(Z_dh, Thi_dh)
    cache['dh'] = {'zeros': [mp.nstr(g, 30) for g in zd],
                   'worst_real': mp.nstr(worstd, 5)}
    save_cache()
GD = [mpf(s) for s in cache['dh']['zeros']]
TIMES['P3_dh_zeros'] = time.time() - tp

tp = time.time()
if 'dh_rect' not in cache:
    Trect_d, nfound_d = pick_Trect(GD, Tmax_dh)
    mn_d, sg_d = real_axis_check(f_dh, mpf('-0.5'), mpf('1.5'))
    Wd, ev_d, mina_d = winding_rect(f_dh, mpf('-0.5'), mpf('1.5'), mpf(0), Trect_d)
    cache['dh_rect'] = {'Trect': mp.nstr(Trect_d, 20), 'nfound': nfound_d,
                        'W': mp.nstr(Wd, 12), 'ev': ev_d,
                        'minabs': mp.nstr(mina_d, 5),
                        'bottom_minabs': mp.nstr(mn_d, 5),
                        'bottom_signs': sorted(sg_d)}
    save_cache()
TIMES['P3_dh_rect'] = time.time() - tp

# hochpraezise Realitaets-Sonden (dps=40)
real_probe = {}
with mp.workdps(40):
    for nm, Zf in [('chi4', Z_chi4), ('dh', Z_dh)]:
        wr = mpf(0)
        for tt in [mpf(5), mpf(15), mpf(25), mpf(35)]:
            z = Zf(tt)
            wr = max(wr, abs(mp.im(z))/abs(z))
        real_probe[nm] = wr

# =========================================================
# P4 — Welten, Unfolding, Observablen, Gate
# =========================================================
tp = time.time()
world_raw = {}
world_raw['W_zeta'] = [g for g in GZ if N_zeta(g) <= BTMAX]
world_raw['W_chi4'] = [g for g in GC if N_chi4(g) <= BTMAX]
ded_all = sorted(list(GZ) + list(GC))
world_raw['W_ded'] = [g for g in ded_all if N_ded(g) <= BTMAX]
world_raw['W_dh'] = [g for g in GD if N_dh(g) <= BTMAX]

syn0 = [invN(N_zeta, mpf(k) - mpf(1)/2, 2*mp.pi) for k in range(1, 31)]
world_raw['W_syn0'] = syn0
phi = (1 + mp.sqrt(5))/2
PHIM = [mp.frac(phi**m) for m in range(1, 9)]
for m in range(1, 9):
    lst = []
    for k in range(1, 31):
        s_mk = mp.sin(2*mp.pi*mp.frac(mpf(k)*PHIM[m - 1]))
        lst.append(invN(N_zeta, mpf(k) - mpf(1)/2 + DELTA*s_mk, 2*mp.pi))
    world_raw['W_jit%d' % m] = lst

NW = {'W_zeta': N_zeta, 'W_chi4': N_chi4, 'W_ded': N_ded, 'W_dh': N_dh}
for w in ['W_syn0'] + ['W_jit%d' % m for m in range(1, 9)]:
    NW[w] = N_zeta

bt = {w: [NW[w](g) for g in world_raw[w]] for w in world_raw}

# Rundlauf-Asserts: syn0/jitter exakt auf Ziellattice; Fenster eingehalten
for k in range(30):
    assert abs(bt['W_syn0'][k] - (mpf(k) + mpf(1)/2)) < mpf('1e-30')
for m in range(1, 9):
    for k in range(30):
        tgt = mpf(k) + mpf(1)/2 + DELTA*mp.sin(2*mp.pi*mp.frac(mpf(k + 1)*PHIM[m - 1]))
        assert abs(bt['W_jit%d' % m][k] - tgt) < mpf('1e-30')
for w in bt:
    assert all(x <= BTMAX + mpf('1e-30') for x in bt[w])

def residuals(w):
    return [NW[w](world_raw[w][i]) - (mpf(i) + mpf(1)/2) for i in range(len(world_raw[w]))]

RES = {w: residuals(w) for w in ['W_zeta', 'W_chi4', 'W_ded', 'W_dh']}

def obs_mp(btl, eps):
    n = len(btl)
    w = [mp.exp(-eps**2*b*b/2) for b in btl]
    Phi = mp.matrix(n, 16)
    for i in range(n):
        for j in range(16):
            Phi[i, j] = w[i]*mp.sin(btl[i]*LOGP[j])
    G = Phi.T*Phi
    O3 = mp.fsum(G[j, j] for j in range(16))
    try:
        Ev = mp.eigsy(G, eigvals_only=True)
        lams = sorted([Ev[i] for i in range(16)], reverse=True)
    except TypeError:
        Ev, _ = mp.eigsy(G)
        lams = sorted([Ev[i] for i in range(16)], reverse=True)
    O2 = lams[0]/O3
    O5 = lams[1]/lams[0]
    Tt = Phi*Phi.T
    sa = mp.fsum(abs(Tt[i, j]) for i in range(n) for j in range(n))
    sd = mp.fsum(abs(Tt[i, i]) for i in range(n))
    O4 = (sa - sd)/sa
    O1 = mp.fsum(LOGP[j]**2*mp.fsum(mp.exp(-eps**2*b*b)*b*b*mp.cos(b*LOGP[j])
                                    for b in btl) for j in range(16))
    return [O1, O2, O3, O4, O5]

def obs_np(btl, eps):
    b = np.array([float(x) for x in btl])
    lp = np.array([float(x) for x in LOGP])
    e = float(eps)
    w = np.exp(-e**2*b**2/2)
    Phi = w[:, None]*np.sin(np.outer(b, lp))
    G = Phi.T @ Phi
    O3 = float(np.trace(G))
    lams = np.sort(np.linalg.eigvalsh(G))[::-1]
    O2 = lams[0]/O3
    O5 = lams[1]/lams[0]
    Tt = Phi @ Phi.T
    sa = float(np.abs(Tt).sum())
    sd = float(np.abs(np.diag(Tt)).sum())
    O4 = (sa - sd)/sa
    O1 = float((lp**2*(np.exp(-e**2*b**2)[:, None]*(b**2)[:, None]
                       *np.cos(np.outer(b, lp))).sum(axis=0)).sum())
    return [O1, float(O2), O3, O4, float(O5)]

ALLW = list(bt.keys())
OBS = {}
OBS_NP = {}
for eps_key, eps in [('0.25', EPS_GATE), ('0.10', EPS_ROB)]:
    OBS[eps_key] = {w: obs_mp(bt[w], eps) for w in ALLW}
    OBS_NP[eps_key] = {w: obs_np(bt[w], eps) for w in ALLW}

# Kreuzpfad-Abgleich (Selbsttest; syn0 = Pflicht, hier: alle Welten)
xcheck = {}
for eps_key in OBS:
    worst = 0.0
    for w in ALLW:
        for i in range(5):
            a = float(OBS[eps_key][w][i]); bnp = OBS_NP[eps_key][w][i]
            rel = abs(a - bnp)/max(abs(a), 1e-300)
            worst = max(worst, rel)
    xcheck[eps_key] = worst

GATE_WORLDS = ['W_zeta', 'W_chi4', 'W_ded', 'W_dh']
EULER_WORLDS = ['W_zeta', 'W_chi4', 'W_ded']

def gate_eval(eps_key):
    O = OBS[eps_key]
    S = [max(abs(O['W_jit%d' % m][i] - O['W_syn0'][i]) for m in range(1, 9))
         for i in range(5)]
    D = {w: [abs(O[w][i] - O['W_syn0'][i]) for i in range(5)] for w in GATE_WORLDS}
    SG = {w: [(1 if (O[w][i] - O['W_syn0'][i]) > 0 else -1) for i in range(5)]
          for w in GATE_WORLDS}
    P = {w: [bool(D[w][i] > 3*S[i]) for i in range(5)] for w in GATE_WORLDS}
    return S, D, SG, P

S25, D25, SG25, P25 = gate_eval('0.25')
S10, D10, SG10, P10 = gate_eval('0.10')

# Verdikt nach §4 (Gate bei eps~=0.25)
gate_carriers = []
for i in range(5):
    if P25['W_zeta'][i] and (P25['W_chi4'][i] or P25['W_ded'][i]):
        eulers_pass = [w for w in EULER_WORLDS if P25[w][i]]
        signs = set(SG25[w][i] for w in eulers_pass)
        if len(signs) == 1:
            gate_carriers.append(i)
any_pass_25 = any(P25[w][i] for w in GATE_WORLDS for i in range(5))
if gate_carriers:
    verdict = 'PASS'
elif not any_pass_25:
    verdict = 'FAIL'
else:
    verdict = 'FAIL'   # Zwischenfall: Separation ohne Euler-These (§4) — differenziert berichten

# Robustheit: Vorzeichen bei 0.10 duerfen Verdikt nicht widersprechen.
robust_conflict = []
if verdict == 'PASS':
    for i in gate_carriers:
        for w in EULER_WORLDS:
            if P25[w][i] and SG25[w][i] != SG10[w][i]:
                robust_conflict.append((i + 1, w))
    if robust_conflict:
        verdict = 'INKONKLUSIV'

# Zwischenfall-Diagnose
interm = {w: [i + 1 for i in range(5) if P25[w][i]] for w in GATE_WORLDS}

TIMES['P4_worlds_obs_gate'] = time.time() - tp

# =========================================================
# P5 — Explorativ (markiert, NICHT gate-wirksam)
# =========================================================
tp = time.time()
syn0_60 = [invN(N_zeta, mpf(k) - mpf(1)/2, 2*mp.pi) for k in range(1, 61)]
bt60 = [N_zeta(g) for g in syn0_60]
TAIL = {}
for eps_key, eps in [('0.25', EPS_GATE), ('0.10', EPS_ROB)]:
    o30 = OBS[eps_key]['W_syn0']
    o60 = obs_mp(bt60, eps)
    TAIL[eps_key] = [(o60[i] - o30[i]) for i in range(5)]
TIMES['P5_tail'] = time.time() - tp
TIMES['TOTAL'] = time.time() - T_START

# =========================================================
# Ausgabe: JSON + Markdown-Tabellen
# =========================================================
def ns(x, d=12):
    return mp.nstr(mpf(x), d)

out = {
    'params': {'primes': PRIMES, 'btmax': 30, 'delta': '0.25',
               'eps_gate': '0.25', 'eps_rob': '0.10',
               'xi_dh': mp.nstr(XI, 30),
               'phim': [mp.nstr(x, 25) for x in PHIM]},
    'Tmax': {'W_zeta': ns(Tmax_zeta), 'W_chi4': ns(Tmax_chi4),
             'W_ded': ns(Tmax_ded), 'W_dh': ns(Tmax_dh)},
    'counts': {w: len(bt[w]) for w in ALLW},
    'csv_max_dev_31': mp.nstr(max_dev_csv, 5),
    'realness_scan': {'chi4': cache['chi4']['worst_real'],
                      'dh': cache['dh']['worst_real']},
    'realness_probe_dps40': {k: mp.nstr(v, 5) for k, v in real_probe.items()},
    'rect_chi4': cache['chi4_rect'],
    'rect_dh': cache['dh_rect'],
    'residual_stats': {w: {'mean': mp.nstr(mp.fsum(RES[w])/len(RES[w]), 8),
                           'maxabs': mp.nstr(max(abs(r) for r in RES[w]), 8)}
                       for w in RES},
    'ordinates': {w: [[mp.nstr(world_raw[w][i], 15), mp.nstr(bt[w][i], 15)]
                      for i in range(len(bt[w]))]
                  for w in ['W_zeta', 'W_chi4', 'W_ded', 'W_dh', 'W_syn0']},
    'obs': {ek: {w: [mp.nstr(v, 15) for v in OBS[ek][w]] for w in ALLW} for ek in OBS},
    'xcheck_rel': xcheck,
    'S': {'0.25': [mp.nstr(v, 12) for v in S25], '0.10': [mp.nstr(v, 12) for v in S10]},
    'D': {'0.25': {w: [mp.nstr(v, 12) for v in D25[w]] for w in GATE_WORLDS},
          '0.10': {w: [mp.nstr(v, 12) for v in D10[w]] for w in GATE_WORLDS}},
    'SIGN': {'0.25': SG25, '0.10': SG10},
    'PASS': {'0.25': P25, '0.10': P10},
    'gate_carriers_O': [i + 1 for i in gate_carriers],
    'intermediate_pass_map': interm,
    'robust_conflict': robust_conflict,
    'verdict': verdict,
    'tail_syn0_60_minus_30': {ek: [mp.nstr(v, 8) for v in TAIL[ek]] for ek in TAIL},
    'times_s': {k: round(v, 2) for k, v in TIMES.items()},
}
json.dump(out, open(_A.results, 'w'), indent=1)

# ---- Markdown-Tabellen fuer den Report ----
L = []
A = L.append
ONAME = ['O1 B~[W]', 'O2 lmax/Tr', 'O3 Tr T~', 'O4 OffDiag', 'O5 l2/l1']
A('## Welten-Uebersicht\n')
A('| Welt | T_max(N_W=30) | #Ordinaten im Fenster | Residuen mean | max|r| |')
A('|---|---|---|---|---|')
for w in ['W_zeta', 'W_chi4', 'W_ded', 'W_dh']:
    A('| %s | %s | %d | %s | %s |' % (w, out['Tmax'][w], out['counts'][w],
      out['residual_stats'][w]['mean'], out['residual_stats'][w]['maxabs']))
A('| W_syn0 | — | %d | 0 (exakt) | 0 (exakt) |' % out['counts']['W_syn0'])
A('| W_jit1..8 | — | je %d | — | — |' % out['counts']['W_jit1'])
A('')
for ek in ['0.25', '0.10']:
    A('## Observablen-Matrix eps~ = %s\n' % ek)
    A('| Welt | ' + ' | '.join(ONAME) + ' |')
    A('|---' * 6 + '|')
    for w in ['W_zeta', 'W_chi4', 'W_ded', 'W_dh', 'W_syn0']:
        A('| %s | ' % w + ' | '.join(mp.nstr(mpf(x), 8) for x in out['obs'][ek][w]) + ' |')
    jmin = [min(mpf(out['obs'][ek]['W_jit%d' % m][i]) for m in range(1, 9)) for i in range(5)]
    jmax = [max(mpf(out['obs'][ek]['W_jit%d' % m][i]) for m in range(1, 9)) for i in range(5)]
    A('| W_jit min | ' + ' | '.join(mp.nstr(v, 8) for v in jmin) + ' |')
    A('| W_jit max | ' + ' | '.join(mp.nstr(v, 8) for v in jmax) + ' |')
    A('')
for ek, Sv, Dv, Pv in [('0.25', S25, D25, P25), ('0.10', S10, D10, P10)]:
    A('## D/S-Matrix eps~ = %s  (Zelle: D_i | D_i/S_i | PASS)\n' % ek)
    A('| Welt | ' + ' | '.join(ONAME) + ' |')
    A('|---' * 6 + '|')
    A('| S_i (Streuskala) | ' + ' | '.join(mp.nstr(s, 6) for s in Sv) + ' |')
    A('| 3*S_i (Schwelle) | ' + ' | '.join(mp.nstr(3*s, 6) for s in Sv) + ' |')
    for w in GATE_WORLDS:
        cells = []
        for i in range(5):
            ratio = Dv[w][i]/Sv[i] if Sv[i] > 0 else mpf('inf')
            cells.append('%s | %s | %s' % (mp.nstr(Dv[w][i], 5), mp.nstr(ratio, 4),
                                           'PASS' if Pv[w][i] else 'fail'))
        A('| %s | ' % w + ' || '.join(cells) + ' |')
    A('')
A('## Ordinaten (roh -> entfaltet), 15 Stellen\n')
for w in ['W_zeta', 'W_chi4', 'W_ded', 'W_dh', 'W_syn0']:
    A('### %s (n = %d)\n' % (w, out['counts'][w]))
    A('| k | b_k (roh) | b~_k (entfaltet) |')
    A('|---|---|---|')
    for i, (raw, unf) in enumerate(out['ordinates'][w]):
        A('| %d | %s | %s |' % (i + 1, raw, unf))
    A('')
open(_A.tables, 'w').write('\n'.join(L))

print('=== Unfolded discrimination battery — Ergebnis ===')
print('Verdikt nach §4:', verdict)
print('Gate-tragende Observablen:', out['gate_carriers_O'])
print('PASS-Landkarte (0.25):', interm)
print('Kreuzpfad worst rel:', xcheck)
print('CSV-Abgleich max|dev| (k=1..31):', out['csv_max_dev_31'])
print('Rect chi4: W =', cache['chi4_rect']['W'], ' gefunden <Trect:', cache['chi4_rect']['nfound'])
print('Rect DH  : W =', cache['dh_rect']['W'], ' gefunden <Trect:', cache['dh_rect']['nfound'])
print('Realness (dps40-Sonden):', out['realness_probe_dps40'])
print('Residuen:', out['residual_stats'])
print('Counts:', out['counts'])
print('Tmax:', out['Tmax'])
print('Tail (syn0 60-30):', out['tail_syn0_60_minus_30'])
print('Zeiten [s]:', out['times_s'])
