#!/usr/bin/env python3
# =============================================================================
# Function-side reconstruction — konsolidiertes Reproduktionsskript
# Paper 9 supporting artefact · Juli 2026 · Status: NUMERISCH
#
# PROVENIENZ:
#   Stage main     : Ursprungsfassung (4 Modelle) + spätere Ergänzung
#                    (A1P, B\6 — im Report tabelliert, im Ur-Skript fehlend)
#                    + Chance-Baseline (aus dem Review ergänzt)
#   Stage ablation : Ablations-Nachlauf [650,990] (im Ur-Artefakt fehlend;
#                    unabhängig nacherzeugt) + Zusatzkontrollen n=8, n=9
#                    (voller Support; n=12 ist support-konfundiert: 25 %)
#
# SOLLWERTE (Selbstcheck; unabhängig verifiziert 11.07.2026):
#   main:     B 98.1/0.0/0.980 · A_RS 34.6/35.4/0.534 · A_53 34.6/31.4/0.675
#             A1 79.8/5.7/0.871 · A1P = B\6 = 85.9/2.2/0.930 · Chance ≈ 22–39 %
#   ablation: Voll{1..12} 98.5 % · Δpp: 6:−14.8 · 10:−9.1 · 9:−8.0 · 8:−7.2
#             7:−3.0 · 12:−0.8 · {6,10}:−20.9
#
# KORREKTUR-HINWEIS: „Schadensranking folgt μ(n)" ist in
# der starken Form FALSIFIZIERT (n=8/9, μ=0, liegen im Semiprim-Band).
# Operativer Träger: Kompositheit (Kombinationstöne aktiver Primfrequenzen).
# 6er-Anomalie (−18.2/Amp, +26 % vor n=10) OFFEN.
#
# Aufruf: python3 function_side_reconstruction.py [main|ablation|all]
# Nicht-Zirkularität: kein RH/GUE/HP-Input; θ(t) Stirling; Nullstellen nur
# als Wahrheitsseite (Intervallgrenzen/Matching-Ziele).
# =============================================================================
import sys, time, csv
import numpy as np
import mpmath as mp
from scipy import stats

H = 0.04
PRIMES = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53]
import os as _os
def _z650():
    _sd = _os.path.dirname(_os.path.abspath(__file__))
    for _p in (_os.path.join(_sd, '..', '..', 'data', 'zeros_650.csv'),
               _os.path.join(_sd, 'data', 'zeros_650.csv'),
               _os.path.join('data', 'zeros_650.csv'),
               _os.path.join(_sd, 'zeros_650.csv'),
               'zeros_650.csv'):
        if _os.path.exists(_p): return _p
    raise SystemExit('zeros_650.csv not found (repo data/, script dir, CWD)')
ZPATH = _z650()
try: open(ZPATH)
except OSError: ZPATH = 'zeros_650.csv'
gam = []
with open(ZPATH) as f:
    r = csv.reader(f); next(r)
    for row in r: gam.append(float(row[1]))
gam = np.array(gam)
EXPECTED_ROWS = 650           # named methods artefact: the extended ordinate list
if len(gam) != EXPECTED_ROWS:
    raise SystemExit(f'zeros_650.csv must contain exactly {EXPECTED_ROWS} rows '
                     f'for this experiment; got {len(gam)}')
mp.mp.dps = 8

def theta_of(t):
    return t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)

def build(T0, T1):
    t = np.arange(T0, T1 + H, H)
    th = theta_of(t)
    t0 = time.time()
    Z = np.array([float(mp.siegelz(x)) for x in t])
    print(f"Z-Gitter [{T0:.0f},{T1:.0f}]: {len(t)} Pkte, {time.time()-t0:.0f}s")
    zs = gam[(gam > T0) & (gam < T1)]
    return t, th, Z, zs

def field(t, th, ns, windowed=True):
    F = np.zeros_like(t); cut = np.sqrt(t/(2*np.pi))
    for n in ns:
        term = 2*n**-0.5*np.cos(th - t*np.log(n)) if n > 1 else 2*np.cos(th)
        F += np.where(n <= cut, term, 0.0) if windowed else term
    return F

def zeros_of(t, F):
    s = np.sign(F); i = np.where(s[:-1]*s[1:] < 0)[0]
    return t[i] - F[i]*H/(F[i+1] - F[i])

def spac(x): return 2*np.pi/np.log(x/(2*np.pi))

def recovery(zs, mz):
    tol = 0.25*spac(zs)
    rec = np.mean(np.min(np.abs(mz[None,:] - zs[:,None]), axis=1) < tol)
    tolm = 0.25*spac(mz)
    spur = np.mean(np.min(np.abs(zs[None,:] - mz[:,None]), axis=1) > tolm)
    return rec, spur, len(mz)

def stage_main():
    T0, T1 = 100.0, 600.0
    t, th, Z, zs = build(T0, T1)
    print(f"wahre Nullstellen: {len(zs)}")
    def areas(F):
        return np.array([np.trapezoid(np.abs(F[(t>=a)&(t<=b)]), t[(t>=a)&(t<=b)])
                         for a, b in zip(zs[:-1], zs[1:])])
    S_true = areas(Z)
    models = {
        'B (RS {1..9})':   field(t, th, range(1,10)),
        'A_RS (Prim,RS)':  field(t, th, PRIMES),
        'A_53 (16P,ohne)': field(t, th, PRIMES, windowed=False),
        'A1 ({1}+P,RS)':   field(t, th, [1]+PRIMES),
        'A1P ({1,4,8,9}+P)': field(t, th, [1,4,8,9]+PRIMES),
        'B\\6 ({1..9}\\6)': field(t, th, [n for n in range(1,10) if n != 6]),
    }
    rng = np.random.default_rng(0)
    print(f"{'Modell':<18}{'Recov':>8}{'Spur':>8}{'#Z':>6}{'Spear':>8}{'Pear':>8}{'Chance':>8}")
    for nm, F in models.items():
        mz = zeros_of(t, F); rec, spur, nz = recovery(zs, mz)
        rho = stats.spearmanr(S_true, areas(F)).statistic
        pear = np.corrcoef(F, Z)[0,1]
        ch = [np.mean(np.min(np.abs(np.sort(rng.uniform(T0,T1,nz))[None,:]
              - zs[:,None]), axis=1) < 0.25*spac(zs)) for _ in range(200)]
        print(f"{nm:<18}{rec*100:>7.1f}%{spur*100:>7.1f}%{nz:>6d}{rho:>8.3f}{pear:>8.3f}{np.mean(ch)*100:>7.1f}%")

def stage_ablation():
    T0, T1 = 650.0, 990.0
    t, th, Z, zs = build(T0, T1)
    cut = np.sqrt(t/(2*np.pi))
    print(f"wahre Nullstellen: {len(zs)} · Support: " +
          ", ".join(f"n={n}:{np.mean(n<=cut)*100:.0f}%" for n in [10,11,12]))
    FULL = list(range(1,13))
    def rec_of(ns):
        mz = zeros_of(t, field(t, th, ns))
        return np.mean(np.min(np.abs(mz[None,:]-zs[:,None]),axis=1) < 0.25*spac(zs))
    base = rec_of(FULL)
    print(f"Voll {{1..12}}: {base*100:.1f}%")
    print(f"{'entfernt':<10}{'Rec':>8}{'Δpp':>8}{'Amp':>8}{'Δ/Amp':>8}")
    for rem in [[6],[10],[9],[8],[7],[12],[6,10]]:
        r = rec_of([n for n in FULL if n not in rem])
        d = (r-base)*100; amp = sum(2*n**-0.5 for n in rem)
        print(f"n={str(rem):<8}{r*100:>7.1f}%{d:>+8.1f}{amp:>8.3f}{d/amp:>8.1f}")

if __name__ == '__main__':
    todo = sys.argv[1:] or ['all']
    for st in (['main','ablation'] if todo == ['all'] else todo):
        print("\n" + "="*70 + f"\nSTAGE {st}\n" + "="*70)
        globals()['stage_'+st]()
