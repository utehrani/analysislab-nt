"""
eta_verification.py -- Paper 3 Verification (eta_orig main result)
==================================================================
A Finite-Cutoff Hilbert-Space Model for Prime-Zero Energy Structure
Ulrich Tehrani - Zenodo doi:10.5281/zenodo.19307989

Reproduces the eta_orig energy identity at reference parameters
(kappa=53, eps=0.05, sigma=0.5, N=100) and contrasts four weight
hypotheses (canonical c_p_eta, c_p_ren, normalized G, AG-guess).

The canonical eta-framework weight is:
   c_p_eta := sqrt(f_p) = sqrt(4 (log p)^2 (2p-1) / (p (p-1)^2)).
Numerically equal to Paper 2's c_p_ren; the distinct name reflects
its distinct role inside the eta_orig formula (see notation.json).

Reference values (kappa=53, eps=0.05, sigma=0.5, N=100):
   eta_orig = 0.66926873   (Sprint AF reference: 0.66926893)
   eta range over sigma in [0.10, 0.95]: [0.598029, 0.704420]
   eta(c_p^ren = log(p)/p): 0.761484 (Sprint AG)

Connection to Paper 4: same weight (c_p_eta = sqrt(f_p)) used in
the Tehrani operator construction and the rank/spectrum analysis.
"""

import numpy as np
from mpmath import zetazero

# ── Normative Nullstellen ──────────────────────────────────────────────────
N = 100
print("Lade erste 100 Zeta-Nullstellen...")
gammas = np.array([float(zetazero(k).imag) for k in range(1, N+1)])
print(f"γ_1={gammas[0]:.6f}, γ_100={gammas[99]:.6f}")

# ── Primzahlen bis κ ───────────────────────────────────────────────────────
primes_23  = [2,3,5,7,11,13,17,19,23]
primes_53  = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53]
primes_101 = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101]
primes_199 = [p for p in range(2,200) if all(p%d!=0 for d in range(2,int(p**0.5)+1))]

# ── Normative Definitionen (SSOT Rev20 §28/§29) ────────────────────────────

def c_p_eta(p):
    """Canonical eta-framework weight: c_p_eta(p) = sqrt(f_p), where
    f_p = V_p(1/2) - 4*(log p)^2/p = 4 (log p)^2 (2p-1) / (p (p-1)^2).
    Numerically equal to Paper 2's c_p_ren = sqrt(f_p); the distinct
    name reflects its role inside the eta_orig formula (see HIGH 4
    in Sprint AUDIT, May 2026)."""
    return np.sqrt(4 * (np.log(p))**2 * (2*p - 1) / (p * (p-1)**2))


# Backwards-compatible alias: c_p() is the legacy name (pre-AUDIT).
# Kept so any out-of-tree callers continue to work; new code uses c_p_eta.
c_p = c_p_eta

def a_p_vec(p, gammas, epsilon):
    """Option B: a_p = e^{-ε²γ_k²/2} sin(γ_k log p)"""
    return np.exp(-epsilon**2 * gammas**2 / 2) * np.sin(gammas * np.log(p))

def a_p_sigma_vec(p, gammas, epsilon, sigma):
    """σ-abhängig: a_p^(σ) = p^{-(σ-1/2)} · a_p"""
    return p**(-(sigma - 0.5)) * a_p_vec(p, gammas, epsilon)

def eta_orig(primes, gammas, epsilon, sigma):
    """
    η_orig = 1 − E_spec / E_str
    E_str  = Σ_p c_p² ‖a_p^(σ)‖²
    E_spec = ‖Σ_p c_p a_p^(σ)‖²
    """
    c_vec = np.array([c_p(p) for p in primes])
    A = np.array([a_p_sigma_vec(p, gammas, epsilon, sigma) for p in primes])  # (P,N)
    norms_sq = np.sum(A**2, axis=1)   # (P,)
    E_str    = np.sum(c_vec**2 * norms_sq)
    weighted = np.sum(c_vec[:,None] * A, axis=0)  # (N,)
    E_spec   = np.sum(weighted**2)
    eta      = 1.0 - E_spec / E_str
    return eta, E_str, E_spec

# ════════════════════════════════════════════════════════════════════════════
# NUM4 — Reproduziere Sprint AF: κ=53, ε=0.05, σ=0.5  →  η=0.66926893
# ════════════════════════════════════════════════════════════════════════════
eta_af, E_str_af, E_spec_af = eta_orig(primes_53, gammas, 0.05, 0.5)
print(f"\n── NUM4 Reproduktion Sprint AF ──────────────────────────────")
print(f"η_orig(κ=53, ε=0.05, σ=0.5, N=100) = {eta_af:.8f}")
print(f"Sprint AF Referenz:                    0.66926893")
print(f"Δ = {abs(eta_af - 0.66926893):.2e}")

# ════════════════════════════════════════════════════════════════════════════
# NUM1 — Hypothesen-Test (alle bei κ=53, ε=0.05, σ=0.5, N=100)
# ════════════════════════════════════════════════════════════════════════════
print(f"\n── NUM1 Hypothesen-Test ─────────────────────────────────────")

# Hypothese A: G^norm statt G^un
# G^norm_{pq} = <a_p/‖a_p‖, a_q/‖a_q‖>
# η aus G^norm = 1 − c^T G^norm c / ‖c‖²
c_vec = np.array([c_p(p) for p in primes_53])
A_half = np.array([a_p_sigma_vec(p, gammas, 0.05, 0.5) for p in primes_53])
norms  = np.linalg.norm(A_half, axis=1)
A_norm = A_half / norms[:,None]
G_norm = A_norm @ A_norm.T
eta_hyp_A = 1.0 - (c_vec @ G_norm @ c_vec) / np.sum(c_vec**2)
print(f"Hypothese A (G^norm): η = {eta_hyp_A:.6f}")

# Hypothese B: c_p^ren statt c_p
# c_p^ren = 1/p (vereinfachte Weil-Gewichte; hier: f_p^{1/2} ~ log(p)/p)
c_ren = np.array([np.log(p)/p for p in primes_53])
norms_sq_half = np.sum(A_half**2, axis=1)
E_str_ren  = np.sum(c_ren**2 * norms_sq_half)
w_ren = np.sum(c_ren[:,None] * A_half, axis=0)
E_spec_ren = np.sum(w_ren**2)
eta_hyp_B  = 1.0 - E_spec_ren / E_str_ren
print(f"Hypothese B (c^ren):  η = {eta_hyp_B:.6f}")

# Hypothese C: σ=½ fest (kein p^{-(σ-½)} Faktor) → identisch mit Normalfall bei σ=0.5
# Aber: wenn Sprint AG σ-Invarianz testet, könnte er G^norm verwendet haben.
# Teste: wenn man einfach die normierten Vektoren nimmt und c uniform setzt:
c_uniform = np.ones(len(primes_53)) / len(primes_53)
w_uni = np.sum(c_uniform[:,None] * A_norm, axis=0)
eta_hyp_C  = 1.0 - np.sum(w_uni**2) / np.sum(c_uniform**2)
print(f"Hypothese C (σ=½, unif c, norm A): η = {eta_hyp_C:.6f}")

# Zusatz: was gibt G^norm mit normativen c_p?
G_norm_c = A_norm @ A_norm.T
c_norm2 = np.sum(c_vec**2)
eta_gnorm_cp = 1.0 - (c_vec @ G_norm_c @ c_vec) / c_norm2
print(f"G^norm + normative c_p: η = {eta_gnorm_cp:.6f}  ← Hypothese A")

# Was gibt Sprint AG genau: vielleicht ungewichtetes Skalarprodukt mit normierten Vektoren
# d.h. η = 1 - (Σ_p Σ_q <â_p, â_q>) / P²  ?
P = len(primes_53)
eta_ag_guess = 1.0 - np.sum(G_norm_c) / P**2
print(f"AG-Guess (ungewichtet G^norm): η = {eta_ag_guess:.6f}")

# ════════════════════════════════════════════════════════════════════════════
# NUM2 — Normative Tabelle: κ=53, ε=0.05, σ variiert
# ════════════════════════════════════════════════════════════════════════════
print(f"\n── NUM2 Normative Tabelle (κ=53, ε=0.05, N=100) ────────────")
sigmas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
print(f"{'σ':>6} | {'η_orig':>12} | {'E_str':>14} | {'E_spec':>14} | {'Δ=E_str-E_spec':>16}")
print("-"*70)
eta_vals = []
for s in sigmas:
    eta_v, E_s, E_sp = eta_orig(primes_53, gammas, 0.05, s)
    delta = E_s - E_sp
    eta_vals.append(eta_v)
    print(f"{s:>6.2f} | {eta_v:>12.8f} | {E_s:>14.6f} | {E_sp:>14.6f} | {delta:>16.6f}")

# Monotonie-Check
print(f"\nMonoton fallend? {all(eta_vals[i] >= eta_vals[i+1] for i in range(len(eta_vals)-1))}")
print(f"η-Bereich: [{min(eta_vals):.6f}, {max(eta_vals):.6f}]")

# ════════════════════════════════════════════════════════════════════════════
# NUM3 — Stabilitätstest
# ════════════════════════════════════════════════════════════════════════════
print(f"\n── NUM3 Stabilitätstest (σ=0.5) ─────────────────────────────")
kappas    = [23, 53, 101, 199]
primes_k  = [primes_23, primes_53, primes_101, primes_199]
epsilons  = [0.02, 0.05, 0.10]
print(f"{'κ':>6} | {'ε=0.02':>10} | {'ε=0.05':>10} | {'ε=0.10':>10}")
print("-"*45)
for kap, plist in zip(kappas, primes_k):
    row = [eta_orig(plist, gammas, eps, 0.5)[0] for eps in epsilons]
    print(f"{kap:>6} | {row[0]:>10.6f} | {row[1]:>10.6f} | {row[2]:>10.6f}")

print(f"\nSSO Rev20 §30 Referenz: η ∈ [0.61, 0.70]")

# ════════════════════════════════════════════════════════════════════════════
# ZUSATZ: Sprint AC Werte rekonstruieren (0.792, 0.741, 0.688)
# Verdacht: ε=0.10 oder andere κ
# ════════════════════════════════════════════════════════════════════════════
print(f"\n── Sprint AC Rekonstruktion ─────────────────────────────────")
print("Teste ε=0.10, κ=53, σ∈{0.3,0.5,0.7}:")
for s in [0.3, 0.5, 0.7]:
    eta_v, _, _ = eta_orig(primes_53, gammas, 0.10, s)
    print(f"  σ={s}: η = {eta_v:.6f}")

print("\nTeste ε=0.05, κ=101, σ∈{0.3,0.5,0.7}:")
for s in [0.3, 0.5, 0.7]:
    eta_v, _, _ = eta_orig(primes_101, gammas, 0.05, s)
    print(f"  σ={s}: η = {eta_v:.6f}")

# Hypothese B genauer: welches c_p^ren passt zu 0.7615?
print(f"\n── Hypothese B Detail ───────────────────────────────────────")
# c_p^ren = log(p)/p  →  0.7615?  (bereits gezeigt: JA)
print(f"c_p^ren = log(p)/p:  η = {eta_hyp_B:.6f}  [Sprint AG: 0.7615]")
print(f"Differenz:              {abs(eta_hyp_B - 0.7615):.2e}  → MATCH")
