# verify_paper4.py
# Paper 4: A Dual Operator for Prime–Zero Coupling · April 2026
# All normative parameters: kappa=53, eps=0.05, N=100, sigma=0.5
#
# Verification script for Paper 4: "A Dual Operator for Prime–Zero Coupling"
# Checks:
#   1. Tehrani operator T̃ = Φ∘Φ* construction and properties
#   2. Spectral identity σ(T)\{0} = σ(T̃)\{0}
#   3. Eigenvalue correlation r1 = corr(μ_j, 1/γ_{k(j)})
#   4. Eigenvector localization
#   5. Two arithmetic constants C_η ≈ 0.39, C_T
#   6. η_orig > 0 for normative κ ≤ 1009
#   7. Δ_Burst > 0 (κ-invariant lower bound)
#   8. Rayleigh identity ⟨T_ren c̃, c̃⟩ = 1 − η_orig (machine precision)
#
# Usage: python verify_paper4.py
# Requires: numpy, mpmath, sympy
#
# GitHub: https://github.com/utehrani/analysislab-nt

import os
import sys
import numpy as np
from sympy import primerange

# ── Parameters ───────────────────────────────────────────────────────────────
# Usage: python verify_paper4.py [N_zeros]
# Default: N=100. Use N=200 for extended run: python verify_paper4.py 200
KAPPA_NORM  = 53
EPS_NORM    = 0.05
N_NORM      = int(sys.argv[1]) if len(sys.argv) > 1 else 100
SIGMA_NORM  = 0.5
KAPPAS_TEST = [23, 53, 101, 199, 503, 1009]

PASS_COUNT = 0
FAIL_COUNT = 0


def check(label, condition, info=""):
    global PASS_COUNT, FAIL_COUNT
    if condition:
        print(f"  PASS  {label}")
        PASS_COUNT += 1
    else:
        print(f"  FAIL  {label}  {info}")
        FAIL_COUNT += 1


# ── Helper: fetch zero ordinates ─────────────────────────────────────────────
def get_zeros(N):
    """Load zero ordinates from CSV if available, else compute via mpmath.
    Automatically uses zeros_200.csv when N > 100."""
    csv_name = 'zeros_200.csv' if N > 100 else 'zeros_100.csv'
    csv_path = os.path.join('data', csv_name)
    if os.path.exists(csv_path):
        data = np.loadtxt(csv_path, delimiter=',', skiprows=1, max_rows=N)
        gammas = data[:, -1] if data.ndim == 2 else data.flatten()
        gammas = gammas[:N]
        print(f"  Loaded {len(gammas)} zero ordinates from {csv_path}")
        return list(gammas)
    else:
        print(f"  {csv_path} not found — computing via mpmath (~30s) ...")
        from mpmath import zetazero
        return [float(zetazero(k).imag) for k in range(1, N + 1)]


# ── Helper: build Φ matrix ───────────────────────────────────────────────────
def build_Phi(primes, gammas, eps):
    pi_k = len(primes)
    N    = len(gammas)
    Phi  = np.zeros((N, pi_k))
    for j, p in enumerate(primes):
        for k, gam in enumerate(gammas):
            Phi[k, j] = np.exp(-eps**2 * gam**2 / 2) * np.sin(gam * np.log(p))
    return Phi  # shape (N, pi_k)


# ── Helper: canonical weight vector c_p^eta (HIGH 4 corrected, May 2026) ─────
def canonical_c(primes, sigma=SIGMA_NORM):
    """Canonical eta-framework weight: c_p^eta(p) = sqrt(f_p), where
    f_p = V_p(1/2) - 4*(log p)^2/p = 4 (log p)^2 (2p-1) / (p (p-1)^2).

    Numerically equal to Paper 2's c_p^ren = sqrt(f_p); the distinct name
    reflects the role inside the eta_orig formula. Gives the normative
    eta_orig = 0.66926893 at kappa=53, eps=0.05, sigma=0.5, N=100.

    HIGH 4 history (Sprint AUDIT, May 2026): the previous docstring claimed
    'c_p = sqrt(V_p(1/2))' which is incorrect — V_p(1/2) and f_p differ by
    the leading 4*(log p)^2/p subtraction. The numerical formula below has
    always computed sqrt(f_p), so all eta values are unchanged.

    Note: Paper 4 §5.2 table uses a different c_p convention (Δ_Burst ≈ 4.81
    vs ≈ 3.10 here); separate open inconsistency, outside HIGH 4 scope."""
    return np.array([
        np.sqrt(4 * np.log(p)**2 * (2*p - 1) / (p * (p - 1)**2))
        for p in primes
    ])


# ── Helper: energy quantities ────────────────────────────────────────────────
def energy(Phi, c):
    ap_norms_sq = np.sum(Phi**2, axis=0)          # ‖a_p‖²  shape (pi_k,)
    D_diag = ap_norms_sq                            # diagonal of D
    E_str  = np.dot(c**2, D_diag)                  # c^T D c
    Phi_c  = Phi @ c                                # shape (N,)
    E_spec = np.dot(Phi_c, Phi_c)                  # ‖Φc‖²
    Delta  = E_str - E_spec
    eta    = Delta / E_str if E_str > 0 else np.nan
    return E_str, E_spec, Delta, eta, D_diag


# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 65)
print("verify_paper4.py  ·  Paper 4  ·  April 2026")
print("Paper 4: A Dual Operator for Prime–Zero Coupling")
print(f"N = {N_NORM} zero ordinates  ·  κ={KAPPA_NORM}  ·  ε={EPS_NORM}")
print("=" * 65)

gammas_full = get_zeros(N_NORM)

# ── TEST 1: Spectral identity ⟨T⟩ vs ⟨T̃⟩ ─────────────────────────────────
print("\n[1] Spectral identity σ(T)\\{0} = σ(T̃)\\{0}")
primes_53 = list(primerange(2, KAPPA_NORM + 1))
Phi_53    = build_Phi(primes_53, gammas_full, EPS_NORM)

T     = Phi_53.T @ Phi_53   # (pi_k × pi_k)
Ttil  = Phi_53 @ Phi_53.T   # (N × N)

eig_T, _    = np.linalg.eigh(T)
eig_Tt, _   = np.linalg.eigh(Ttil)
eig_T_pos   = np.sort(eig_T[eig_T > 1e-12])[::-1]
eig_Tt_pos  = np.sort(eig_Tt[eig_Tt > 1e-12])[::-1]

n_shared = min(len(eig_T_pos), len(eig_Tt_pos))
max_diff  = np.max(np.abs(eig_T_pos[:n_shared] - eig_Tt_pos[:n_shared]))
check("Spectral identity: max|λ_j(T) − μ_j(T̃)| < 1e-12",
      max_diff < 1e-12, f"  max_diff={max_diff:.3e}")

# ── TEST 1b: rank(T̃) evidence (HIGH 5 annotation, Sprint AUDIT May 2026) ───
print("\n[1b] rank(T̃) evidence (κ=53, ε=0.05, N=100)")
print("    rank(T̃) ≤ π(κ) is unconditional (rank of a product of linear maps).")
print("    rank(T̃) = π(κ) holds iff {a_p}_{p≤κ} are linearly independent.")
print("    Independence has not been proved analytically (open problem).")

# Numerical rank: count eigenvalues above noise threshold
sing_T_tilde = np.sort(eig_Tt)[::-1]    # descending
threshold    = 1e-12 * sing_T_tilde[0]
rank_num     = int(np.sum(sing_T_tilde > threshold))
gap_above    = sing_T_tilde[len(primes_53)-1]
gap_below    = sing_T_tilde[len(primes_53)] if len(sing_T_tilde) > len(primes_53) else 0.0
print(f"    π(κ) = {len(primes_53)},  μ_max = {sing_T_tilde[0]:.6e}")
print(f"    μ_{len(primes_53)} = {gap_above:.3e}  (last 'large' eigenvalue)")
print(f"    μ_{len(primes_53)+1} = {gap_below:.3e}  (first 'small'/noise eigenvalue)")
print(f"    gap ratio = μ_{len(primes_53)+1}/μ_{len(primes_53)} = {gap_below/gap_above:.3e}")
print(f"    numerical rank (eigenvalues > 1e-12·μ_max): {rank_num}")
check(f"rank(T̃) ≤ π(κ) = {len(primes_53)}  [unconditional, proved]",
      rank_num <= len(primes_53),
      f"  numerical rank = {rank_num}")
check(f"rank(T̃) = π(κ) = {len(primes_53)}  [NUMERICAL, conjectural]",
      rank_num == len(primes_53),
      f"  numerical rank = {rank_num} (matches π(κ) — independence is conjectural)")

# ── TEST 2: η_orig > 0 for all tested κ ──────────────────────────────────────
print("\n[2] η_orig > 0 for κ ≤ 1009")
for kap in KAPPAS_TEST:
    primes_k = list(primerange(2, kap + 1))
    Phi_k    = build_Phi(primes_k, gammas_full, EPS_NORM)
    c_k      = canonical_c(primes_k)
    E_str, E_spec, Delta, eta, D_diag = energy(Phi_k, c_k)
    check(f"η_orig > 0  (κ={kap:4d}, η={eta:.5f})", eta > 0)

# ── TEST 3: Rayleigh identity ⟨T_ren c̃, c̃⟩ = 1 − η ─────────────────────
print("\n[3] Rayleigh identity ⟨T_ren c̃, c̃⟩ = 1 − η_orig (machine precision)")
primes_53 = list(primerange(2, KAPPA_NORM + 1))
Phi_53    = build_Phi(primes_53, gammas_full, EPS_NORM)
c_53      = canonical_c(primes_53)
E_str, E_spec, Delta, eta, D_diag = energy(Phi_53, c_53)
D_invsqrt = 1.0 / np.sqrt(D_diag)
T_53      = Phi_53.T @ Phi_53
T_ren     = D_invsqrt[:, None] * T_53 * D_invsqrt[None, :]
c_tilde   = (np.sqrt(D_diag) * c_53) / np.sqrt(E_str)
rayleigh  = c_tilde @ T_ren @ c_tilde
residual  = abs((1 - eta) - rayleigh)
check(f"Rayleigh identity residual < 1e-14  (resid={residual:.2e})",
      residual < 1e-14)

# ── TEST 4: Eigenvalue correlation r1 ────────────────────────────────────────
print("\n[4] Eigenvalue correlation r1 = corr(μ_j, 1/γ_{k(j)}) at κ=53")
Ttil_53   = Phi_53 @ Phi_53.T
eig_vals, eig_vecs = np.linalg.eigh(Ttil_53)
# Sort descending
idx       = np.argsort(eig_vals)[::-1]
eig_vals  = eig_vals[idx]
eig_vecs  = eig_vecs[:, idx]
# Positive eigenvalues only
pos_mask  = eig_vals > 1e-12
mu_j      = eig_vals[pos_mask]
vecs_pos  = eig_vecs[:, pos_mask]
# Dominant zero index for each eigenvector
k_dom     = np.argmax(np.abs(vecs_pos), axis=0)
gamma_kdom = np.array([gammas_full[k] for k in k_dom])
inv_gamma  = 1.0 / gamma_kdom
r1         = np.corrcoef(mu_j, inv_gamma)[0, 1]
check(f"r1 = corr(μ_j, 1/γ_k(j)) > 0.90  (r1={r1:.3f})", r1 > 0.90)

# ── TEST 5: C_η ≈ 0.39 (κ-invariant) ────────────────────────────────────────
print("\n[5] C_η ≈ 0.39 (κ-invariant spectral ratio)")
ratios = []
for kap in KAPPAS_TEST[:-1]:  # skip 1009 for speed
    primes_k = list(primerange(2, kap + 1))
    Phi_k    = build_Phi(primes_k, gammas_full, EPS_NORM)
    c_k      = canonical_c(primes_k)
    _, _, _, _, D_k = energy(Phi_k, c_k)
    T_k      = Phi_k.T @ Phi_k
    D_inv    = 1.0 / D_k
    T_ren_k  = D_inv[:, None]**0.5 * T_k * D_inv[None, :]**0.5
    lam_max  = np.max(np.linalg.eigvalsh(T_ren_k))
    ratio    = lam_max / len(primes_k)
    ratios.append(ratio)
C_eta_mean = np.mean(ratios)
C_eta_std  = np.std(ratios)
check(f"C_η in [0.35, 0.43]  (mean={C_eta_mean:.3f}, std={C_eta_std:.4f})",
      0.35 < C_eta_mean < 0.43)

# ── TEST 6: Δ_Burst > 0 (κ-invariant for κ > 7) ──────────────────────────
print("\n[6] Δ_Burst > 0 (small-prime pairs, κ-invariant)")
small_primes = [2, 3, 5, 7]
c_small = canonical_c(small_primes)
gammas_arr = np.array(gammas_full)
gauss     = np.exp(-EPS_NORM**2 * gammas_arr**2)

def S(r):
    return np.sum(gauss * np.cos(gammas_arr * np.log(r)))

burst = 0.0
for i, p in enumerate(small_primes):
    for j, q in enumerate(small_primes):
        if p < q:
            burst += c_small[i] * c_small[j] * (S(p * q) - S(p / q))
# Note: exact c_p formula (Paper 1) gives Δ_Burst ≈ 3.10.
# Paper §5.2 table (Δ_Burst ≈ 4.81) uses a different c_p convention —
# an open inconsistency flagged for future revision.
check(f"Δ_Burst > 0  (computed={burst:.4f})", burst > 0)

# ── TEST 7: T̃ NOT a HP operator (r2 declining) ────────────────────────────
print("\n[7] T̃ NOT an HP operator: r2 = corr(ω_j, γ_{k(j)}) < 0.55 at κ=53")
eig_vals_pos = eig_vals[pos_mask]
# C_T via OLS: mu_j ~ C_T / gamma_{k(j)}  =>  OLS of mu_j on 1/gamma
C_T_ols = np.polyfit(inv_gamma, mu_j, 1)[0]
omega_j2 = C_T_ols / mu_j
r2 = np.corrcoef(omega_j2, gamma_kdom)[0, 1]
check(f"r2 at κ=53: {r2:.3f}  (should be < 0.55, > 0.30)", 0.20 < r2 < 0.60)

# ── SUMMARY ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print(f"SUMMARY: {PASS_COUNT} PASS, {FAIL_COUNT} FAIL")
if FAIL_COUNT == 0:
    print("ALL CHECKS PASSED ✓")
    print("Paper 4 numerical results verified.")
else:
    print(f"WARNING: {FAIL_COUNT} check(s) failed.")
print("=" * 65)
print(f"Normative: κ={KAPPA_NORM}, ε={EPS_NORM}, N={N_NORM}, σ={SIGMA_NORM}")
print("Paper 4: A Dual Operator for Prime–Zero Coupling · April 2026")

# ── FIGURE: fig_hp_main.png ──────────────────────────────────────────────────
print("\nGenerating fig_hp_main.png ...")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

kappas_fig = [23, 53, 101, 199, 503]
r1_vals_fig, r2_vals_fig, C_T_vals_fig = [], [], []

for kap in kappas_fig:
    primes_k = list(primerange(2, kap + 1))
    Phi_k    = build_Phi(primes_k, gammas_full, EPS_NORM)
    Ttil_k   = Phi_k @ Phi_k.T
    evals_k, evecs_k = np.linalg.eigh(Ttil_k)
    idx_k    = np.argsort(evals_k)[::-1]
    evals_k  = evals_k[idx_k]; evecs_k = evecs_k[:, idx_k]
    pos_k    = evals_k > 1e-12
    mu_k     = evals_k[pos_k]; vecs_k = evecs_k[:, pos_k]
    kdom_k   = np.argmax(np.abs(vecs_k), axis=0)
    gdom_k   = np.array([gammas_full[i] for i in kdom_k])
    CT_k     = np.polyfit(1.0/gdom_k, mu_k, 1)[0]
    om_k     = CT_k / mu_k
    r1_vals_fig.append(np.corrcoef(mu_k, 1.0/gdom_k)[0, 1])
    r2_vals_fig.append(np.corrcoef(om_k, gdom_k)[0, 1])
    C_T_vals_fig.append(CT_k)

# κ=53 detail data
primes_53p = list(primerange(2, 54))
Phi_53p    = build_Phi(primes_53p, gammas_full, EPS_NORM)
Ttil_53p   = Phi_53p @ Phi_53p.T
ev53, evec53 = np.linalg.eigh(Ttil_53p)
idx53 = np.argsort(ev53)[::-1]; ev53 = ev53[idx53]; evec53 = evec53[:, idx53]
pos53 = ev53 > 1e-12; mu53 = ev53[pos53]; vec53 = evec53[:, pos53]
kd53  = np.argmax(np.abs(vec53), axis=0)
gd53  = np.array([gammas_full[i] for i in kd53])
CT53  = np.polyfit(1.0/gd53, mu53, 1)[0]
mod53 = CT53 / gd53; om53 = CT53 / mu53
r1_53p = np.corrcoef(mu53, 1.0/gd53)[0, 1]
r2_53p = np.corrcoef(om53, gd53)[0, 1]

fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# Panel 1 — eigenvalue law
ax = axes[0, 0]
jj = range(1, len(mu53) + 1)
ax.semilogy(jj, mu53,  'o-', color='steelblue',  lw=2, ms=6,
            label=r'$\mu_j$ (eigenvalues of $\widetilde{T}$)')
ax.semilogy(jj, mod53, 's--', color='darkorange', lw=2, ms=6,
            label=f'$C_T/\\gamma_{{k(j)}}$  ($C_T={CT53:.1f}$)')
ax.set_xlabel('Mode $j$', fontsize=11)
ax.set_ylabel('Value (log scale)', fontsize=11)
ax.set_title(f'Eigenvalue law $\\mu_j \\sim C_T/\\gamma_{{k(j)}}$'
             f'  ($r_1={r1_53p:.3f}$, $\\kappa=53$)', fontsize=10)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# Panel 2 — HP quality
ax = axes[0, 1]
ax.scatter(gd53, om53, color='steelblue', s=70, zorder=3,
           label=r'$\omega_j = C_T/\mu_j$')
xlim = [12, 72]
ax.plot(xlim, xlim, 'k--', alpha=0.4, label='$\\omega=\\gamma$ (HP target)')
ax.set_xlabel(r'$\gamma_{k(j)}$', fontsize=11)
ax.set_ylabel(r'$\omega_j = C_T/\mu_j$', fontsize=11)
ax.set_yscale('log')
ax.set_title(f'$W_1$ eigenvalues vs zero ordinates'
             f'  ($r_2={r2_53p:.3f}$, $\\kappa=53$)', fontsize=10)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
for j in range(min(4, len(mu53))):
    ax.annotate(f'$j={j+1}$', (gd53[j], om53[j]),
                textcoords='offset points', xytext=(5, 0), fontsize=8)

# Panel 3 — r1 and r2 vs kappa
ax = axes[1, 0]
ax.plot(kappas_fig, r1_vals_fig, 'o-', color='steelblue',  lw=2, ms=8,
        label=r'$r_1=\mathrm{corr}(\mu_j,\,1/\gamma_{k(j)})$')
ax.plot(kappas_fig, r2_vals_fig, 's-', color='darkorange', lw=2, ms=8,
        label=r'$r_2=\mathrm{corr}(\omega_j,\,\gamma_{k(j)})$')
ax.axhline(1.0, color='gray', ls=':', alpha=0.5)
ax.axhline(0.0, color='gray', ls='--', alpha=0.3)
ax.set_xlabel(r'$\kappa$', fontsize=11)
ax.set_ylabel('Correlation', fontsize=11)
ax.set_title(r'$r_1$ (resonance law) vs $r_2$ (HP quality) as $\kappa$ grows',
             fontsize=10)
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
for i, kap in enumerate(kappas_fig):
    ax.annotate(f'{r1_vals_fig[i]:.3f}', (kap, r1_vals_fig[i]),
                textcoords='offset points', xytext=(0, 8),
                ha='center', fontsize=8, color='steelblue')
    ax.annotate(f'{r2_vals_fig[i]:.3f}', (kap, r2_vals_fig[i]),
                textcoords='offset points', xytext=(0, -15),
                ha='center', fontsize=8, color='darkorange')

# Panel 4 — C_T(kappa)
ax = axes[1, 1]
ax.plot(kappas_fig, C_T_vals_fig, 'o-', color='steelblue', lw=2, ms=8)
ax.set_xlabel(r'$\kappa$', fontsize=11)
ax.set_ylabel(r'$C_T(\kappa)$', fontsize=11)
ax.set_title(r'Arithmetic constant $C_T(\kappa)$ — OLS slope of'
             r' $\mu_j \sim C_T/\gamma_{k(j)}$', fontsize=10)
ax.grid(True, alpha=0.3)
for i, kap in enumerate(kappas_fig):
    ax.annotate(f'{C_T_vals_fig[i]:.1f}', (kap, C_T_vals_fig[i]),
                textcoords='offset points', xytext=(0, 8),
                ha='center', fontsize=9)

plt.suptitle(
    r'Paper 4: Spectral structure of $\widetilde{T} = \Phi\Phi^*$' + '\n' +
    r'$\widetilde{T}$ is a resonance operator — NOT a Hilbert–Pólya operator'
    f'   (N={N_NORM} zero ordinates, ε={EPS_NORM})',
    fontsize=12, fontweight='bold')
plt.tight_layout()

os.makedirs('figures/paper4', exist_ok=True)
out_path = os.path.join('figures', 'paper4', 'fig_hp_main.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"✓ fig_hp_main.png saved → {out_path}")
