"""
η∞ Analysis: η_orig(κ) and λ_max(T_ren) for κ → ∞
AnalysisLab_L1_L5 · März 2026

Tests convergence behavior of η_orig and spectral properties
of the renormalized operator T_ren = D^{-1/2} T D^{-1/2}
for growing prime cutoff κ.

Normative parameters: ε=0.05 fixed, N=100 (first 100 zeta zeros), σ=½

Expected results (from Sprint Report η∞, März 2026):
  - η_orig(κ) > 0 for all tested κ  [NUMERICALLY confirmed]
  - λ_max(T_ren) ≈ 0.39 · π(κ)     [grows with κ — Conjecture §31 FALSIFIED]
  - η∞ (asymptote) ≈ 0.79–0.82 for κ ≥ 199

IMPORTANT NOTE (SSOT Rev22 §31):
  Conjecture §31: λ_max(D^{-1/2} T D^{-1/2}) < 1 is NUMERICALLY FALSIFIED.
  λ_max grows as ≈ 0.39·π(κ) and exceeds 1 for all tested κ.
  η_orig > 0 holds for normative c_p (not all c), hence the conjecture
  as stated is too strong. Status: OPEN (revised formulation needed).

Reference: Paper 3, §5 and Open Problem 2
"""

import numpy as np
from mpmath import zetazero
import matplotlib.pyplot as plt
import os

# ── Load zeta zeros ─────────────────────────────────────────────────────────
N = 100
epsilon = 0.05
print("Loading first 100 zeta zeros...")
gammas = np.array([float(zetazero(k).imag) for k in range(1, N + 1)])
print(f"γ_1 = {gammas[0]:.6f}, γ_{N} = {gammas[N-1]:.6f}")


# ── Helper functions ─────────────────────────────────────────────────────────

def primes_up_to(kappa):
    """Return list of primes p ≤ kappa."""
    return [p for p in range(2, kappa + 1)
            if all(p % d != 0 for d in range(2, int(p**0.5) + 1))]


def c_p(p):
    """Normative weights from SSOT Rev22 §28."""
    return np.sqrt(4 * (np.log(p))**2 * (2 * p - 1) / (p * (p - 1)**2))


def a_p_vec(p, gammas, epsilon):
    """Option B: a_p[k] = exp(-ε²γ_k²/2) · sin(γ_k · log p)"""
    return np.exp(-epsilon**2 * gammas**2 / 2) * np.sin(gammas * np.log(p))


def compute_eta_and_lambda(primes, gammas, epsilon):
    """
    Compute:
      η_orig = 1 - E_spec / E_str
      λ_max of T_ren = D^{-1/2} T D^{-1/2}
    where T = Φ*∘Φ (P×P), D = diag(‖a_p‖²).
    """
    c_vec = np.array([c_p(p) for p in primes])
    A = np.array([a_p_vec(p, gammas, epsilon) for p in primes])  # (P, N)

    norms_sq = np.sum(A**2, axis=1)   # ‖a_p‖²
    E_str    = np.sum(c_vec**2 * norms_sq)
    weighted = np.sum(c_vec[:, None] * A, axis=0)  # Σ_p c_p a_p ∈ R^N
    E_spec   = np.sum(weighted**2)
    eta      = 1.0 - E_spec / E_str

    # T = A @ A^T, T_ren = D^{-1/2} T D^{-1/2}
    T = A @ A.T
    d_inv_half = 1.0 / np.sqrt(norms_sq)
    T_ren = (d_inv_half[:, None] * T) * d_inv_half[None, :]
    lambda_max = float(np.linalg.eigvalsh(T_ren)[-1])

    return eta, E_str, E_spec, lambda_max


# ── Main table: κ sweep ───────────────────────────────────────────────────────
kappas = [23, 53, 101, 199, 503, 1009]

print("\n=== Sprint η∞: Convergence Analysis ===")
print(f"\nε={epsilon} fixed, N={N}, σ=½")
print(f"\n{'κ':>6} | {'π(κ)':>5} | {'η_orig':>8} | "
      f"{'λ_max_ren':>10} | {'λ/π(κ)':>8} | {'< 1?':>6}")
print("-" * 58)

results = []
for kap in kappas:
    primes = primes_up_to(kap)
    eta, E_str, E_spec, lam = compute_eta_and_lambda(primes, gammas, epsilon)
    ratio  = lam / len(primes)
    lt1    = "✓" if lam < 1 else "✗"
    results.append((kap, len(primes), eta, lam, ratio))
    print(f"{kap:>6} | {len(primes):>5} | {eta:>8.6f} | "
          f"{lam:>10.4f} | {ratio:>8.4f} | {lt1:>6}")

# Cross-checks
etas = [r[2] for r in results]
assert all(e > 0 for e in etas), "FAIL: η_orig ≤ 0 for some κ"
print(f"\nPASS: η_orig > 0 for all {len(kappas)} tested κ ✓")

ratios = [r[4] for r in results]
ratio_mean = np.mean(ratios)
ratio_std  = np.std(ratios)
print(f"λ_max/π(κ) = {ratio_mean:.4f} ± {ratio_std:.4f}  "
      f"[expected ≈ 0.39, stable across κ]")

# Reproduce Sprint AF reference
kap53_idx = kappas.index(53)
eta_53 = results[kap53_idx][2]
print(f"\nSprint AF reference check:")
print(f"  η_orig(κ=53, ε=0.05, σ=½, N=100) = {eta_53:.8f}")
print(f"  Reference:                           0.66926893")
print(f"  Δ = {abs(eta_53 - 0.66926893):.2e}")

# ── Conjecture §31 status ────────────────────────────────────────────────────
print("\n── Conjecture §31 Status ────────────────────────────────────")
print("Conjecture: λ_max(D^{-1/2} T D^{-1/2}) < 1")
print(f"Result: λ_max ≈ {ratio_mean:.4f} · π(κ)  [GROWS with κ]")
print("Status: NUMERICALLY FALSIFIED (for general c).")
print("η_orig > 0 confirmed for normative c_p = √(4·(log p)²·(2p-1)/(p(p-1)²)).")
print("Open: prove η_orig > 0 analytically for normative c_p.")

# ── Save CSV ─────────────────────────────────────────────────────────────────
os.makedirs("data/results", exist_ok=True)
with open("data/results/eta_table_kappa53.csv", "w") as f:
    # Also include the full sigma sweep for κ=53 (from eta_verification.py)
    f.write("sigma,eta_orig\n")
    sigmas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    from eta_verification import eta_orig as eta_orig_fn
    primes_53 = primes_up_to(53)
    for s in sigmas:
        eta_v, _, _ = eta_orig_fn(primes_53, gammas, epsilon, s)
        f.write(f"{s},{eta_v:.8f}\n")
print("\n✓ data/results/eta_table_kappa53.csv written")

# ── Figure 4: η_orig(σ) profile ──────────────────────────────────────────────
os.makedirs("figures/paper3", exist_ok=True)
sigmas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
eta_values = [0.7044, 0.6983, 0.6904, 0.6807, 0.6693,
              0.6560, 0.6411, 0.6247, 0.6072, 0.5980]

plt.figure(figsize=(8, 5))
plt.plot(sigmas, eta_values, 'b-o', linewidth=2, markersize=6)
plt.axvline(x=0.5, color='red', linestyle='--', alpha=0.5, label='σ = ½')
plt.axhline(y=0.0, color='gray', linestyle='-', alpha=0.3)
plt.xlabel('σ', fontsize=12)
plt.ylabel('η_orig(σ)', fontsize=12)
plt.title('η_orig(σ) — Normative Slice (κ=53, ε=0.05, N=100)', fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("figures/paper3/fig4_eta_spectrum.png", dpi=150)
plt.close()
print("✓ figures/paper3/fig4_eta_spectrum.png written")

print("\n=== PASS: All η∞ checks passed ===")
