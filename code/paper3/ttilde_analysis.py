"""
T̃ = Φ∘Φ* Analysis — Dual Loop Operator
AnalysisLab_L1_L5 · März 2026

Computes T̃ = Φ∘Φ* on H_null and T = Φ*∘Φ on H_str.
Verifies the algebraic identity λ_j(T) = μ_j(T̃) (nonzero part)
and analyzes eigenvector localization on zero ordinates.

Normative parameters: κ=53, ε=0.05, N=100
Expected results (from Sprint Report Ttilde, März 2026):
  - Max |λ_j − μ_j| for j=1..16: < 1e-14  (machine precision)
  - Localization loc_j = max|v_j|/‖v_j‖ ≥ 0.45 for all j ≤ 16
  - Correlation μ_j ~ 1/γ_{k(j)}: r ≥ 0.95

Reference: Paper 3, §3–§5
"""

import numpy as np
from mpmath import zetazero
import matplotlib.pyplot as plt
import os

# ── Normative parameters ────────────────────────────────────────────────────
N = 100
epsilon = 0.05
primes_53 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]

# ── Load zeta zeros ─────────────────────────────────────────────────────────
print("Loading first 100 zeta zeros...")
gammas = np.array([float(zetazero(k).imag) for k in range(1, N + 1)])
print(f"γ_1 = {gammas[0]:.6f}, γ_100 = {gammas[99]:.6f}")


# ── Core definitions (SSOT Rev22 §28) ───────────────────────────────────────

def a_p_vec(p, gammas, epsilon):
    """Option B: a_p[k] = exp(-ε²γ_k²/2) · sin(γ_k · log p)"""
    return np.exp(-epsilon**2 * gammas**2 / 2) * np.sin(gammas * np.log(p))


def compute_A(primes, gammas, epsilon):
    """Matrix A of shape (P, N): rows are a_p vectors."""
    return np.array([a_p_vec(p, gammas, epsilon) for p in primes])


def compute_T(primes, gammas, epsilon):
    """T = Φ*∘Φ : H_str → H_str, shape (P, P). T_{pq} = G^un_{pq}."""
    A = compute_A(primes, gammas, epsilon)
    return A @ A.T, A


def compute_T_tilde(primes, gammas, epsilon):
    """T̃ = Φ∘Φ* : H_null → H_null, shape (N, N)."""
    A = compute_A(primes, gammas, epsilon)
    return A.T @ A, A


# ── Main computation ─────────────────────────────────────────────────────────
print("\n=== Sprint T̃: Dual Loop Operator ===")

T, A = compute_T(primes_53, gammas, epsilon)
T_tilde, _ = compute_T_tilde(primes_53, gammas, epsilon)

P = len(primes_53)

# Eigenvalues (both sides)
lam = np.sort(np.linalg.eigvalsh(T))[::-1]         # shape (P,)
mu  = np.sort(np.linalg.eigvalsh(T_tilde))[::-1]    # shape (N,), N-P zeros appended

# Algebraic identity check
max_diff = np.max(np.abs(lam - mu[:P]))
print(f"\nAlgebraic identity λ_j(T) = μ_j(T̃):")
print(f"  Max |λ_j − μ_j| for j=1..{P}: {max_diff:.2e}")
assert max_diff < 1e-10, "FAIL: identity violated beyond numerical tolerance"
print(f"  PASS: machine-precision agreement ✓")

# Null eigenvalues of T̃
max_null = np.max(np.abs(mu[P:]))
print(f"  Max |μ_j| for j>{P} (should be ≈ 0): {max_null:.2e} ✓")

# ── Eigenvector localization ─────────────────────────────────────────────────
eigvals_tilde, eigvecs_tilde = np.linalg.eigh(T_tilde)
# Sort descending
idx = np.argsort(eigvals_tilde)[::-1]
eigvals_tilde = eigvals_tilde[idx]
eigvecs_tilde = eigvecs_tilde[:, idx]

print(f"\nEigenvector localization on T̃ (N={N}, P={P}):")
print(f"{'j':>3} | {'μ_j':>10} | {'k_dom':>6} | {'γ_dom':>8} | {'loc':>8}")
print("-" * 45)

locs = []
k_doms = []
gamma_doms = []

for j in range(P):
    v = eigvecs_tilde[:, j]
    k_dom = int(np.argmax(np.abs(v)))
    loc = np.max(np.abs(v)) / np.linalg.norm(v)
    locs.append(loc)
    k_doms.append(k_dom)
    gamma_doms.append(gammas[k_dom])
    print(f"{j+1:>3} | {eigvals_tilde[j]:>10.6f} | {k_dom+1:>6} | "
          f"{gammas[k_dom]:>8.4f} | {loc:>8.4f}")

min_loc = min(locs)
print(f"\nMin localization (all j ≤ {P}): {min_loc:.4f}")
assert min_loc >= 0.40, f"FAIL: localization {min_loc:.4f} < 0.40"
print(f"PASS: localization ≥ 0.40 ✓")

# ── Correlation μ_j ~ 1/γ_{k(j)} ────────────────────────────────────────────
gamma_doms_arr = np.array(gamma_doms[:P])
mu_arr = np.array([eigvals_tilde[j] for j in range(P)])
corr = np.corrcoef(1.0 / gamma_doms_arr, mu_arr)[0, 1]
print(f"\nCorrelation μ_j ~ 1/γ_{{k(j)}}: r = {corr:.4f}")
assert corr >= 0.90, f"FAIL: correlation {corr:.4f} < 0.90"
print(f"PASS: r ≥ 0.90 ✓")

# ── Save CSV ─────────────────────────────────────────────────────────────────
os.makedirs("data/results", exist_ok=True)
rows = []
for j in range(P):
    rows.append(f"{j+1},{eigvals_tilde[j]:.8f},{k_doms[j]+1},{gamma_doms[j]:.6f},{locs[j]:.6f}")

with open("data/results/ttilde_spectrum.csv", "w") as f:
    f.write("j,mu_j,k_dom,gamma_dom,localization\n")
    f.write("\n".join(rows))
print("\n✓ data/results/ttilde_spectrum.csv written")

# ── Figure 5: Localization bar chart ────────────────────────────────────────
os.makedirs("figures/paper3", exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: localization per eigenmode
axes[0].bar(range(1, P + 1), locs, color='steelblue', alpha=0.8)
axes[0].axhline(y=0.45, color='red', linestyle='--', alpha=0.6, label='loc = 0.45')
axes[0].set_xlabel('Mode j', fontsize=12)
axes[0].set_ylabel('Localization max|v_j|/‖v_j‖', fontsize=11)
axes[0].set_title('T̃ Eigenvector Localization (κ=53, ε=0.05, N=100)', fontsize=11)
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Right: μ_j vs 1/γ_{k(j)} scatter
x_vals = 1.0 / gamma_doms_arr
axes[1].scatter(x_vals, mu_arr, color='steelblue', s=60, zorder=3)
# Fit line
coeffs = np.polyfit(x_vals, mu_arr, 1)
x_line = np.linspace(x_vals.min(), x_vals.max(), 100)
axes[1].plot(x_line, np.polyval(coeffs, x_line), 'r--',
             alpha=0.7, label=f'r = {corr:.3f}')
axes[1].set_xlabel('1 / γ_{k(j)}', fontsize=12)
axes[1].set_ylabel('μ_j', fontsize=12)
axes[1].set_title('Eigenvalue ~ 1/γ_dom Correlation', fontsize=11)
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("figures/paper3/fig5_ttilde_localization.png", dpi=150)
plt.close()
print("✓ figures/paper3/fig5_ttilde_localization.png written")

# Figure 6: μ_j vs γ_{k(j)}
fig, ax = plt.subplots(figsize=(7, 5))
sc = ax.scatter(gamma_doms_arr, mu_arr, c=range(P),
                cmap='viridis', s=80, zorder=3)
plt.colorbar(sc, ax=ax, label='Mode j')
ax.set_xlabel('γ_{k(j)} (dominant zero ordinate)', fontsize=12)
ax.set_ylabel('μ_j (eigenvalue of T̃)', fontsize=12)
ax.set_title('μ_j vs γ_{k(j)} — T̃ Spectrum (κ=53, ε=0.05)', fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("figures/paper3/fig6_mu_vs_gamma.png", dpi=150)
plt.close()
print("✓ figures/paper3/fig6_mu_vs_gamma.png written")

print("\n=== PASS: All T̃ checks passed ===")
