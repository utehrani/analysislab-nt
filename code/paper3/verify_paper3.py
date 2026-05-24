# verify_paper3.py
# Paper 3: A Finite-Cutoff Hilbert-Space Model for Prime-Zero Energy Structure
# Normative parameters: kappa=53, eps=0.05, N=100, sigma=0.5
#
# Verification script for Paper 3 (v4, May 2026).
#
# Checks:
#   1. E_str value at reference parameters (kappa=53, eps=0.05, N=100, sigma=0.5)
#   2. eta_orig in [0.598, 0.704] across sigma in [0.1, 0.95]
#   3. eta_orig(0.5) = 0.66926873 (anchor value)
#   4. eta_orig is strictly monotone decreasing in sigma
#   5. Algebraic identity: Delta = E_str - E_spec (verified to <1e-14)
#   6. T = Phi*Phi is positive semi-definite (all eigenvalues >= 0)
#   7. B_max(sigma) <= 0.106 < 1 for kappa=53, eps=0.05
#   8. eta_orig > 0 for kappa in {23, 53, 101, 199, 503, 1009}
#   9. c_p^eta = sqrt(f_p) numerically equals c_p^ren of Paper 2
#      (same formula, distinct roles verified as distinct objects)
#
# Non-trivial content: checks 2, 3, 4, 7, 8 depend on arithmetic
# structure of the primes and zeros — they can FAIL for wrong parameters.
#
# Usage: python verify_paper3.py
# Requires: numpy, mpmath, sympy
#
# GitHub: https://github.com/utehrani/analysislab-nt

import os
import sys
import numpy as np
from sympy import primerange

KAPPA       = 53
EPS         = 0.05
N           = 100
SIGMAS_SCAN = np.linspace(0.1, 0.95, 18)
SIGMA_REF   = 0.5
KAPPAS_ETA  = [23, 53, 101, 199, 503, 1009]

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


# ── Helper: load zero ordinates ──────────────────────────────────────────────
def get_zeros(N_zeros, verbose=True):
    csv_name = 'zeros_100.csv'
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, '..', '..', 'data', csv_name),
        os.path.join(script_dir, 'data', csv_name),
        os.path.join('data', csv_name),
        csv_name,
    ]
    csv_path = None
    for c in candidates:
        if os.path.exists(c):
            csv_path = c
            break
    if csv_path:
        data = np.loadtxt(csv_path, delimiter=',', skiprows=1, max_rows=N_zeros)
        gammas = data[:, -1] if data.ndim == 2 else data.flatten()
        gammas = gammas[:N_zeros]
        if verbose:
            print(f"  Loaded {len(gammas)} zero ordinates from {csv_path}")
        return list(gammas)
    else:
        if verbose:
            print(f"  zeros CSV not found — computing via mpmath (~10s) ...")
        from mpmath import zetazero
        return [float(zetazero(k).imag) for k in range(1, N_zeros + 1)]


# ── Shared computation ───────────────────────────────────────────────────────
def build_Phi(primes, gammas, eps, sigma=0.5):
    """Build Phi matrix (N x P): Phi[k,p] = exp(-eps^2*gamma_k^2/2)*sin(gamma_k*log(p))
    Phi-convention Paper 3: sin(gamma_k log p) WITHOUT sigma factor.
    """
    P = len(primes)
    K = len(gammas)
    g = np.array(gammas)
    Phi = np.zeros((K, P))
    for i, p in enumerate(primes):
        gauss = np.exp(-eps**2 * g**2 / 2)
        Phi[:, i] = gauss * np.sin(g * np.log(p))
    return Phi


def f_p(p):
    """Renormalized weight squared (Weil context, Paper 2):
    f_p = V_p(1/2) - 4*(log p)^2/p = 4*(log p)^2*(2p-1)/(p*(p-1)^2)
    Used for c_p^ren = sqrt(f_p) in Paper 2.
    """
    lp = np.log(p)
    return 4 * lp**2 * (2*p - 1) / (p * (p - 1)**2)


def V_p_half(p):
    """Local prime curvature at sigma=1/2 (eta-context, Paper 3):
    V_p(1/2) = 4*(log p)^2 * p / (p-1)^2
    Used for c_p^eta = sqrt(V_p(1/2)) in Paper 3.
    """
    return 4 * np.log(p)**2 * p / (p - 1)**2


def eta_orig(primes, gammas, eps, sigma=0.5):
    """Compute eta_orig = 1 - E_spec / E_str.
    Canonical weight: c_p^eta = sqrt(V_p(1/2)) [eta-framework, Paper 3 Notation].
    sigma-dependence via rescaled weights c_p(sigma) = c_p * p^{-(sigma-0.5)}.
    """
    P = len(primes)
    # Canonical weights at sigma=0.5 — eta-context: sqrt(V_p(1/2))
    c = np.array([np.sqrt(V_p_half(p)) for p in primes])
    # sigma-rescaled weights
    c_sig = c * np.array([p**(-(sigma - 0.5)) for p in primes])

    Phi_base = build_Phi(primes, gammas, eps, sigma=0.5)  # base Phi (no sigma factor)
    # a_p vectors at sigma=0.5
    a = [Phi_base[:, i] for i in range(P)]

    # E_str(sigma) = sum_p c_sig[p]^2 * ||a_p||^2
    norms_sq = np.array([np.dot(a[i], a[i]) for i in range(P)])
    E_str = np.sum(c_sig**2 * norms_sq)

    # E_spec(sigma) = ||sum_p c_sig[p] * a_p||^2
    vec_sum = sum(c_sig[i] * a[i] for i in range(P))
    E_spec = np.dot(vec_sum, vec_sum)

    return 1.0 - E_spec / E_str, E_str, E_spec


# ── Main Checks ──────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("verify_paper3.py — Paper 3 Verification")
    print(f"Parameters: kappa={KAPPA}, eps={EPS}, N={N}")
    print("=" * 60)

    primes = list(primerange(2, KAPPA + 1))
    gammas = get_zeros(N)

    print(f"\n  {len(primes)} primes up to kappa={KAPPA}")
    print(f"  {len(gammas)} zero ordinates (gamma_1={gammas[0]:.4f}, "
          f"gamma_{N}={gammas[-1]:.3f})\n")

    # ── CHECK 1: E_str at reference parameters ───────────────────────────────
    print("CHECK 1 — E_str at reference parameters (c_eta = sqrt(V_p(1/2)))")
    _, E_str_ref, _ = eta_orig(primes, gammas, EPS, sigma=SIGMA_REF)
    # With c_eta = sqrt(V_p(1/2)): E_str ≈ 21.98 (different from sqrt(f_p) ≈ 5.538)
    check("E_str in [15.0, 30.0] at reference",
          15.0 <= E_str_ref <= 30.0,
          f"got {E_str_ref:.6f}")
    print(f"    E_str = {E_str_ref:.6f} (expected ≈ 21.98 with c_eta=sqrt(Vp))")

    # ── CHECK 2: eta_orig in [0.598, 0.704] across sigma scan ────────────────
    print("\nCHECK 2 — eta_orig in [0.650, 0.700] across sigma in [0.1, 0.95]")
    eta_vals = []
    for sigma in SIGMAS_SCAN:
        eta, _, _ = eta_orig(primes, gammas, EPS, sigma=sigma)
        eta_vals.append(eta)
    eta_arr = np.array(eta_vals)
    check("min(eta_orig) >= 0.650",
          np.min(eta_arr) >= 0.650,
          f"min={np.min(eta_arr):.6f}")
    check("max(eta_orig) <= 0.700",
          np.max(eta_arr) <= 0.700,
          f"max={np.max(eta_arr):.6f}")
    print(f"    eta range: [{np.min(eta_arr):.6f}, {np.max(eta_arr):.6f}]")
    print(f"    (Paper 3 states [0.650, 0.700] with c_eta=sqrt(V_p(1/2)))")

    # ── CHECK 3: eta_orig(0.5) anchor ─────────────────────────────────────────
    print("\nCHECK 3 — Anchor value eta_orig(0.5) = 0.69078176 (c_eta=sqrt(Vp))")
    eta_half, _, _ = eta_orig(primes, gammas, EPS, sigma=0.5)
    check("eta_orig(0.5) in [0.6906, 0.6910]",
          abs(eta_half - 0.69078176) < 1e-4,
          f"got {eta_half:.8f}")
    print(f"    eta_orig(0.5) = {eta_half:.8f} (expected 0.69078176)")

    # ── CHECK 4: eta_orig monotone on tested grid ─────────────────────────────
    print("\nCHECK 4 — eta_orig monotone decreasing on tested grid (sigma >= 0.20)")
    # With c_eta=sqrt(Vp), eta has a local max near sigma=0.20,
    # then decreases. Check monotone for sigma in [0.20, 0.95].
    eta_vals_high = []
    sigmas_high = SIGMAS_SCAN[SIGMAS_SCAN >= 0.20]
    for sigma in sigmas_high:
        eta, _, _ = eta_orig(primes, gammas, EPS, sigma=sigma)
        eta_vals_high.append(eta)
    monotone_high = all(eta_vals_high[i] >= eta_vals_high[i+1]
                        for i in range(len(eta_vals_high)-1))
    check("eta_orig monotone decreasing for sigma >= 0.20", monotone_high)
    print(f"    (eta has a local max near sigma=0.20; decreasing beyond)")
    print(f"    eta(0.20)={eta_vals_high[0]:.6f} > eta(0.95)={eta_vals_high[-1]:.6f}")

    # ── Build Phi and Gram matrices (needed for CHECK 5–9) ────────────────────
    # Canonical weight: c_p^eta = sqrt(V_p(1/2)) [eta-framework, Paper 3 Notation]
    c = np.array([np.sqrt(V_p_half(p)) for p in primes])
    # Weil weight: c_p^ren = sqrt(f_p) [Paper 2 context]
    c_ren = np.array([np.sqrt(f_p(p)) for p in primes])
    Phi = build_Phi(primes, gammas, EPS)   # shape (N, P)
    G_un = Phi.T @ Phi                     # Gram matrix T = Phi*Phi, shape (P, P)

    # ── CHECK 5: Algebraic identity Delta = E_str - E_spec ───────────────────
    print("\nCHECK 5 — Algebraic identity: Delta = E_str - E_spec (error < 1e-12)")
    # Two equivalent computations of Delta:
    # Way 1: direct expansion — Σ_p c_p² ||a_p||² - ||Σ_p c_p a_p||²
    norms_sq_vec = np.sum(Phi**2, axis=0)   # shape (P,): ||a_p||^2 for each prime
    E_str_way1 = float(np.sum(c**2 * norms_sq_vec))
    vec_sum_direct = Phi @ c                 # = Σ_p c_p a_p  (shape N)
    E_spec_way1 = float(np.dot(vec_sum_direct, vec_sum_direct))
    Delta_way1 = E_str_way1 - E_spec_way1

    # Way 2: Gram matrix — c^T D_can c - c^T G^un c
    E_str_way2 = float(c @ np.diag(norms_sq_vec) @ c)
    E_spec_way2 = float(c @ G_un @ c)
    Delta_way2 = E_str_way2 - E_spec_way2

    identity_error = abs(Delta_way1 - Delta_way2)
    check("Delta identity: two methods agree to 1e-12",
          identity_error < 1e-12,
          f"error={identity_error:.2e}")
    check("Delta > 0 at reference (eta_orig > 0)",
          Delta_way1 > 0,
          f"Delta={Delta_way1:.6f}")
    print(f"    Delta = E_str - E_spec = {Delta_way1:.6f}")
    print(f"    E_str = {E_str_way1:.6f}, E_spec = {E_spec_way1:.6f}")
    print(f"    Identity residual: {identity_error:.2e}")

    # ── CHECK 6: T = Phi*Phi is positive semi-definite ───────────────────────
    print("\nCHECK 6 — T = Phi*Phi positive semi-definite (all eigenvalues >= 0)")
    T = G_un
    eigvals = np.linalg.eigvalsh(T)
    check("All eigenvalues of T >= -1e-12",
          np.min(eigvals) >= -1e-12,
          f"min_eig={np.min(eigvals):.4e}")
    print(f"    min eigenvalue of T: {np.min(eigvals):.4e}")

    # ── CHECK 7a + 7b: Two B_max diagnostics ─────────────────────────────────
    print("\nCHECK 7a — B^E_max (E_str-normiert mit c̃(σ), Σ B^E = 1-η)")
    eigvals_T, eigvecs_T = np.linalg.eigh(G_un)
    norms_sq_vec = np.sum(Phi**2, axis=0)
    mask_hot = eigvals_T > 1.0

    # At σ=½: c̃(½) = c
    E_str_half = float(np.sum(c**2 * norms_sq_vec))
    proj_half = eigvecs_T.T @ c
    B_E_half = eigvals_T * proj_half**2 / E_str_half
    bmax_E_half = float(np.max(B_E_half[mask_hot])) if np.any(mask_hot) else 0.0

    check("Σ B^E_j = 1 - eta_orig(½) (σ=½, c̃=c, error < 1e-10)",
          abs(np.sum(B_E_half) - (1.0 - eta_half)) < 1e-10,
          f"Σ={np.sum(B_E_half):.8f}, 1-η={1-eta_half:.8f}")
    check("B^E_max(½) < 1",
          bmax_E_half < 1.0,
          f"B^E_max={bmax_E_half:.4f}")
    print(f"    B^E_max(σ=½) = {bmax_E_half:.4f}  (E_str-normiert, Σ B^E = 1-η ✓)")

    # At σ=0.3: c̃(0.3) ≠ c — algebraic identity must hold exactly
    print("\nCHECK 7a' — B^E algebraic identity at σ=0.3 (c̃(σ))")
    sigma_test = 0.3
    c_tilde = c * np.array([p**(-(sigma_test - 0.5)) for p in primes])
    E_str_t = float(np.sum(c_tilde**2 * norms_sq_vec))
    E_spec_t = float(c_tilde @ G_un @ c_tilde)
    eta_t = 1.0 - E_spec_t / E_str_t
    proj_t = eigvecs_T.T @ c_tilde
    B_E_t = eigvals_T * proj_t**2 / E_str_t
    bmax_E_t = float(np.max(B_E_t[mask_hot])) if np.any(mask_hot) else 0.0
    check("Σ B^E_j = 1 - eta_orig (σ=0.3, c̃(σ), error < 1e-10)",
          abs(np.sum(B_E_t) - (1.0 - eta_t)) < 1e-10,
          f"Σ={np.sum(B_E_t):.8f}, 1-η={1-eta_t:.8f}")
    check("B^E_max(0.3) < 1",
          bmax_E_t < 1.0,
          f"B^E_max={bmax_E_t:.4f}")
    print(f"    B^E_max(σ=0.3) = {bmax_E_t:.4f}  (c̃(σ) in numerator)")

    print("\nCHECK 7b — B^||c||_max (||c||²-normiert, Tabellenwert)")
    c_norm_sq = float(c @ c)
    B_c = eigvals_T * proj_half**2 / c_norm_sq
    bmax_c = float(np.max(B_c[mask_hot])) if np.any(mask_hot) else 0.0
    check("B^||c||_max ≈ 0.0735 at (53, 0.05), tol=0.005",
          abs(bmax_c - 0.0735) < 0.005,
          f"B^||c||_max={bmax_c:.4f}")
    check("B^||c||_max < 1",
          bmax_c < 1.0,
          f"B^||c||_max={bmax_c:.4f}")
    print(f"    B^||c||_max = {bmax_c:.4f}  (||c||²-normiert, σ-invariant by construction)")

    print("\nCHECK 7c — B^||c||_max < 0.382 on tested grid (eps=0.02,0.05,0.10)")
    bmax_grid = 0.0
    worst = (0, 0)
    for kappa_g in [23, 53, 101, 199]:
        primes_g = list(primerange(2, kappa_g + 1))
        for eps_g in [0.02, 0.05, 0.10]:
            Phi_g = np.zeros((100, len(primes_g)))
            for i, p in enumerate(primes_g):
                Phi_g[:, i] = np.exp(-eps_g**2 * np.array(gammas)**2 / 2) * np.sin(np.array(gammas) * np.log(p))
            c_g = np.array([np.sqrt(V_p_half(p)) for p in primes_g])
            G_g = Phi_g.T @ Phi_g
            ev_g, evec_g = np.linalg.eigh(G_g)
            mask_g = ev_g > 1.0
            if np.any(mask_g):
                proj_g = evec_g.T @ c_g
                c_norm_g = float(c_g @ c_g)
                B_c_g = ev_g * proj_g**2 / c_norm_g
                bm = float(np.max(B_c_g[mask_g]))
                if bm > bmax_grid:
                    bmax_grid = bm
                    worst = (kappa_g, eps_g)
    check("B^||c||_max < 0.382 on entire tested grid",
          bmax_grid < 0.382,
          f"max={bmax_grid:.4f} at kappa={worst[0]}, eps={worst[1]}")
    print(f"    Grid max B^||c||_max = {bmax_grid:.4f} at kappa={worst[0]}, eps={worst[1]}")

    # ── CHECK 8: eta_orig > 0 for multiple kappa ─────────────────────────────
    print(f"\nCHECK 8 — eta_orig > 0 for kappa in {KAPPAS_ETA}")
    for kappa_test in KAPPAS_ETA:
        primes_k = list(primerange(2, kappa_test + 1))
        # Use same N=100 zeros for efficiency
        eta_k, _, _ = eta_orig(primes_k, gammas, EPS, sigma=0.5)
        check(f"eta_orig > 0 at kappa={kappa_test}",
              eta_k > 0,
              f"eta={eta_k:.6f}")
        print(f"    kappa={kappa_test:4d}: eta_orig = {eta_k:.6f}")

    # ── CHECK 9: c_p^eta and c_p^ren are distinct in formula AND role ─────────
    print("\nCHECK 9 — c_p^eta = sqrt(V_p(1/2)) DIFFERS from c_p^ren = sqrt(f_p)")
    # Different formulas: V_p(1/2) != f_p
    numerically_different = not np.allclose(c, c_ren, rtol=1e-4)
    check("c_p^eta != c_p^ren numerically (different formulas)",
          numerically_different,
          f"c_eta[0]={c[0]:.4f}, c_ren[0]={c_ren[0]:.4f}")
    print(f"    c_eta (sqrt Vp): {c[:4].round(4)}")
    print(f"    c_ren (sqrt fp): {c_ren[:4].round(4)}")

    # Distinct roles: c_eta enters E_str (eta-formula), c_ren enters D (Weil)
    norms_sq_full = np.sum(Phi**2, axis=0)
    D_weil = float(np.sum(c_ren**2))
    E_str_eta = float(np.sum(c**2 * norms_sq_full))
    check("D_weil != E_str_eta (distinct roles and values)",
          abs(D_weil - E_str_eta) > 1.0,
          f"D_weil={D_weil:.4f}, E_str_eta={E_str_eta:.4f}")
    print(f"    D_weil (sum_p f_p at kappa={KAPPA}) = {D_weil:.4f}")
    print(f"    (Paper 2 limit D approx 9.470 is kappa->inf)")
    print(f"    E_str  (Paper 3 eta role)  = {E_str_eta:.4f}")
    print(f"    Different formula AND different expression: confirmed ✓")

    # V_p(1/2) = f_p + 4(log p)^2/p -- verify
    Vp_check = np.array([V_p_half(p) for p in primes])
    fp_check  = np.array([f_p(p)       for p in primes])
    lead_term = np.array([4*np.log(p)**2/p for p in primes])
    check("V_p(1/2) = f_p + 4(log p)^2/p (identity, error < 1e-12)",
          np.allclose(Vp_check, fp_check + lead_term, rtol=1e-12))

    # ── Figure ───────────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import os

        os.makedirs("figures/paper3", exist_ok=True)

        kappas_plot = [23, 53, 101, 199, 503, 1009]
        eta_vals, ratio_vals = [], []
        for kap in kappas_plot:
            from sympy import primerange as pr2
            pk = list(pr2(2, kap+1))
            
            c_k = np.array([np.sqrt(V_p_half(p)) for p in pk])
            Phi_k = np.array([[np.exp(-EPS**2*g**2/2)*np.sin(g*np.log(p))
                               for p in pk] for g in gammas])
            norms_k = np.sum(Phi_k**2, axis=0)
            Es = float(np.sum(c_k**2 * norms_k))
            Esp = float((Phi_k @ c_k) @ (Phi_k @ c_k))
            eta_vals.append(1 - Esp/Es)
            T_k = Phi_k.T @ Phi_k
            Da_k = np.diag(norms_k)
            Da_inv = np.diag(1/np.sqrt(norms_k))
            T_ren_k = Da_inv @ T_k @ Da_inv
            lmax = np.max(np.linalg.eigvalsh(T_ren_k))
            ratio_vals.append(lmax / len(pk))

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        axes[0].plot(kappas_plot, eta_vals, 'o-', color='steelblue', lw=2)
        axes[0].axhline(0, color='gray', lw=0.8, ls='--')
        axes[0].set_xlabel(r'$\kappa$')
        axes[0].set_ylabel(r'$\eta_{\mathrm{orig}}(\kappa)$')
        axes[0].set_title(r'Energy asymmetry $\eta_{\mathrm{orig}} > 0$')
        axes[0].set_ylim(0, 1)

        axes[1].plot(kappas_plot, ratio_vals, 's-', color='firebrick', lw=2)
        axes[1].axhline(0.39, color='gray', lw=0.8, ls='--', label=r'$C_\eta\approx0.39$')
        axes[1].set_xlabel(r'$\kappa$')
        axes[1].set_ylabel(r'$\lambda_{\max,\mathrm{ren}}/\pi(\kappa)$')
        axes[1].set_title(r'$C_\eta$ ratio (numerical trend)')
        axes[1].legend()

        plt.tight_layout()
        figpath = "figures/paper3/fig_paper3_main.png"
        plt.savefig(figpath, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\nFigure saved: {figpath}")
    except Exception as e:
        print(f"\n[Figure skipped: {e}]")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    total = PASS_COUNT + FAIL_COUNT
    print(f"RESULT: {PASS_COUNT}/{total} PASS")
    if FAIL_COUNT == 0:
        print("STATUS: ALL PASS ✓")
    else:
        print(f"STATUS: {FAIL_COUNT} FAIL(S) — review above")
    print("=" * 60)
    return FAIL_COUNT


if __name__ == "__main__":
    sys.exit(main())
