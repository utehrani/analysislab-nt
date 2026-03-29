# AnalysisLab: Curvature, Energy and the Riemann Zeta Function

Three companion papers developing a geometric and operator-theoretic approach
to the Riemann zeta function, building a connected framework from local
curvature through the Weil functional to a finite-dimensional Hilbert-space
model of prime–zero energy structure.

**Author:** Ulrich Tehrani  
**License:** MIT  
**DOIs:**
 [Paper 1](https://doi.org/10.5281/zenodo.19025598) ·
[Paper 2](https://doi.org/10.5281/zenodo.19106992) ·
[Paper 3](https://doi.org/10.5281/zenodo.19307989)

---

## The Three Papers

### Paper 1 — A Curvature Decomposition of the Explicit Formula
**DOI:** [10.5281/zenodo.19025598](https://doi.org/10.5281/zenodo.19025598)  
**File:** [`papers/paper1/curvature_note_v7.pdf`](papers/paper1/)  
**Scripts:** [`code/paper1/`](code/paper1/)

**What it proves:**  
The second logarithmic derivative of the completed zeta function decomposes as

```
H_xi(sigma, t) = H_local(sigma, kappa) + H_dual(sigma, t, kappa)
```

where `H_local` contains all prime contributions and `H_dual` the zero
contributions. Three key results:

| Result | Status |
|--------|--------|
| `V_p(sigma) >= 0` for all primes p, sigma > 0 | **PROVED** (Lemma 1) |
| `H_local(1/2, kappa) ~ 2*(log kappa)^2 → ∞` | **PROVED** (Lemma 2) |
| `H_local(sigma, kappa) → C(sigma) < ∞` for sigma > 1/2 | **PROVED** (Lemma 3) |

The critical line `sigma = 1/2` is the unique phase boundary: divergence below,
convergence above. The connection to the Li coefficients (whose positivity is
equivalent to RH) is identified.

**Reproduce:**
```bash
python code/paper1/verify_v6.py
```

**Figures generated:** `fig1_H_local_divergence.png`, `fig2_sigma_profile.png`

---

### Paper 2 — From Local Curvature to the Weil Functional
**DOI:** [10.5281/zenodo.19106992](https://doi.org/10.5281/zenodo.19106992)  
**File:** [`papers/paper2/weil_v4.pdf`](papers/paper2/)  
**Scripts:** [`code/paper2/`](code/paper2/)

**Imports from Paper 1:** `H_local` divergence (R1)

**What it proves:**  
Explicit admissible test functions `g*_{sigma,eps}` for the Weil explicit
formula satisfy

```
W(g* * g̃*) = Z(g*) - H_local(sigma, kappa) + O(eps)
```

This isolates the remaining open problem as a single spectral inequality
`Z(g*_{1/2,eps}) >= H_local(1/2, kappa)`.

| Result | Status |
|--------|--------|
| `c_p^ren = f_p^{1/2} > 0` (renormalized weights) | **PROVED** |
| `D = sum_p (c_p^ren)^2 ≈ 9.471 < ∞` | **PROVED** (Convergence Lemma) |
| `W(g*,g*) = Z(g*) - H_local + O(eps)` | **PROVED** (Theorem 3.1) |
| `H_local^ren` has Mertens-scale sawtooth structure | **PROVED** |
| Bridge constant `D/(2*pi) ≈ 1.507` | **NUMERICAL** |

**Reproduce:**
```bash
python code/paper2/verify_paper2.py
```

**Figures generated:** `fig3_Weil_decomposition.png`, `fig4_sawtooth_Hren.png`

---

### Paper 3 — A Finite-Cutoff Hilbert-Space Model for Prime–Zero Energy Structure
**DOI:** [10.5281/zenodo.19307989](https://doi.org/10.5281/zenodo.19307989)  
**File:** [`papers/paper3/paper3_v1.pdf`](papers/paper3/)  
**Scripts:** [`code/paper3/`](code/paper3/)

**Imports from Paper 1:** `H_local` divergence (R1) — motivational only  
**Imports from Paper 2:** `W(g*,g*)` identity (R2), `D ≈ 9.471` (R3) — context only  
**Formal core (§§2–5) is self-contained.**

**What it constructs:**  
Two finite-dimensional real Hilbert spaces connected by a linear map and
a self-adjoint loop operator:

```
H_str  = l^2(primes ≤ kappa)              dim = pi(kappa)
H_null = l^2({gamma_1, ..., gamma_N})     dim = N
Phi: H_str → H_null,   Phi(e_p) = a_p
                        a_p[k] = exp(-eps^2 * gamma_k^2 / 2) * sin(gamma_k * log p)
T = Phi* ∘ Phi : H_str → H_str            T_{pq} = G^un_{pq}
```

| Result | Status |
|--------|--------|
| `T` self-adjoint, positive semi-definite | **PROVED** |
| `G^norm` is sigma-invariant | **PROVED** |
| Algebraic identity `Delta = E_str - E_spec` | **PROVED** |
| `eta_orig ∈ [0.598, 0.704]` for (kappa=53, eps=0.05) | **NUMERICAL** |
| `eta_orig > 0` for canonical weights | **NUMERICAL** (analytically open) |
| `lambda_max(T_ren) ≈ 0.39 * pi(kappa)` [grows with kappa] | **NUMERICAL** |
| Conjecture: `lambda_max < 1` (universal) | **FALSIFIED** (numerically) |
| `lambda_j(T) = mu_j(T̃)` to machine precision | **PROVED** (algebraic) |

**Reproduce:**
```bash
# Core energy identity and eta_orig table
python code/paper3/eta_verification.py

# Dual loop operator T̃ = Phi*Phi*, eigenvector localization
python code/paper3/ttilde_analysis.py

# eta_orig convergence as kappa → ∞, lambda_max scaling
python code/paper3/eta_inf_analysis.py
```

**Expected outputs:**

| Script | Key result | Reference |
|--------|-----------|-----------|
| `eta_verification.py` | `eta_orig = 0.66926873` | Sprint AF: 0.66926893 ✓ |
| `ttilde_analysis.py` | `loc_min ≥ 0.45`, `r(mu_j,1/gamma) = 0.950` | Sprint T̃ ✓ |
| `eta_inf_analysis.py` | `eta > 0` all kappa, `lambda/pi(kappa) ≈ 0.39` | Sprint η∞ ✓ |

**Figures generated:** `fig4_eta_spectrum.png`, `fig5_ttilde_localization.png`, `fig6_mu_vs_gamma.png`

---

## How the Papers Connect

```
Paper 1                        Paper 2                        Paper 3
───────────────────────────    ───────────────────────────    ───────────────────────────
H_xi = H_local + H_dual        W(g*,g*) = Z(g*)              Two Hilbert spaces:
                                          - H_local + O(eps)   H_str (primes)
Lemma 2: H_local(1/2,k)  ──R1──>  Feeds H_local into         H_null (zeros)
          ~ 2*(log k)^2            Weil functional        ──R1──>  sigma=1/2 as
                                                                    distinguished origin
Lemma 3: H_local(sigma,k)      c_p^ren = f_p^{1/2}             of H_null [motivation]
          → finite for           D = 9.471 < inf
          sigma > 1/2       ──R3──>  Context for eta  ──R2──>  Weil analogy in §10.2
                                     normalization                [context only]
V_p(sigma) >= 0               Open: Z >= H_local ?            T = Phi* Phi (explicit,
[local positivity]             [spectral inequality]            not postulated)
                                                               eta_orig > 0 [open]
```

**Mathematical thread:**
The divergence `H_local(1/2, kappa) ~ 2*(log kappa)^2` (Paper 1) motivates
`sigma = 1/2` as the distinguished origin of `H_null` (Paper 3). The Weil
functional identity (Paper 2) provides the outer framework within which
the energy decomposition `Delta = E_str - E_spec` (Paper 3) lives as a
qualitative structural analogy. The precise normalization-level identification
between the Weil framework and the eta-framework remains an open problem.

---

## Normative Parameters

All Paper 3 results use these normative parameters unless stated otherwise:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `kappa` | 53 | prime cutoff (16 active primes) |
| `eps` | 0.05 | Gaussian damping |
| `N` | 100 | zero ordinates used |
| `sigma` | 0.5 | evaluation point |
| `c_p` | `sqrt(8)*log(p)/p` | canonical weight vector |

---

## Repository Structure

```
analysislab-nt/
├── README.md
├── requirements.txt
├── LICENSE                    MIT
│
├── papers/
│   ├── paper1/
│   │   └── curvature_note_v7.tex     LaTeX source, Paper 1 v7
│   ├── paper2/
│   │   └── weil_v4.tex               LaTeX source, Paper 2 v4
│   └── paper3/
│       ├── paper3_v1.tex             LaTeX source, Paper 3 v1
│       └── paper3_v1.pdf             Compiled PDF, 14 pages
│
├── code/
│   ├── paper1/
│   │   └── verify_v6.py    H_local divergence, sigma profile, figures
│   ├── paper2/
│   │   └── verify_paper2.py  D=9.471, sawtooth, bridge constant, figures
│   └── paper3/
│       ├── eta_verification.py   eta_orig energy identity (main result)
│       ├── ttilde_analysis.py    T̃ = Phi*Phi*, eigenvector localization
│       └── eta_inf_analysis.py   convergence kappa → ∞, lambda_max scaling
│
├── data/
│   ├── zeros_100.csv              First 100 zeta zero ordinates gamma_k
│   ├── zeros_200.csv              First 200 zeta zero ordinates gamma_k
│   └── results/
│       ├── eta_table_kappa53.csv  eta_orig(sigma) normative table
│       └── ttilde_spectrum.csv    T̃ eigenvalues and localization data
│
└── figures/
    ├── paper1/
    │   ├── fig1_H_local_divergence.png   Lemma 2 divergence
    │   └── fig2_sigma_profile.png         Phase boundary
    ├── paper2/
    │   ├── fig3_Weil_decomposition.png   D convergence, f_p weights
    │   └── fig4_sawtooth_Hren.png        Mertens sawtooth
    └── paper3/
        ├── fig4_eta_spectrum.png          eta_orig(sigma) profile
        ├── fig5_ttilde_localization.png   T̃ eigenvector localization
        └── fig6_mu_vs_gamma.png           mu_j vs gamma_{k(j)} correlation
```

---

## Setup and Run

```bash
# Install dependencies
pip install -r requirements.txt

# Run all verifications
python code/paper1/verify_v6.py
python code/paper2/verify_paper2.py
python code/paper3/eta_verification.py
python code/paper3/ttilde_analysis.py
python code/paper3/eta_inf_analysis.py
```

**Requirements:** Python 3.x, NumPy ≥ 1.24, mpmath ≥ 1.3,
matplotlib ≥ 3.5, scipy ≥ 1.9

---

## Open Problems (as of March 2026)

| Problem | Statement | Paper |
|---------|-----------|-------|
| **OP 6.1** | Prove `eta_orig > 0` analytically for canonical `c_p = sqrt(8)*log(p)/p` | Paper 3 |
| **OP 6.2** | Prove `<T_ren * c̃, c̃> < 1` for canonical weights (Rayleigh-quotient bound) | Paper 3 |
| **OP 6.3** | Analytical cancellation bound via Abel summation | Paper 3 |
| **OP 6.4** | Analytical explanation for numerical sigma-invariance of `B_max` | Paper 3 |
| **OP 6.5** | Normalization-level identification of `(E_str, E_spec)` with `(Z_ren, H_ren)` | Papers 2–3 |
| **Spectral** | Unconditional lower bound `Z(g*_{1/2,eps}) >= C*(log kappa)^2` | Paper 2 |

**Closed (falsified):** Conjecture §31 `lambda_max(D^{-1/2} T D^{-1/2}) < 1`
universally is **numerically false**: `lambda_max ≈ 0.39 * pi(kappa) → ∞`.
The weaker statement for canonical `c_p` remains open (OP 6.2).

---

## Citation

```bibtex
@misc{tehrani2026curvature,
  author    = {Tehrani, Ulrich},
  title     = {A Curvature Decomposition of the Explicit Formula
               for the Riemann Zeta Function},
  year      = {2026},
  doi       = {10.5281/zenodo.19025598},
  publisher = {Zenodo}
}

@misc{tehrani2026weil,
  author    = {Tehrani, Ulrich},
  title     = {From Local Curvature to the Weil Functional:
               An Explicit Construction},
  year      = {2026},
  doi       = {10.5281/zenodo.19106992},
  publisher = {Zenodo}
}

@misc{tehrani2026hilbert,
  author    = {Tehrani, Ulrich},
  title     = {A Finite-Cutoff Hilbert-Space Model
               for Prime--Zero Energy Structure},
  year      = {2026},
  publisher = {Zenodo}
}
```

---

*AnalysisLab_L1_L5 · March 2026 · MIT License*