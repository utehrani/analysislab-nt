# AnalysisLab: Curvature, Energy and the Riemann Zeta Function

Four papers on ζ(s): curvature decomposition, Weil functional, a
Hilbert-space model for prime–zero energy structure, and the Tehrani
operator T̃ = ΦΦ*.

**Author:** Ulrich Tehrani  
**License:** MIT  
**DOIs:** [Paper 1](https://doi.org/10.5281/zenodo.19025598) ·
[Paper 2](https://doi.org/10.5281/zenodo.19106992) ·
[Paper 3](https://doi.org/10.5281/zenodo.19307989) ·
[Paper 4](https://doi.org/10.5281/zenodo.19364703)

---

## The Four Papers

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
convergence above.

**Reproduce:**
```bash
python code/paper1/verify_v6.py
```

---

### Paper 2 — From Local Curvature to the Weil Functional
**DOI:** [10.5281/zenodo.19106992](https://doi.org/10.5281/zenodo.19106992)  
**File:** [`papers/paper2/weil_v5.pdf`](papers/paper2/)  
**Scripts:** [`code/paper2/`](code/paper2/)

**Imports from Paper 1:** `H_local` divergence (R1)

**What it proves:**  
Explicit admissible test functions `g*_{sigma,eps}` for the Weil explicit
formula satisfy

```
W(g* * g̃*) = Z(g*) - H_local(sigma, kappa) + O(eps)
```

| Result | Status |
|--------|--------|
| `c_p^ren = f_p^{1/2} > 0` (renormalized weights) | **PROVED** |
| `D = sum_p (c_p^ren)^2 ≈ 9.471 < ∞` | **PROVED** (Convergence Lemma) |
| `W(g*,g*) = Z(g*) - H_local + O(eps)` | **PROVED** (Theorem 3.1) |
| Bridge constant `D/(2*pi) ≈ 1.507` | **NUMERICAL** |

**Reproduce:**
```bash
python code/paper2/verify_paper2.py
```

---

### Paper 3 — A Finite-Cutoff Hilbert-Space Model for Prime–Zero Energy Structure
**DOI:** [10.5281/zenodo.19307989](https://doi.org/10.5281/zenodo.19307989)  
**File:** [`papers/paper3/paper3_v2.pdf`](papers/paper3/)  
**Scripts:** [`code/paper3/`](code/paper3/)

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
| Algebraic identity `Delta = E_str - E_spec` | **PROVED** |
| `eta_orig ∈ [0.598, 0.704]` for (kappa=53, eps=0.05) | **NUMERICAL** |
| `eta_orig > 0` for canonical weights, kappa ≤ 1009 | **NUMERICAL** |
| `lambda_max(T_ren) ≈ 0.39 * pi(kappa)` [grows with kappa] | **NUMERICAL** |
| Conjecture: `lambda_max < 1` (universal) | **FALSIFIED** (numerically) |
| `lambda_j(T) = mu_j(T̃)` to machine precision | **PROVED** (algebraic) |

**Reproduce:**
```bash
python code/paper3/eta_verification.py
python code/paper3/ttilde_analysis.py
python code/paper3/eta_inf_analysis.py
```

---

### Paper 4 — A Dual Operator for Prime–Zero Coupling and a Conditional Proof of Energy Asymmetry (The Tehrani Operator)
**DOI:** [10.5281/zenodo.19364703](https://doi.org/10.5281/zenodo.19364703)
**File:** [`papers/paper4/paper4_v4.pdf`](papers/paper4/)  
**Scripts:** [`code/paper4/`](code/paper4/)

**What it introduces:**  
The Tehrani operator T̃ = ΦΦ*, the dual of T = Φ*Φ, acting on the
finite-dimensional zero space H_null. Built directly from prime resonance
vectors — no new postulate.

```
T̃ = Phi ∘ Phi* : H_null → H_null
T̃_{kl} = sum_{p≤kappa} (a_p)_k * (a_p)_l
```

| Result | Status |
|--------|--------|
| T̃ self-adjoint, positive semi-definite | **PROVED** |
| Spectral identity σ(T)\{0} = σ(T̃)\{0} | **PROVED** |
| T̃ is NOT a Hilbert–Pólya operator | **NUMERICAL** (negative result) |
| W₁ = C_T · T̃⁺ self-adjoint | **PROVED** |
| Lemma M3 (Abel Summation Principle) | **PROVED** |
| Theorem EXPLICIT: M_k(κ) = O(π(κ)/γ_k) | **PROVED** (PNT only, no RH) |
| μ_j ~ C_T/γ_{k(j)}, r₁ = 0.950 | **NUMERICAL** |
| Eigenvector localization 0.45–0.99 | **NUMERICAL** |
| η_orig > 0 for κ ≤ 1009 | **NUMERICAL** |
| Δ_Burst ≈ 4.81 > 0, κ-invariant | **NUMERICAL** |
| Theorem 5.3: η_orig(κ) → η∞ > 0 | **CONDITIONAL** (OF-EXPLICIT-1') |

**Two arithmetic constants (distinct — never confuse):**

```
C_η ≈ 0.39  — λ_max,ren / π(κ)  [κ-independent]
C_T ≈ 2π(κ) — OLS slope μ_j ~ C_T/γ_{k(j)}  [grows with κ]
```

**Open problem OF-EXPLICIT-1' — genuine open content (b):**  
Joint equidistribution of {(γ_{k₁}log p, γ_{k₂}log p) mod 2π}.  
Via Weyl's criterion: (1/π(κ)) Σ_{p≤κ} e^{inγ_k log p} → 0 as κ→∞.  
Concerns Re(s)=1, not Re(s)=½. Not RH.  
*(Note: part (a) ζ(1+inγ_k)≠0 is settled by Hadamard 1896 — corrected in v4)*

**Reproduce:**
```bash
python code/paper4/verify_paper4.py
```

---

## How the Papers Connect

```
Paper 1              Paper 2              Paper 3              Paper 4
─────────────────    ─────────────────    ─────────────────    ─────────────────
H_xi = H_local       W(g*,g*) =           T = Phi*Phi          T̃ = Phi Phi*
       + H_dual       Z(g*)-H_local        eta_orig > 0         Resonance operator
                      + O(eps)             [Numerical]          Spectral structure

H_local(1/2,k)  ──R1──>  Weil bridge  ──R2──>  Geometry      ──R3──>  Conditional
  ~ 2(log k)²                                   H_str,H_null           η_orig > 0
                                                                         → η∞ > 0
```

**Mathematical thread:**  
Local curvature divergence at σ=½ (Paper 1) motivates σ=½ as distinguished
origin in H_null (Paper 3). The Weil identity (Paper 2) provides the outer
framework. The Tehrani operator T̃ (Paper 4) encodes the prime-mediated
coupling between zero ordinates and yields a conditional analytic proof
that energy asymmetry persists in the limit κ→∞.

---

## Normative Parameters

All results use these normative parameters unless stated otherwise:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `kappa` | 53 | prime cutoff (16 active primes) |
| `eps` | 0.05 | Gaussian damping |
| `N` | 100 | zero ordinates used |
| `sigma` | 0.5 | evaluation point |
| `c_p` | `sqrt(V_p(1/2))` | canonical weight vector |

---

## Repository Structure

```
analysislab-nt/
├── README.md
├── requirements.txt
├── LICENSE                         MIT
│
├── papers/
│   ├── paper1/
│   │   └── curvature_note_v7.tex   LaTeX source, Paper 1 v7
│   ├── paper2/
│   │   └── weil_v5.tex             LaTeX source, Paper 2 v5
│   ├── paper3/
│   │   ├── paper3_v2.tex           LaTeX source, Paper 3 v2
│   │   └── paper3_v2.pdf           Compiled PDF
│   └── paper4/
│       ├── paper4_v4.tex           LaTeX source, Paper 4 v4
│       └── paper4_v4.pdf           Compiled PDF
│
├── code/
│   ├── paper1/
│   │   └── verify_v6.py            H_local divergence, sigma profile
│   ├── paper2/
│   │   └── verify_paper2.py        D=9.471, sawtooth, bridge constant
│   ├── paper3/
│   │   ├── eta_verification.py     η_orig energy identity (main result)
│   │   ├── ttilde_analysis.py      T̃ = ΦΦ*, eigenvector localization
│   │   └── eta_inf_analysis.py     convergence κ→∞, λ_max scaling
│   └── paper4/
│       └── verify_paper4.py        spectral identity, η_orig, HP test,
│                                   generates figures/paper4/fig_hp_main.png
│
├── data/
│   ├── zeros_100.csv               First 100 Riemann zeta zero ordinates γ_k
│   ├── zeros_200.csv               First 200 Riemann zeta zero ordinates γ_k
│   └── results/                    Script outputs (CSV, intermediate data)
│       ├── eta_table_kappa53.csv   η_orig(σ) normative table
│       └── ttilde_spectrum.csv     T̃ eigenvalues and localization data
│
└── figures/
    ├── paper1/
    │   ├── fig1_H_local_divergence.png   Lemma 2 divergence at σ=½
    │   └── fig2_sigma_profile.png        Phase boundary σ=½
    ├── paper2/
    │   ├── fig3_Weil_decomposition.png   D convergence, f_p weights
    │   └── fig4_sawtooth_Hren.png        Mertens sawtooth
    ├── paper3/
    │   ├── fig4_eta_spectrum.png          η_orig(σ) profile
    │   ├── fig5_ttilde_localization.png   T̃ eigenvector localization
    │   └── fig6_mu_vs_gamma.png           μ_j vs γ_{k(j)} correlation
    └── paper4/
        └── fig_hp_main.png               T̃ spectral structure and HP test
```

---

## Setup and Run

```bash
# Install dependencies
pip install -r requirements.txt

# Run all verifications (from repo root)
python code/paper1/verify_v6.py
python code/paper2/verify_paper2.py
python code/paper3/eta_verification.py
python code/paper3/ttilde_analysis.py
python code/paper3/eta_inf_analysis.py
python code/paper4/verify_paper4.py
```

**Notes:**
- `verify_paper4.py` accepts an optional argument for the number of zero ordinates:
  ```bash
  python code/paper4/verify_paper4.py        # N=100 (default, uses zeros_100.csv)
  python code/paper4/verify_paper4.py 200    # N=200 (uses zeros_200.csv)
  ```
  Loads from `data/zeros_100.csv` or `data/zeros_200.csv` automatically.
  Falls back to mpmath if the CSV is missing.
- All scripts write figures to `figures/paperN/` and data to `data/results/`.
- Run from the repository root so that relative paths resolve correctly.

**Requirements:** Python 3.x, NumPy ≥ 1.24, mpmath ≥ 1.3,
matplotlib ≥ 3.5, sympy ≥ 1.11

---

## Open Problems (as of April 2026)

| Problem | Statement | Paper |
|---------|-----------|-------|
| **OF-EXPLICIT-1'(b)** | Joint equidistribution of {(γ_{k₁}log p, γ_{k₂}log p) mod 2π} — the genuine open condition (Weyl criterion) | Paper 4 |
| **OP 6.2** | Prove η_orig > 0 analytically without equidistribution assumption | Paper 4 |
| **OF-HP-3** | Does e^{-ε·γ_k} damping (instead of Gaussian) yield r₂→1? | Paper 4 |
| **L6** | Transfer σ=½ → σ≠½ on Γ-curvature level (the RH step) | Series |

**Closed / Settled:**  
- **OF-EXPLICIT-1'(a)** `ζ(1+inγ_k) ≠ 0 for all n≥1, k≤N`: **SETTLED** —
  follows directly from Hadamard–de la Vallée Poussin (1896),
  which establishes ζ(1+it) ≠ 0 for all real t ≠ 0.
  The values t = nγ_k are real and positive; no further proof needed.
  *(Corrected in v4, April 2026)*  
- Conjecture §31 `λ_max(D^{-1/2}TD^{-1/2}) < 1` universally: **FALSIFIED**  
- HP-question for W₁ = C_T·T̃⁺: **CLOSED** (r₂→0.16, structural failure)

---

## Citation

**Paper 1:**  
Tehrani, U. (2026). A Curvature Decomposition of the Explicit Formula.
Zenodo. https://doi.org/10.5281/zenodo.19025598

**Paper 2:**  
Tehrani, U. (2026). From Local Curvature to the Weil Functional.
Zenodo. https://doi.org/10.5281/zenodo.19106992

**Paper 3:**  
Tehrani, U. (2026). A Finite-Cutoff Hilbert-Space Model
for Prime–Zero Energy Structure. Zenodo.
https://doi.org/10.5281/zenodo.19307989

**Paper 4:**  
Tehrani, U. (2026). A Dual Operator for Prime–Zero Coupling
and a Conditional Proof of Energy Asymmetry. Zenodo.
https://doi.org/10.5281/zenodo.19364703

---

*v4/v2/v5 (April 2026) — Paper 4 v4 + Paper 3 v2 + Paper 2 v5: Weil framing clarified · MIT License*
