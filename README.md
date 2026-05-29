# AnalysisLab: Curvature, Energy and the Riemann Zeta Function

Seven papers on ζ(s): curvature decomposition, Weil functional, a
Hilbert-space model for prime–zero energy structure, the Tehrani
operator T̃(σ) = Φ(σ)Φ(σ)*, a spectral trace formula with smoothed
zero sums, positive curvature at the critical line, and a conditional
selection principle identifying σ = ½ as a near-minimum.

**Author:** Ulrich Tehrani  
**License:** MIT  
**DOIs:** [Paper 1](https://doi.org/10.5281/zenodo.19025598) ·
[Paper 2](https://doi.org/10.5281/zenodo.19106992) ·
[Paper 3](https://doi.org/10.5281/zenodo.19307989) ·
[Paper 4](https://doi.org/10.5281/zenodo.19364703) ·
[Paper 5](https://doi.org/10.5281/zenodo.19508547) ·
[Paper 6](https://doi.org/10.5281/zenodo.19665790) ·
[Paper 7](https://doi.org/10.5281/zenodo.20440671)

---

## The Seven Papers

### Paper 1 — A Curvature Decomposition of the Explicit Formula
**DOI:** [10.5281/zenodo.19025598](https://doi.org/10.5281/zenodo.19025598)  
**File:** [`papers/paper1/`](papers/paper1/)  
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
| `V_p(sigma) >= 0` for all primes p, sigma > 0 | **PROVED** |
| `H_local(1/2, kappa) ~ 2*(log kappa)^2 → ∞` | **PROVED** |
| `H_local(sigma, kappa) → C(sigma) < ∞` for sigma > 1/2 | **PROVED** |

The critical line `sigma = 1/2` is the unique phase boundary: divergence below,
convergence above.

**Reproduce:**
```bash
python code/paper1/verify_paper1.py
```

---

### Paper 2 — From Local Curvature to the Weil Functional
**DOI:** [10.5281/zenodo.19106992](https://doi.org/10.5281/zenodo.19106992)  
**File:** [`papers/paper2/`](papers/paper2/)  
**Scripts:** [`code/paper2/`](code/paper2/)

**Imports from Paper 1:** `H_local` divergence

**What it proves:**  
Explicit admissible test functions `g*_{sigma,eps}` for the Weil explicit
formula, with renormalized prime weights and convergent diagonal energy.

| Result | Status |
|--------|--------|
| Admissibility of `g*` in `S_ad` | **PROVED** |
| `c_p^ren = f_p^{1/2} > 0` (renormalized weights) | **PROVED** |
| `D = sum_p (c_p^ren)^2` converges | **PROVED** |
| Prime-side localisation `K_ε` targets `V_p(σ)` exactly (Observation 3.2) | **PROVED** |
| `D ≈ 9.470` at reference parameters | **NUMERICAL** |
| Finite-grid stability of `Z - H_local` | **NUMERICAL** |

**Reproduce:**
```bash
python code/paper2/verify_paper2.py
```

---

### Paper 3 — A Finite-Cutoff Hilbert-Space Model for Prime–Zero Energy Structure
**DOI:** [10.5281/zenodo.19307989](https://doi.org/10.5281/zenodo.19307989)  
**File:** [`papers/paper3/`](papers/paper3/)  
**Scripts:** [`code/paper3/`](code/paper3/)

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
| `eta_orig ∈ [0.650, 0.700]` for (kappa=53, eps=0.05), `c_p^eta = sqrt(V_p(1/2))` | **NUMERICAL** |
| `eta_orig(1/2) = 0.69078` at reference parameters | **NUMERICAL** |
| `eta_orig > 0` for kappa ∈ {23,53,101,199,503,1009} | **NUMERICAL** |
| `lambda_max(T_ren) ≈ 0.39 * pi(kappa)` [grows with kappa] | **NUMERICAL** |
| Conjecture: `lambda_max < 1` (universal) | **FALSIFIED** (numerically) |
| `lambda_j(T) = mu_j(T̃)` to machine precision | **PROVED** (algebraic) |

**Reproduce:**
```bash
python code/paper3/verify_paper3.py
python code/paper3/ttilde_analysis.py
```

---

### Paper 4 — A Dual Operator for Prime–Zero Coupling and a Conditional Proof of Energy Asymmetry
**DOI:** [10.5281/zenodo.19364703](https://doi.org/10.5281/zenodo.19364703)  
**File:** [`papers/paper4/`](papers/paper4/)  
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
| `rank(T̃) ≤ min{N, π(κ)}` | **PROVED** |
| T̃ is NOT a Hilbert–Pólya operator | **NUMERICAL** (negative result) |
| W₁ = C_T · T̃⁺ self-adjoint | **PROVED** |
| Δ = Δ_Burst + Δ_Cross + Δ_Stream (exact decomposition) | **PROVED** |
| Lemma M3 (Abel Summation Principle) | **PROVED** |
| M_k(κ) = O(π(κ)/γ_k) at fixed γ_k (PNT only, no RH) | **PROVED** |
| μ_j ≈ C_T/γ_{k(j)}, r₁ = 0.950 (finite-grid OLS) | **NUMERICAL** |
| `eta_ren > 0` on tested grid κ ∈ {23,53,101,199,503,1009} | **NUMERICAL** |
| `rho_max = 0.583 < 1` on tested grid | **NUMERICAL** |
| `eta_ren > 0` under `(E_rem)` and `Delta_Burst > 0` | **CONDITIONAL/NUMERICAL** |

**Reproduce:**
```bash
python code/paper4/verify_paper4.py
```

---

### Paper 5 — Spectral Trace Formula and Smoothed Zero Sums: A Prime–Zero Duality Framework
**DOI:** [10.5281/zenodo.19508547](https://doi.org/10.5281/zenodo.19508547)  
**File:** [`papers/paper5/`](papers/paper5/)  
**Scripts:** [`code/paper5/`](code/paper5/)

**What it introduces:**  
A σ-dependent operator family T̃(σ) = Φ(σ)∘Φ(σ)* and an exact algebraic
trace formula. The smoothed zero-sum theory exposes a structural negative
bias at σ=½ via the Bias Conjecture.

```
T̃(σ) = Phi(σ) ∘ Phi(σ)* : H_null → H_null
Phi(σ)_{k,p} = exp(-eps²*gamma_k²/2) * sin(sigma * gamma_k * log p)
Tr(T̃(σ)) = D_SEL − O(σ)   [proved algebraically]
```

| Result | Status |
|--------|--------|
| Trace formula: `Tr(T̃(σ)) = D_SEL − O(σ)` | **PROVED** (algebraic) |
| `D_SEL = (1/2) · A(ε,N) · π(κ) = 10.985` | **PROVED** |
| Decomposition: `B = Σ_p (log p)² Re(Z_p^{(2)})` | **PROVED** |
| Structural Reduction: dominant term of Z̃_p^GW(ε) formally identified | **PROVED** |
| Sign transfer: `Re(Z̃_p^GW) < 0 → B_int < 0` | **CONDITIONAL** |
| `B = −19342.5 < 0` (κ=53, ε=0.05, N=100) | **NUMERICAL** |
| `Re(Z_p) < 0` for 14 of 16 primes p ≤ 53 at ε=0.05 | **NUMERICAL** |
| `η_ren(κ=53) = 0.66927`, η_∞ ≈ 0.81 | **NUMERICAL** |
| Three spectral signatures at σ=½ | **NUMERICAL** |
| Bias Conjecture: `Z_p^∞(ε) < 0` for each fixed prime p | **OPEN** |

**Reproduce:**
```bash
python code/paper5/verify_paper5.py
```

---

### Paper 6 — Positive Curvature of the Spectral Trace at the Critical Line
**DOI:** [10.5281/zenodo.19665790](https://doi.org/10.5281/zenodo.19665790)  
**File:** [`papers/paper6/`](papers/paper6/)  
**Scripts:** [`code/paper6/`](code/paper6/)

**Imports from Paper 5:** trace formula Tr(T̃(σ)) = D_SEL − O(σ),
B-decomposition, B = −19342.5 at reference parameters.

**What it proves:**  
Starting from the trace formula of Paper 5, Paper 6 establishes a
four-term Guinand–Weil decomposition of the smoothed zero sum Z̃_p^GW(ε)
with quantitative control of each component. The algebraic identity
O″(½) = −2B connects the second derivative of the oscillatory trace
component to the direct curvature sum. At the reference parameters,
O″(½) > 0 follows from B = −19342.5 < 0.

| Result | Status |
|--------|--------|
| Main term negativity: `Main_p(ε) < 0` for all p, all ε>0 | **PROVED** |
| Other-prime error: `Err_other ≤ 0` | **PROVED** |
| Truncation error: `|R_{p,100}|/|Main_p| ≤ 3×10⁻⁶¹` | **NUMERICAL** |
| Curvature–bias identity: `O''(½) = −2B` | **PROVED** |
| Gamma term subleading: `\|Γ_p(ε)\|/\|Main_p(ε)\| ~ ε` | **NUMERICAL** |
| Proxy ratio consistent with `r_p^∞ ≤ ½` as ε → 0 | **NUMERICAL** |
| Sign-crossover localised in (0.020, 0.025) | **NUMERICAL** |
| Integrated bias: `B_int^+(0.05,100) = −42.21 < 0` | **NUMERICAL** |
| Positive curvature: `O''(½) = +38685 > 0` | **NUMERICAL** |
| Strict local minimum of O at σ=½ (under stationarity) | **CONDITIONAL** |

**Five open problems** frame the completion path: (1) analytic first-derivative
cancellation estimate |S_κ(γ)| ≤ C₀·P(κ)/(γ(log γ)^A) — the trivial bound is
explicit, the nontrivial logarithmic saving is open; (2) asymptotic pointwise bias
and r_p^∞ = ½ analytically via GW bridge; (3) uniformity of positive curvature;
(4) integrated-to-direct transfer; (5) bridge to Weil positivity.

**Reproduce:**
```bash
python code/paper6/verify_paper6.py
```

---

### Paper 7 — Conditional Stationarity and Positive Curvature of the Spectral Trace at the Critical Line
**DOI:** [10.5281/zenodo.20440671](https://doi.org/10.5281/zenodo.20440671)  
**File:** [`papers/paper7/`](papers/paper7/)  
**Scripts:** [`code/paper7/`](code/paper7/)

**Imports from Papers 5–6:** trace formula Tr(T̃(σ)) = D_SEL − O(σ),
curvature identity O″(½) = −2B, Guinand–Weil decomposition.

**What it proves:**  
Under two explicit hypotheses (subleading archimedean + proxy-transfer,
both weaker than RH), Paper 7 proves that the asymptotic constant-term
ratio r_p^∞ = ½ for every prime p, and that the integrated bias
B_int^∞(κ,ε) < 0 for all sufficiently small ε. An unconditional
three-term decomposition of O′(½) is introduced, and antisymmetry
A(½+δ) = −A(½−δ) is proved. Under an additional logarithmic-saving
hypothesis on prime exponential sums (Cancellation Hypothesis A**),
a conditional derivative bound on |O′(½)| is established. The full
conditional selection — strict local minimum of O at σ = ½ — requires
B < 0 and exact stationarity O′(½) = 0 as additional inputs.

| Result | Status |
|--------|--------|
| Asymptotic ratio: `r_p^∞ = ½` for every prime p | **CONDITIONAL** (SAA+PT) |
| Integrated bias: `B_int^∞ < 0` for small ε | **CONDITIONAL** (SAA+PT) |
| Sign: `Z_{p,∞}^+(ε) < 0` for small ε, all p | **CONDITIONAL** (SAA+PT) |
| Three-term decomposition of O′(½) | **DEFINITION** |
| Antisymmetry: `A(½+δ) = −A(½−δ)` | **PROVED** |
| Derivative bound: `|O′(½)| ≤ C·W·P` | **CONDITIONAL** (A**) |
| Local minimum at σ=½ | **CONDITIONAL** (B<0 + O′=0) |
| O′(½) = +2.4751 ≠ 0, σ* ≈ 0.4999 | **NUMERICAL** |
| Cancellation ratio: `|O′|/(W·P) ≈ 0.00194` | **NUMERICAL** |

**Six open problems:** (1) Weighted Bias Bridge; (2) Cancellation
Hypothesis analytically; (3) GW bridge with off-critical control;
(4) scaling of ε_int(κ); (5) near-minimum behaviour; (6) Weil-transfer
operator and positivity.

**Reproduce:**
```bash
python code/paper7/verify_paper7.py
```

---

## How the Papers Connect

```
Paper 1              Paper 2              Paper 3              Paper 4
─────────────────    ─────────────────    ─────────────────    ─────────────────
H_xi = H_local       W(g*,g*) =           T = Phi*Phi          T̃ = Phi Phi*
       + H_dual       Z(g*)-H_local        eta_orig > 0         Resonance operator
                      + O(eps)             [Numerical]          Spectral structure

H_local(1/2,k)  ────>  Weil bridge  ────>  Geometry      ────>  Conditional
  ~ 2(log k)²                                   H_str,H_null           η_orig > 0
                                                                         → η∞ > 0
                                                                              │
                                                                              ▼
                                                                         Paper 5
                                                                    ─────────────────
                                                                    T̃(σ) family
                                                                    Tr = D_SEL − O(σ)
                                                                    Bias Conjecture
                                                                         │
                                                                         ▼
                                                                    Paper 6
                                                                ─────────────────
                                                                O''(½) = −2B
                                                                O''(½) > 0 at ref.
                                                                [5 open problems]
                                                                     │
                                                                     ▼
                                                                Paper 7
                                                            ─────────────────
                                                            r_p^∞ = ½ (COND.)
                                                            O'(½) = +2.4751
                                                            Selection hierarchy
                                                            [6 open problems]
```

**Mathematical thread:**  
Local curvature divergence at σ=½ (Paper 1) motivates σ=½ as distinguished
origin in H_null (Paper 3). The Weil identity (Paper 2) provides the outer
framework. The Tehrani operator T̃ (Paper 4) encodes the prime-mediated
coupling between zero ordinates. Paper 5 introduces the σ-dependent family
T̃(σ), proves the exact trace formula Tr(T̃(σ)) = D_SEL − O(σ), and opens
the smoothed zero-sum route toward the Bias Conjecture.
Paper 6 uses the trace formula of Paper 5 to establish the positive-curvature
statement O''(½) > 0 at reference parameters (numerical) via the algebraic
identity O''(½) = -2B (proved).
Paper 7 addresses two open problems of Paper 6: it proves r_p^∞ = ½
conditionally and establishes a three-layer selection hierarchy
(SAA+PT → A** → B<0+O′=0) for σ = ½ as a local minimum.

---

## Reference Parameters

All results use these reference parameters unless stated otherwise:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `kappa` | 53 | prime cutoff (16 active primes) |
| `eps` | 0.05 | Gaussian damping |
| `N` | 100 | zero ordinates used |
| `sigma` | 0.5 | evaluation point |
| `c_p` | `sqrt(V_p(½))` (Papers 1–3) / `sqrt(f_p)` (Papers 4–6) | weight convention |

---

## Repository Structure

```
analysislab-nt/
├── README.md
├── requirements.txt
├── LICENSE                         MIT
│
├── papers/
│   ├── paper1/                     LaTeX source + PDF
│   ├── paper2/
│   ├── paper3/
│   ├── paper4/
│   ├── paper5/
│   ├── paper6/
│   └── paper7/
│
├── code/
│   ├── paper1/
│   │   └── verify_paper1.py            H_local divergence, sigma profile
│   ├── paper2/
│   │   └── verify_paper2.py        D=9.470, sawtooth, diagonal energy
│   ├── paper3/
│   │   ├── verify_paper3.py        η_orig, E_str, B-diagnostics, c_p^eta
│   │   └── ttilde_analysis.py      T̃ = ΦΦ*, eigenvector localization
│   ├── paper4/
│   │   └── verify_paper4.py        spectral identity, η_ren, HP test
│   ├── paper5/
│   │   └── verify_paper5.py        trace formula, B-decomposition,
│   │                               η_ren, Re(Z_p), 3 signatures
│   ├── paper6/
│   │   └── verify_paper6.py        curvature identity, B, B_int^+,
│   │                                   r_p ratio, sign-crossover localisation
│   └── paper7/
│       └── verify_paper7.py        r_p^∞ convergence, O'/O'' sign checks,
│                                   three-term decomposition, near-minimum σ*
│
├── data/
│   ├── zeros_100.csv               First 100 Riemann zeta zero ordinates γ_k
│   ├── zeros_200.csv               First 200 Riemann zeta zero ordinates γ_k
│   └── results/                    Script outputs (CSV, intermediate data)
│       └── ttilde_spectrum.csv     T̃ eigenvalues and localization data
│
└── figures/
    ├── paper1/
    │   ├── fig1_H_local_divergence.png   H_local divergence at σ=½
    │   └── fig2_sigma_profile.png        Phase boundary σ=½
    ├── paper2/
    │   └── fig3_Weil_decomposition.png   D convergence, f_p weights
    ├── paper3/
    │   ├── fig_paper3_main.png            η_orig(σ) profile, E_str, B-diagnostics
    │   ├── fig5_ttilde_localization.png   T̃ eigenvector localization
    │   └── fig6_mu_vs_gamma.png           μ_j vs γ_{k(j)} correlation
    ├── paper4/
    │   └── fig_hp_main.png               T̃ spectral structure and HP test
    ├── paper5/
    │   └── fig_paper5_main.png           trace formula, Re(Z_p),
    │                                      η_ren convergence, 3 signatures
    └── paper6/
        └── fig_paper6_main.png           r_p(ε) grid, sign-crossover,
                                          B_int vs ε, B vs κ scaling
```

---

## Setup and Run

```bash
# Install dependencies
pip install -r requirements.txt

# Run all verifications (from repo root)
python code/paper1/verify_paper1.py
python code/paper2/verify_paper2.py
python code/paper3/verify_paper3.py
python code/paper3/ttilde_analysis.py
python code/paper4/verify_paper4.py
python code/paper5/verify_paper5.py
python code/paper6/verify_paper6.py
python code/paper7/verify_paper7.py
```

**Notes:**
- All scripts accept an optional argument for the number of zero ordinates:
  ```bash
  python code/paper6/verify_paper6.py        # N=100 (default, uses zeros_100.csv)
  python code/paper6/verify_paper6.py 200    # N=200 (uses zeros_200.csv)
  ```
  Loads from `data/zeros_100.csv` or `data/zeros_200.csv` automatically.
  Falls back to mpmath if the CSV is missing.
- All scripts write figures to `figures/paperN/` and data to `data/results/`.
- **Run from the repository root** so that relative paths resolve correctly.

**Requirements:** Python 3.x, NumPy ≥ 1.24, mpmath ≥ 1.3,
matplotlib ≥ 3.5, scipy ≥ 1.9, sympy ≥ 1.14

---

## Open Problems (as of May 2026)

| Problem | Statement | Paper |
|---------|-----------|-------|
| **Bias Conjecture** | `Z_p^∞(ε) < 0` for each fixed prime p, all small ε | Paper 5 |
| **η_∞ identity** | `η_∞ = 1 − m₁(∞)` algebraically | Paper 5 |
| **Remainder control** | Prove (E_rem): |Δ_Cross+Δ_Stream| ≤ ρ·Δ_Burst analytically | Paper 4 |
| **Analytic positivity** | Prove η_ren > 0 without (E_rem) | Paper 4 |
| **Weighted Bias Bridge** | B_int^∞ < 0 → B < 0 via γ_k²-weighted transfer + finite-N truncation | Paper 7 OP 9.1 |
| **Cancellation Hypothesis (A**)** | \|M_k(κ)\| ≤ C₀·P(κ)/(log γ_k)^A analytically; Vinogradov–Korobov class | Paper 7 OP 9.2 |
| **GW bridge + off-critical control** | Bound E_p^off without RH; unconditional proxy transfer | Paper 7 OP 9.3 |
| **Scaling of ε_int(κ)** | Asymptotic behaviour as κ→∞; is inf_κ ε_int(κ) > 0? | Paper 7 OP 9.4 |
| **Near-minimum behaviour** | Under what stationarity condition does σ=½ become exact minimum? | Paper 7 OP 9.5 |
| **Weil-transfer operator** | Finite-dim Weil form W^Weil_{κ,ε}; connection to Lagarias framework | Paper 7 OP 9.6 |
| **Uniformity** | `O''(½) > 0` beyond reference parameters | Paper 6 |
| **rank(T̃) = π(κ)** | Requires linear independence of {a_p} in H_null; only ≤ min{N,π(κ)} proved | Paper 4 |
| **Remainder control** | Prove (E_rem): \|Δ_Cross+Δ_Stream\| ≤ ρ·Δ_Burst analytically | Paper 4 |
| **Analytic positivity** | Prove η_ren > 0 without (E_rem) | Paper 4 |
| **Bias Conjecture** | `Z_p^∞(ε) < 0` for each fixed prime p, all small ε | Paper 5 |
| **η_∞ identity** | `η_∞ = 1 − m₁(∞)` algebraically | Paper 5 |

**Closed / Settled:**
- `ζ(1+iγ_k) ≠ 0`: **SETTLED** (Hadamard 1896)
- `λ_max < 1` universally: **FALSIFIED**
- HP-question for W₁ = C_T·T̃⁺: **CLOSED** (r₂ → 0.16)
- `r_p^∞ = ½` for every prime: **CONDITIONAL** (Paper 7 Thm 3.1 under SAA+PT)
- First-derivative cancellation bound: **CONDITIONAL** (Paper 7 Prop 7.1 under A**)
- Asymptotic pointwise bias `Z_{p,∞}^+ < 0`: **CONDITIONAL** (Paper 7 Cor 3.7 under SAA+PT)

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

**Paper 5:**  
Tehrani, U. (2026). Spectral Trace Formula and Smoothed Zero Sums:
A Prime–Zero Duality Framework. Zenodo.
https://doi.org/10.5281/zenodo.19508547

**Paper 6:**  
Tehrani, U. (2026). Positive Curvature of the Spectral Trace at
the Critical Line. Zenodo.
https://doi.org/10.5281/zenodo.19665790

**Paper 7:**  
Tehrani, U. (2026). Conditional Stationarity and Positive Curvature
of the Spectral Trace at the Critical Line. Zenodo.
https://doi.org/10.5281/zenodo.20440671

---

*Papers 1–7 · v4.0.0 · May 2026 · MIT License*
