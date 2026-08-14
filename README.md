# AnalysisLab: An Explicit Prime–Zero Coupling Operator

Nine papers on an operator built explicitly from prime phasors:
curvature decomposition, Weil functional, a Hilbert-space model for
prime–zero energy structure, the Tehrani operator T̃(σ) = Φ(σ)Φ(σ)*, a
spectral trace formula with smoothed zero sums, positive curvature at
the critical line, a conditional selection principle identifying σ = ½
as a near-minimum, an unconditional proof that the full Guinand–Weil
second-moment bias is negative with certified off-line robustness, and
finally a classification of what this coupling can and cannot hear —
including a proved symmetry barrier that bounds the reach of the whole
construction.

The series is complete. It develops an operator-theoretic framework and
its limits; the Riemann Hypothesis is the historical motivation, not a
component of the results. No paper in the series claims a proof of it,
and none uses RH, GUE, Montgomery pair correlation or a Hilbert–Pólya
postulate as an input.

**Author:** Ulrich Tehrani  
**License:** code MIT; papers CC BY 4.0 (as published on Zenodo)  
**DOIs:** [Paper 1](https://doi.org/10.5281/zenodo.19025598) ·
[Paper 2](https://doi.org/10.5281/zenodo.19106992) ·
[Paper 3](https://doi.org/10.5281/zenodo.19307989) ·
[Paper 4](https://doi.org/10.5281/zenodo.19364703) ·
[Paper 5](https://doi.org/10.5281/zenodo.19508547) ·
[Paper 6](https://doi.org/10.5281/zenodo.19665790) ·
[Paper 7](https://doi.org/10.5281/zenodo.20440671) ·
[Paper 8](https://doi.org/10.5281/zenodo.20792123) ·
[Paper 9](https://doi.org/10.5281/zenodo.21899170)

---

## The Nine Papers

### Paper 1 — A Curvature Decomposition of the Explicit Formula
**DOI:** [10.5281/zenodo.19025598](https://doi.org/10.5281/zenodo.19025598)  
**File:** [`papers/paper1/`](papers/paper1/)  
**Scripts:** [`code/paper1/`](code/paper1/)

**Imports:** none — this is the entry point of the series.

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
| `H_local(1/2, kappa) ≍ (log kappa)^2 → ∞` | **PROVED** |
| Sharper ratio `H_local(1/2,kappa) / [2(log kappa)^2] → 1` | **NUMERICAL** (verified to kappa = 10⁶) |
| `H_local(sigma, kappa) → C(sigma) < ∞` for sigma > 1/2 | **PROVED** |

The critical line `sigma = 1/2` is the unique phase boundary: divergence below,
convergence above.

**One open problem:** whether the singular curvature kernel of this
formulation embeds into the admissible Weil test-function framework —
taken up by Paper 2.

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

**Three open problems:** the Weil-normalised embedding; the spectral
inequality against H_local(½,κ); off-diagonal control of the zero sum.

**Reproduce:**
```bash
python code/paper2/verify_paper2.py
```

---

### Paper 3 — A Finite-Cutoff Hilbert-Space Model for Prime–Zero Energy Structure
**DOI:** [10.5281/zenodo.19307989](https://doi.org/10.5281/zenodo.19307989)  
**File:** [`papers/paper3/`](papers/paper3/)  
**Scripts:** [`code/paper3/`](code/paper3/)

**Imports from Papers 1–2:** the prime-cutoff curvature decomposition and
the local weights V_p(σ), f_p.

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
| Spectral bound `B_max < 1` on the tested grid; `B^{‖c‖}_max < 0.382` | **NUMERICAL** |
| `lambda_j(T) = mu_j(T̃)` to machine precision | **PROVED** (algebraic) |

**Five open problems:** analytic proof of η_orig > 0 for the canonical
weight vector; the spectral bound for canonical weights; the analytical
cancellation bound; an analytic explanation of B_max; the renormalised
shell energy.

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

**Imports from Paper 3:** the coupling map Φ, the spaces H_str and H_null,
and the loop operator T = Φ*∘Φ, of which this paper studies the dual.

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
| No support for a Hilbert–Pólya reading of T̃: `r₂ = corr(ω_j, γ_{k(j)})` falls 0.50 → 0.16 | **NUMERICAL** (diagnostic; no obstruction theorem) |
| W₁ = C_T · T̃⁺ self-adjoint | **PROVED** |
| Δ = Δ_Burst + Δ_Cross + Δ_Stream (exact decomposition) | **PROVED** |
| Lemma M3 (Abel Summation Principle) | **PROVED** |
| M_k(κ) = O(π(κ)/γ_k) at fixed γ_k (PNT only, no RH) | **PROVED** |
| μ_j ≈ C_T/γ_{k(j)}, r₁ = 0.950 (finite-grid OLS) | **NUMERICAL** |
| `eta_ren > 0` on tested grid κ ∈ {23,53,101,199,503,1009} | **NUMERICAL** |
| `rho_max = 0.583 < 1` on tested grid | **NUMERICAL** |
| `eta_ren > 0` under `(E_rem)` and `Delta_Burst > 0` | **CONDITIONAL/NUMERICAL** |

**Seven open problems:** rank equality rank(T̃) = π(κ); a proof of
(E_rem); analytic positivity of η_ren without it; a circularity-free
spectral function; the arithmetic origin of C_η; asymptotics of C_T; and a
monotone bound from prime exponential-sum control.

**Reproduce:**
```bash
python code/paper4/verify_paper4.py
```

---

### Paper 5 — Spectral Trace Formula and Smoothed Zero Sums: A Prime–Zero Duality Framework
**DOI:** [10.5281/zenodo.19508547](https://doi.org/10.5281/zenodo.19508547)  
**File:** [`papers/paper5/`](papers/paper5/)  
**Scripts:** [`code/paper5/`](code/paper5/)

**Imports from Papers 3–4:** the coupling map Φ and the operator pair
T = Φ*∘Φ, T̃ = ΦΦ*, here extended to the σ-dependent family.

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

**Five open problems:** the Bias Conjecture; stationarity of O′(½); the
Weighted Bias Bridge; the weighted spectral measure; the Guinand–Weil
bridge.

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
| Truncation error: `\|R_{p,100}\|/\|Main_p\| ≤ 3×10⁻⁶¹` at the reference parameters | **NUMERICAL** (Thm 4.3, labelled numerical in the paper) |
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
Under two explicit hypotheses (subleading archimedean + proxy-transfer;
the proxy-transfer hypothesis is implied by RH and formally weaker in
content, the subleading one is not known to follow from RH), Paper 7
proves that the asymptotic constant-term
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

### Paper 8 — Unconditional Negativity of the Second-Moment Bias and Off-Line Robustness
**DOI:** [10.5281/zenodo.20792123](https://doi.org/10.5281/zenodo.20792123)  
**File:** [`papers/paper8/`](papers/paper8/)  
**Scripts:** [`code/paper8/`](code/paper8/)

**Imports from Papers 5–6:** trace formula Tr(T̃(σ)) = D_SEL − O(σ),
curvature identity O″(½) = −2B.

**What it proves:**  
Paper 8 addresses two open problems of Paper 7 unconditionally, by
bypassing the Weighted Bias Bridge rather than proving it. It
distinguishes the full Guinand–Weil second-moment object B_GW^∞ (over
all non-trivial zeros, free of any zero-location hypothesis) from the
on-line ordinate proxy B_line^∞ that carries the operator curvature
O″(½) = −2B_line^∞. Differentiating the Gaussian explicit-formula
identity in the smoothing parameter yields a second-moment identity in
which the negative diagonal prime self-interaction dominates; from it
B_GW^∞(53,ε) < 0 is proved for all 0 < ε ≤ 0.05 (closed-form for
ε ≤ 0.004, Arb-certified to the reference width ε = 0.05), together with
a uniform statement B_GW^∞(κ,ε) < 0 for every κ ≥ 202 and
0 < ε ≤ ε₀(κ) = 1/(4eκ√log κ). A β-uniform off-line magnitude estimate,
combined with the Platt–Trudgian verified height H₀ = 3·10¹², bounds the
off-line correction below the certified margin, transferring the sign to
B_line^∞ and giving O″(½) > 0 for ε_off ≤ ε ≤ 0.05. The archimedean
term is controlled by an explicit bound I₄ ≤ 3561.1 (four-fold
integration by parts), recertified in Arb ball arithmetic. No RH, GUE,
Montgomery, or Hilbert–Pólya input is used; the off-line bound rests on
a finite-height verification, not on RH.

| Result | Status |
|--------|--------|
| Second-moment explicit-formula identity | **PROVED** |
| Archimedean bound: `Γ_p^(2) = O(1)`, `I₄ ≤ 3561.1` | **PROVED / CERTIFIED** |
| Eigenterm extraction + cross-prime control | **PROVED** |
| Reference negativity: `B_GW^∞(53,ε) < 0`, `0 < ε ≤ 0.05` | **PROVED / CERTIFIED** |
| Uniform negativity: `B_GW^∞(κ,ε) < 0`, `κ ≥ 202` | **PROVED** |
| Off-line robustness (magnitude form) | **PROVED** |
| Operator curvature: `O″(½) = −2B_line^∞ > 0`, `ε_off ≤ ε ≤ 0.05` | **PROVED** |
| `B(53,0.05,100) = −19 342.5`, `O″(½) = +38 685.1` | **NUMERICAL** |
| Off-critical defect coefficient `C₂(γ) = −½G″(γ)` sign-indefinite | **PROVED (scope)** |

**Scope:** the negativity B_GW^∞ < 0 is a second-moment bias / curvature
statement, **not** an RH criterion: the off-critical continuation of the
curvature test function yields a sign-indefinite δ²-defect, so the
curvature sign alone does not furnish a Weil-positivity / RH criterion.

**Five open problems:** (1) Cancellation Hypothesis A** analytically;
(2) Weighted Bias Bridge (bypassed here, not proved); (3) exact
stationarity O′(½) = 0; (4) Weil-transfer operator and positivity;
(5) scaling of the negativity window.

**Reproduce:**
```bash
python code/paper8/cert_paper8.py     # Arb interval certificate
python code/paper8/verify_paper8.py   # numerical gate (30/0)
```

---

### Paper 9 — Scaling Laws, Class Invariance, and the Limits of Audibility of an Explicit Prime–Zero Coupling Operator
**DOI:** [10.5281/zenodo.21899170](https://doi.org/10.5281/zenodo.21899170)  
**File:** [`papers/paper9/`](papers/paper9/)  
**Scripts:** [`code/paper9/`](code/paper9/)

**Imports from Papers 3–8:** the coupling map Φ(σ), the trace formula
Tr(T̃(σ)) = D_SEL − O(σ), the curvature identity O″(½) = −2B, the
energy-asymmetry weights.

**What it proves:**  
Paper 9 closes the series by characterising the object rather than
extending it: what does this coupling hear, and what can it provably not
hear? The answer is a classification along three axes, and the order
matters. A symmetry barrier comes first: every symmetric functional of
the orbit of a zero under the functional equation is even in the signed
distance δ from the critical line, so the observables of this series are
blind to the orientation of that displacement — unconditionally, with no
regularity hypothesis. On the prime side, the weighted point
configuration separates from three control objects, kept apart
throughout: a smooth measure whose own zeta function is zero-free, its
rescaling to the exact prime-count mass, and a discrete equal-mass
quadrature world; the margins are certified in Arb, and for the discrete
pair no functional of the counting data with bounded-Lipschitz constant
below 35.51 can factor the separation. On the prime input side, both
scaling laws transfer verbatim across the Beurling density class,
carrying the same limiting function, the same certified separation and
the same leading profile — a within-class non-distinction, proved, not a
separation. On the ordinate side, a pre-registered unfolded
discrimination test over five worlds and five observables finds no
separation beyond the imposed unit-density skeleton.

| Result | Status |
|--------|--------|
| Kernel identity: `T_pq = ½[G(log(p/q)) − G(log(pq))]` | **PROVED** |
| Canonicity of the phasor kernel (trigonometric class) | **PROVED** |
| Non-convergence of `r(x)`; separation `Δ₀ ≥ 0.0840112` | **PROVED / CERTIFIED** |
| Divergence of the energy-asymmetry functional in mean | **PROVED** |
| Ordinate tail: limit parameters are `(κ,ε)`, not `(κ,N)` | **PROVED** |
| Class transfer over `D_δ`: same `R`, `Δ₀`, `F`, same band | **PROVED** |
| Sideband transfer for a smooth oscillating measure | **PROVED** (smooth family) |
| Curvature margins: `≥ 5343.90` (smooth), `≥ 5416.28` (discrete) | **CERTIFIED** |
| Bounded-Lipschitz distance `d_BL = 152.529455`, two-sided | **CERTIFIED** |
| Factorisation barrier: resolution lower bound `L ≥ 35.51` | **PROVED from certified inputs** |
| Functional-equation symmetry barrier | **PROVED** |
| Five-world unfolded battery: `max D_i/S_i = 0.630` vs. threshold 3 | **NUMERICAL** |

**Scope:** audible is used relatively throughout — a feature is called
audible if some observable of the stated class separates the world
carrying it from the stated comparison worlds. The certified prime-side
separation identifies the point configuration relative to those
controls; it does not by itself identify multiplicativity as the cause,
because the controls vary point geometry and multiplicative structure
together. The factorisation barrier covers one pair of weighted worlds.
No inference from any discriminator to the location of a zero is drawn
anywhere in the paper.

**Seven open problems:** pointwise finality of the divergence; the
minimal order break that destroys the class laws; the exceptional-set
measure under intermediate conditions; the family form of the
factorisation barrier; an identification design for the prime-side
separation; pair and correlation statistics on the unfolded axis; and
the combinatorial gap in the canonicity lemma.

**Reproduce:**
```bash
python code/paper9/verify_paper9.py                 # 277 checks
python code/paper9/certify_smooth_controls.py       # Arb: margins, d_BL, ordering
python code/paper9/cert_paper9_ratio.py             # Arb: two-point separation Δ₀
python code/paper9/unfolded_discrimination.py       # five-world battery
python code/paper9/function_side_reconstruction.py  # reconstruction experiment
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
                                                                 │
                                                                 ▼
                                                            Paper 8
                                                            ─────────────────
                                                            B_GW^∞ < 0 (PROVED)
                                                            Arb-certified
                                                            off-line robustness
                                                            O″(½) > 0 (transfer)
                                                            [not an RH criterion]
                                                                │
                                                                ▼
                                                            Paper 9
                                                        ─────────────────
                                                        symmetry barrier
                                                        class invariance
                                                        certified margins
                                                        [limits of audibility]
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
Paper 8 addresses two open problems of Paper 7 unconditionally: it
bypasses the Weighted Bias Bridge by proving the full Guinand–Weil
second-moment negativity B_GW^∞ < 0 directly (closed-form + Arb-certified),
and controls hypothetical off-critical zeros by a finite-height magnitude
estimate, transferring the sign to the operator curvature O″(½) > 0. The
curvature sign is explicitly **not** claimed as an RH criterion.
Paper 9 turns the question around and asks what the construction of
Papers 3–8 can resolve at all. It proves a symmetry barrier that bounds
the whole class of observables, certifies what the prime side does
separate, and shows that both scaling laws are invariant across the
Beurling density class. The series therefore ends not with a criterion
but with a map of the instrument: the questions are either answered or
provably unanswerable within it.

---

## Reference Parameters

All results use these reference parameters unless stated otherwise:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `kappa` | 53 | prime cutoff (16 active primes) |
| `eps` | 0.05 | Gaussian damping |
| `N` | 100 | zero ordinates used |
| `sigma` | 0.5 | evaluation point |
| `c_p^η` | `sqrt(V_p(½))` | η-framework weight (Papers 1–3; imported by Paper 9 for η_orig) |
| `c_p^ren` | `sqrt(f_p)` | renormalised Weil weight (Papers 4–6; imported by Paper 9 for η_ren) |
| `ε` window | `[0.04, 0.07]` | window on which Papers 8–9 state their certified margins |
| `ε̃` | 0.25 | unfolded scale of the Paper-9 battery; derived, not chosen, from the effective bandwidth at the reference parameters (robustness checked at 0.10) |

The two prime weights are numerically distinct and are never interchanged:
η_orig(½) = 0.69078176 uses `c_p^η`, η_ren(½) = 0.66926873 uses `c_p^ren`.
Papers 7–8 work with the curvature sums B and the Guinand–Weil object rather
than with a weight vector. Paper 9 uses both conventions explicitly and says at
each point which one is in force; the phasor argument carries no trace parameter
in the imported η-quantities, that is `sin(γ_k log p)`, which does not coincide
with `sin(σγ_k log p)` at σ = ½.

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
│   ├── paper7/
│   ├── paper8/
│   └── paper9/
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
│   ├── paper7/
│   │   └── verify_paper7.py        r_p^∞ convergence, O'/O'' sign checks,
│   │                               three-term decomposition, near-minimum σ*
│   ├── paper8/
│   │   ├── cert_paper8.py          Arb interval certificate, B_GW^∞ < 0
│   │   └── verify_paper8.py        second-moment identity, I₄ bound,
│   │                               eigenterm, off-line transfer (30/0)
│   └── paper9/
│       ├── verify_paper9.py        anchors, profile band, data contracts,
│       │                           textual anchors (277 checks)
│       ├── certify_smooth_controls.py  Arb: curvature margins, d_BL,
│       │                               certified ordering of |E_π|
│       ├── cert_paper9_ratio.py    Arb: two-point separation Δ₀
│       ├── unfolded_discrimination.py  five-world battery, completeness
│       ├── function_side_reconstruction.py  reconstruction experiment
│       ├── gen_zeros.py            ordinate regeneration protocol
│       └── ess_gate.py             effective gate size (Kish ESS)
│
├── data/
│   ├── zeros_100.csv               First 100 Riemann zeta zero ordinates γ_k (all papers, N=100)
│   ├── zeros_200.csv               First 200 Riemann zeta zero ordinates γ_k (Papers 5–7, N>100)
│   ├── zeros_650.csv               First 650 Riemann zeta zero ordinates γ_k (Papers 8–9, N>100)
│   ├── paper9_contract.json        Forbidden-pattern list read by verify_paper9.py
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
    ├── paper6/
    │   └── fig_paper6_main.png           r_p(ε) grid, sign-crossover,
    │                                     B_int vs ε, B vs κ scaling
    ├── paper8/
    │   └── fig_paper8_main.png           B_line(ε) proxy, B_line<0 with
    │                                     margins, O''(½)>0, |B_line| vs κ
    └── paper9/
        └── fig_paper9_main.png           certified ordering of |E_π|,
                                          profile band and exact envelope,
                                          curvature margins, D_i/S_i matrix

    (Paper 7 generates no figure; its verify script is purely numerical.
     The Paper 8 figure is produced by verify_paper8.py for the repository
     but is intentionally not embedded in the paper PDF.)
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
python code/paper8/cert_paper8.py      # Arb interval certificate (needs python-flint)
python code/paper8/verify_paper8.py
python code/paper9/verify_paper9.py
python code/paper9/certify_smooth_controls.py   # Arb (needs python-flint)
python code/paper9/cert_paper9_ratio.py         # Arb (needs python-flint)
python code/paper9/unfolded_discrimination.py
python code/paper9/function_side_reconstruction.py
python code/paper9/gen_zeros.py
python code/paper9/ess_gate.py
```

**Notes:**
- All scripts accept an optional argument for the number of zero ordinates:
  ```bash
  python code/paper6/verify_paper6.py        # N=100 (default, uses zeros_100.csv)
  python code/paper6/verify_paper6.py 200    # N=200 (uses zeros_200.csv)
  python code/paper8/verify_paper8.py        # N=100 (default, uses zeros_100.csv)
  python code/paper8/verify_paper8.py 650    # N=650 (uses zeros_650.csv)
  ```
  Papers 5–7 load `data/zeros_100.csv` or `data/zeros_200.csv`; Paper 8 loads
  `data/zeros_100.csv` or `data/zeros_650.csv`. All fall back to mpmath if the
  CSV is missing. The Paper-9 scripts do **not** fall back: each enforces the
  row count its statement rests on (100 for the normative layer, 650 for the
  extended list) and stops with a message naming the contract if it is not met.
  An explicitly given `--zeros`, `--tex` or `--contract` path is authoritative
  and never silently replaced by a copy next to the script.
- The `verify_paperN.py` scripts are numerical gates (mpmath; candidate values,
  not certificates). `cert_paper8.py` is different: it is a rigorous interval
  certificate that **proves** Theorem 6.2 in Arb ball arithmetic (midpoint–radius
  with directed rounding, proven Gaussian tail bounds). It requires
  `python-flint` (FLINT 3 / Arb layer) — listed in `requirements.txt` — and on
  success prints `CERTIFICATE VALID` together with the parameter-hash and the
  certificate source-hash recorded in the paper.
- All scripts write figures to `figures/paperN/` and data to `data/results/`.
  `verify_paper9.py` reads its forbidden-pattern list from
  `data/paper9_contract.json`; the figure block is wrapped so that a missing
  matplotlib prints a notice and leaves the check count unchanged.
- **Run from the repository root** so that relative paths resolve correctly.

**Requirements:** Python 3.x, NumPy ≥ 1.24, mpmath ≥ 1.3,
matplotlib ≥ 3.5, scipy ≥ 1.9, sympy ≥ 1.14, python-flint ≥ 0.8
(the last for the three Arb certificate scripts: `cert_paper8.py`,
`certify_smooth_controls.py`, `cert_paper9_ratio.py`)

---

## Open Problems (as of August 2026, after Paper 9)

A selection, ordered by the paper that states them; each paper's own list is
authoritative. None is used as a hypothesis anywhere in the series.

| Problem | Statement | Source |
|---------|-----------|--------|
| **Remainder control** | Prove (E_rem): \|Δ_Cross+Δ_Stream\| ≤ ρ·Δ_Burst analytically | Paper 4 |
| **Analytic positivity** | Prove η_ren > 0 without (E_rem) | Paper 4 |
| **rank(T̃) = π(κ)** | Requires linear independence of {a_p} in H_null; only ≤ min{N,π(κ)} proved | Paper 4 |
| **Bias Conjecture** | `Z_p^∞(ε) < 0` for each fixed prime p, all small ε — proved conditionally in Paper 7, bypassed in Paper 8 | Paper 5 |
| **η_∞ identity** | `η_∞ = 1 − m₁(∞)` algebraically | Paper 5 |
| **Uniformity at fixed ε** | `O''(½) > 0` for growing κ at the reference width; Paper 8 proves uniformity in the small-ε regime `ε ≤ ε₀(κ)`, not at ε = 0.05 | Paper 6 / Paper 8 OP 9.5 |
| **Weighted Bias Bridge** | `B_int^∞ < 0 ⇒ B_GW^∞ < 0` — Paper 8 bypasses the implication rather than proving it | Paper 7 OP 9.1 / Paper 8 OP 9.2 |
| **Cancellation Hypothesis (A\*\*)** | \|M_k(κ)\| ≤ C₀·P(κ)/(log γ_k)^A analytically; Vinogradov–Korobov class. Asks for less cancellation than RH, but is not known to follow from it | Paper 7 OP 9.2 / Paper 8 OP 9.1 |
| **Exact stationarity** | `O′(½) = 0`; at the reference parameters O′(½) = +2.4751 ≠ 0 | Paper 7 OP 9.5 / Paper 8 OP 9.3 |
| **Weil-transfer operator** | Finite-dimensional Weil form W^Weil_{κ,ε}; connection to the Lagarias framework | Paper 7 OP 9.6 / Paper 8 OP 9.4 |
| **Scaling of ε_int(κ)** | Asymptotic behaviour as κ→∞; is inf_κ ε_int(κ) > 0? | Paper 7 OP 9.4 |
| **Pointwise divergence** | Does `η_orig(κ) → −∞` without averaging? | Paper 9 OP 9.1 |
| **Minimal order break** | Smallest deviation from `ϑ_𝔅(x) ∼ x` that destroys the class laws | Paper 9 OP 9.2 |
| **Exceptional measure** | Any exceptional-set control under `D_δ`, and of what strength | Paper 9 OP 9.3 |
| **Family factorisation barrier** | Barrier uniform over a family of control worlds | Paper 9 OP 9.4 |
| **Pair statistics** | Do pair or higher correlations of unfolded ordinates separate arithmetic worlds under a pre-registered criterion? | Paper 9 OP 9.5 |
| **Identification design** | Controls matched in multiplicative structure while varying point geometry, or conversely | Paper 9 OP 9.6 |
| **Canonicity beyond finite spectra** | Do (H1)–(H2) force a single frequency for countable spectra? | Paper 9 OP 9.7 |

**Settled during the series.** Questions that were open at some point and are now
decided, with the status they carry:

- Asymptotic constant-term ratio `r_p^∞ = ½` for every prime: **CONDITIONAL**
  (Paper 7 Thm 3.1, under the subleading and proxy-transfer assumptions)
- Asymptotic pointwise bias `Z_{p,∞}^+ < 0`: **CONDITIONAL** (Paper 7 Cor 3.7,
  same assumptions)
- First-derivative cancellation bound: **CONDITIONAL** (Paper 7 Prop 7.1, under A\*\*)
- Direct second-moment negativity `B_GW^∞ < 0` (κ = 53, 0 < ε ≤ 0.05; uniform for
  κ ≥ 202 and ε ≤ ε₀(κ)): **PROVED / CERTIFIED** (Paper 8 Thm 6.2, 7.3)
- Off-critical control of `E_p^off` and transfer of the curvature sign — Paper 7's
  GW-bridge problem: **PROVED** (Paper 8 Thm 8.2 in magnitude form, via the
  Platt–Trudgian verified height; Cor 8.3 for the transfer)
- Hilbert–Pólya reading of W₁ = C_T·T̃⁺: **NUMERICAL** — the tested diagnostics give
  no support for it (r₂ → 0.16 as the cutoff grows); no structural obstruction
  theorem is proved (Paper 4)
- Class invariance of both scaling laws over the Beurling class `D_δ`: **PROVED**
  (Paper 9 Thm 4.2) — a within-class non-distinction, not a separation
- Separation of the prime configuration from the three stated controls:
  **CERTIFIED** (Paper 9 Thm 5.4, margins ≥ 5343.90 and ≥ 5416.28)
- Orientation of an off-critical displacement: **PROVED UNHEARABLE** for symmetric
  functionals of the zero orbit (Paper 9 Thm 6.1) — a boundary of the construction,
  not a gap in it
- Separation of arithmetic worlds on the unfolded ordinate axis by the registered
  battery: **NOT FOUND** (Paper 9 § 6, max D_i/S_i = 0.630 against a threshold of 3)
  — a measurement with declared scope, not a no-go theorem

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

**Paper 8:**  
Tehrani, U. (2026). Unconditional Negativity of the Second-Moment Bias
and Off-Line Robustness. Zenodo.
https://doi.org/10.5281/zenodo.20792123

**Paper 9:**  
Tehrani, U. (2026). Scaling Laws, Class Invariance, and the Limits of
Audibility of an Explicit Prime–Zero Coupling Operator. Zenodo.
https://doi.org/10.5281/zenodo.21899170

---

*Papers 1–9 · v6.0.0 · August 2026 · code MIT · papers CC BY 4.0*
