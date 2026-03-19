# AnalysisLab: Curvature and the Weil Functional

A numerical and analytical study of the curvature decomposition of the
explicit formula for the Riemann zeta function, and its connection to
the Weil positivity framework.

---

## The core decomposition

The second logarithmic derivative of the completed zeta function decomposes as

```
H_ξ(σ, t) = H_local(σ, κ) + H_dual(σ, t, κ)
```

where `H_local` collects local prime contributions and `H_dual` encodes
the spectral remainder. This note studies the structure of `H_local`:
its positivity, divergence rate, and connection to the Weil functional.

To the best of our knowledge, this curvature formulation does not appear
explicitly in the literature (J.-F. Burnol, personal communication).

---

## Verify everything yourself

All numerical claims in the papers can be independently verified.

**Paper 1 (curvature note v6):**

```bash
python code/paper1/verify_v6.py
```

Expected output:

```
=== Verification: curvature_note_v6 ===

Lemma 1: V_p(σ) ≥ 0
  p=2, σ=0.5: V_p = 5.545  ✓
  p=3, σ=0.5: V_p = 3.625  ✓
  p=5, σ=0.5: V_p = 3.213  ✓

Lemma 2: H_local(½, κ) ~ C·(log κ)²
  κ=100:   H = 43.98,  C = 2.074  ✓
  κ=1000:  H = 95.86,  C = 2.009  ✓
  κ=10000: H = 169.47, C = 1.998  ✓

Lemma 3: H_local(σ > ½, κ) converges
  σ=0.7: H(1000) = 18.49, H(10000) = 21.46  ✓

Numerical illustration:
  H_local(0.3, 1300) = 914.76  ≈ 915  ✓
  C ≈ 2  ✓

=== All checks passed ✓ ===
```

**Reproduce all figures:**

```bash
python code/paper1/plot_paper1.py   # generates figures/paper1/
python code/paper2/plot_paper2.py   # generates figures/paper2/
```

---

## Repository structure

```
analysislab-nt/
│
├── papers/
│   ├── paper1/
│   │   ├── curvature_note_v6.tex
│   │   └── curvature_note_v6.pdf
│   └── paper2/
│       ├── weil_decomposition_v1.tex    (in preparation)
│       └── weil_decomposition_v1.pdf    (in preparation)
│
├── code/
│   ├── paper1/
│   │   ├── verify_v6.py        ← verify all claims in paper 1
│   │   └── plot_paper1.py      ← reproduce all figures in paper 1
│   └── paper2/
│       ├── verify_paper2.py    ← verify all claims in paper 2
│       └── plot_paper2.py      ← reproduce all figures in paper 2
│
├── data/
│   ├── zeros_100.csv           ← first 100 nontrivial zeros of ζ(s)
│   └── results/
│       ├── paper1_table.csv    ← H_local values, C constants
│       └── paper2_table.csv    ← c_ren, D, Z vs H_local
│
└── figures/
    ├── paper1/
    │   ├── fig1_H_local_divergence.png
    │   ├── fig2_sigma_profile.png
    │   └── fig3_Vp_profile.png
    └── paper2/
        ├── fig4_c_ren_profile.png
        └── fig5_Z_vs_H_local.png
```

---

## Requirements

```bash
pip install mpmath sympy matplotlib numpy
```

Python 3.8 or later. No other dependencies.

---

## Papers

**Paper 1:** *A Curvature Decomposition of the Explicit Formula for the
Riemann Zeta Function*
— Ulrich Tehrani, March 2026
— [Zenodo DOI: 10.5281/zenodo.19025598. ](https://10.5281/zenodo.19025598)

**Paper 2:** *Renormalized Weights for the Weil Explicit Formula:
Construction and Numerical Evidence*
— Ulrich Tehrani, March 2026 (in preparation)

---

## Open question

The central open problem studied here is whether an unconditional
lower bound

```
Z(g*_{½,κ}) ≥ H_local(½, κ) ~ 2·(log κ)²
```

can be established for the zero-sum `Z` of the admissible test function
`g*`, without assuming RH. See Paper 1, §7 and Paper 2, §8.

---

## License

MIT License — see LICENSE file.

---

*utehrani · March 2026*