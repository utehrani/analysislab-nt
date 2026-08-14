#!/usr/bin/env python3
# ess_gate.py v1.0 — effective gate size (Kish ESS) for the section-6 battery.
# Revision, Aug 2026: names the one-line derivation of n ~ 7.1 / 17.6.
# w_k = exp(-eps~^2 b~_k^2 / 2) on the 30 unfolded ordinates b~_k = k - 1/2.
import numpy as np
b = np.arange(30) + 0.5
for eps in (0.25, 0.10):
    w = np.exp(-eps**2 * b**2 / 2)
    print(f"eps~={eps}:  n_eff = (sum w)^2 / sum w^2 = {w.sum()**2/np.square(w).sum():.2f}")
