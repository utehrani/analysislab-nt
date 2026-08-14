#!/usr/bin/env python3
# certify_smooth_controls.py v2.0 — smooth/discrete margin certificate (AnalysisLab_L1_L5)
# One script reproduces: anchor table, certified zero enclosures, T-DISC
# interval certificate (Arb / python-flint), T-INV explicit constants, d_BL
# (LP + rigorous dual upper bound), scatter/mass-wedge balance.
#
# Changelog v2.0 (revision, Aug 2026 — external review; error record E-6):
#   - ordinate chain: acb.zeta_zeros(1,100) isolation at working precision,
#     midpoint cross-check exact in double precision against the normative
#     zeros_100.csv (same protocol as cert_paper9_ratio.py v1.2); the former
#     dps-30 regeneration file is no longer an input.  Hardy-Z sign changes
#     are kept as a redundant SECOND check on the certified balls.
#   - three named constants replace the single target: the certified margins
#     SMOOTH_MARGIN = 5343.90 (downward-rounded lower bound, record E-6) and
#     DISCRETE_MARGIN = 5416.28, plus the historical pre-registered target
#     PREREGISTERED = 5000 (printed only).  Directed asserts run against the
#     two margins; L* is computed from DISCRETE_MARGIN (prints 35.51).
#   - path cascade for zeros_100.csv: script dir, then CWD (no fixed mount).
#
# Requires: python-flint (FLINT/Arb), mpmath, numpy, scipy, sympy.

import sys, time, math
import numpy as np
import flint
from flint import arb, acb
from sympy import primerange
from fractions import Fraction

T0 = time.time()
flint.ctx.prec = 192          # ~57 decimal digits working precision
DELTA_ZERO = arb(10)**(-20)   # zero enclosure radius (certified below)
DELTA_A    = arb(10)**(-22)   # a_j enclosure radius (certified below)
EPS_LO, EPS_HI = Fraction(4,100), Fraction(7,100)   # exact rational window
N_COVER   = 800               # subintervals covering [0.04, 0.07]
SMOOTH_MARGIN   = arb("5343.90")   # certified lower bound, DOWNWARD-rounded (E-6)
DISCRETE_MARGIN = arb("5416.28")   # certified lower bound, discrete comparison
PREREGISTERED   = 5000             # historical pre-registered target (printed only)
PASS = True
def fail(msg):
    global PASS
    print("  CERT-FAIL ", msg); PASS = False

# ---------------------------------------------------------------- 0. zeros ---
import os
def _find_csv(name='zeros_100.csv'):
    sd = os.path.dirname(os.path.abspath(__file__))
    for p in (os.path.join(sd, '..', '..', 'data', name),
              os.path.join(sd, 'data', name),
              os.path.join('data', name),
              os.path.join(sd, name), name):
        if os.path.exists(p): return p
    print(f"  CERT-FAIL  {name} not found (repo data/, script dir, CWD)")
    sys.exit(2)
ref = np.loadtxt(_find_csv(), delimiter=',', skiprows=1)[:, 1]
if len(ref) != 100:
    fail(f"normative layer must contain exactly 100 rows; got {len(ref)}")
    sys.exit(2)
try:
    zs = acb.zeta_zeros(1, 100)
except AttributeError:
    fail("this python-flint build does not expose acb.zeta_zeros"); sys.exit(2)
if len(zs) != 100: fail(f"zeta-zero isolation returned {len(zs)} of 100"); sys.exit(2)
gam  = [z.imag for z in zs]                 # certified balls at working precision
gamf = np.array([float(g.mid()) for g in gam])
bad = int(np.sum(gamf != ref))
if bad: fail(f"midpoint cross-check failed for {bad} of 100"); sys.exit(2)
print("[0] zeros: acb.zeta_zeros(1,100) isolation; midpoint cross-check exact "
      "in double precision against normative zeros_100.csv (all 100)")

# certified enclosures: Hardy Z sign change on [g-d, g+d]
def hardy_Z(t):
    s = acb(arb(1)/4, t/2)
    theta = s.lgamma().imag - t/2*arb.pi().log()
    return (acb(0, theta).exp() * acb(arb(1)/2, t).zeta()).real

t0 = time.time()
for k, g in enumerate(gam, 1):
    prod = hardy_Z(g - DELTA_ZERO) * hardy_Z(g + DELTA_ZERO)
    if not (prod.upper() < 0):
        fail(f"zero enclosure k={k}: no certified sign change")
gam = [g.union(g - DELTA_ZERO).union(g + DELTA_ZERO) for g in gam]  # widen to certified balls
print(f"[0] redundant second check: Hardy-Z sign change around each certified "
      f"ball (radius 1e-20, on-line), {time.time()-t0:.1f}s")

# ---------------------------------------------------- 1. worlds P, A, A' ----
primes = list(primerange(2, 54)); assert len(primes) == 16
li  = lambda t: t.li()                       # arb li (= Ei(log t), t>1)
li2, li53 = li(arb(2)), li(arb(53))
M = li53 - li2
print(f"[1] M = li(53)-li(2) = {M}  (mass wedge M/16 = {M/16})")

# a_j: high-precision midpoints from mpmath, enclosure certified in Arb.
import mpmath as mp
mp.mp.dps = 32
mli2 = mp.li(2); mM = mp.li(53) - mli2
def quantiles(off):
    out = []
    for j in range(1, 17):
        tgt = mli2 + (mp.mpf(j) - off)/16 * mM
        out.append(mp.nstr(mp.findroot(lambda t: mp.li(t) - tgt, 2 + 51*j/16.0), 28))
    return out
A_str  = quantiles(mp.mpf('0.5'))
Ap_str = quantiles(mp.mpf('0.25'))
Am_str = quantiles(mp.mpf('0.75'))

A_balls = []
for j, s in enumerate(A_str, 1):
    aj = arb(s)
    tgt = li2 + (arb(2*j - 1)/32)*M
    lo, hi = li(aj - DELTA_A) - tgt, li(aj + DELTA_A) - tgt
    if not (lo.upper() < 0 and hi.lower() > 0):     # li strictly increasing on t>1
        fail(f"a_{j} enclosure not certified")
    A_balls.append(aj.union(aj - DELTA_A).union(aj + DELTA_A))
Af = np.array([float(s) for s in A_str])
print("[1] all 16 a_j enclosures certified (li monotone, radius 1e-22)")
print("    a_j =", np.round(Af, 4))

# ------------------------------------------- 2. anchor self-check (float64) --
Pf = np.array(primes, float)
def B_disc_f(X, eps):
    w = np.exp(-eps**2*gamf**2)*gamf**2
    lx = np.log(X)
    return float(np.sum(lx**2 * (w[None, :]*np.cos(np.outer(lx, gamf))).sum(axis=1)))
def B_smooth_f(eps):
    w = np.exp(-eps**2*gamf**2)*gamf**2
    z = 1 + 1j*gamf
    F = lambda u: np.exp(z*u)*(u/z - 1/z**2)
    return float(np.sum(w*((F(np.log(53.0)) - F(np.log(2.0))).real)))

TARGETS = {  # series reference anchors (float64-verified)
 0.03: (-145405.08, -12867.45, -1583.28), 0.04: (-46187.20, -2900.28, -794.60),
 0.05: (-19342.55,   -458.63,  -241.15), 0.06: ( -9826.75,    +11.51,  -80.13),
 0.07: ( -5424.69,     -8.17,   -80.78)}
print("[2] anchor self-check (float64 vs series reference anchors):")
ok = True
for eps, (tP, tA, tG) in TARGETS.items():
    vP, vA, vG = B_disc_f(Pf, eps), B_disc_f(Af, eps), B_smooth_f(eps)
    o = (abs(vP-tP) < 0.01 and abs(vA-tA) < 0.01 and abs(vG-tG) < 0.01)
    ok &= o
    print(f"    eps={eps:.2f}  B[P]={vP:12.2f}  B[A]={vA:10.2f}  B_glatt={vG:9.2f}  "
          f"{'PASS' if o else 'FAIL'}")
vApl, vAmi = B_disc_f(np.array([float(s) for s in Ap_str]), 0.05), \
             B_disc_f(np.array([float(s) for s in Am_str]), 0.05)
ok &= abs(vApl - (-2976.68)) < 0.01 and abs(vAmi - 1528.96) < 0.01
print(f"    A'+ = {vApl:.2f}, A'- = {vAmi:.2f}  (targets -2976.68 / +1528.96)")
if not ok: fail("ANCHOR SELF-CHECK FAILED"); sys.exit(1)
print("    ANCHOR: PASS (incl. B[P](0.05) = series anchor -19342.55)")

# --------------------------- 3. eps-free coefficients C_k, I_k, A_k (Arb) ----
logp  = [arb(p).log() for p in primes]
logp2 = [l*l for l in logp]
u2, u53 = arb(2).log(), arb(53).log()

C = []; Aq = []; I = []
for g in gam:
    C.append(sum((l2*(g*l).cos() for l2, l in zip(logp2, logp)), arb(0)))
    la  = [a.log() for a in A_balls]
    Aq.append(sum(((l*l)*(g*l).cos() for l in la), arb(0)))
    z  = acb(1, g)
    F  = lambda u: (z*u).exp()*(u/z - 1/(z*z))
    I.append((F(u53) - F(u2)).real)
g2 = [g*g for g in gam]
dS = [gg*(c - i) for gg, c, i in zip(g2, C, I)]     # coeffs of Delta_smooth
dA = [gg*(c - a) for gg, c, a in zip(g2, C, Aq)]    # coeffs of Delta_atom
print(f"[3] eps-free coefficients C_k, I_k (closed form), A_k computed as Arb balls")

# ------------------- 4. T-DISC: uniform interval certificate on [0.04,0.07] --
def uniform_sup(coeffs, ncov):
    """certified sup over eps in [EPS_LO, EPS_HI] of sum coeffs_k * exp(-eps^2 g_k^2),
    by covering with ncov subintervals (eps as an Arb ball per subinterval)."""
    sup = None
    W = EPS_HI - EPS_LO
    for i in range(ncov):
        a = EPS_LO + W*i/ncov; b = EPS_LO + W*(i+1)/ncov
        e = (arb(a.numerator)/a.denominator).union(arb(b.numerator)/b.denominator)
        e2 = e*e
        val = sum((c*(-e2*gg).exp() for c, gg in zip(coeffs, g2)), arb(0))
        u = val.upper()
        sup = u if sup is None else max(sup, u)
    return sup

t0 = time.time()
supS = uniform_sup(dS, N_COVER)
supA = uniform_sup(dA, N_COVER)
print(f"[4] T-DISC certificate ({N_COVER} covering intervals, {time.time()-t0:.1f}s):")
print(f"    sup_eps (B[P]-B_glatt) <= {supS.str(12)}   "
      f"=> |B[P]-B_glatt| >= {SMOOTH_MARGIN.str(6, radius=False)} (downward-rounded)")
print(f"    sup_eps (B[P]-B[A])    <= {supA.str(12)}   "
      f"=> |B[P]-B[A]|    >= {DISCRETE_MARGIN.str(6, radius=False)} (downward-rounded)")
if not (supS < -SMOOTH_MARGIN):
    fail("smooth margin below certified DOWNWARD bound 5343.90")
if not (supA < -DISCRETE_MARGIN):
    fail("discrete margin below certified bound 5416.28")
print(f"    CERTIFIED: |B[P]-B_glatt| >= {SMOOTH_MARGIN} (downward)  and  "
      f"|B[P]-B[A]| >= {DISCRETE_MARGIN} uniformly on [0.04, 0.07]")
print(f"    historical pre-registered target {PREREGISTERED}: exceeded by both")

# certified grid values for the report (verdict window)
print("    certified point values (Arb balls):")
for eps in [Fraction(4,100), Fraction(5,100), Fraction(6,100), Fraction(7,100)]:
    e = arb(eps.numerator)/eps.denominator; e2 = e*e
    w  = [(-e2*gg).exp() for gg in g2]
    bp = sum((wk*gg*c for wk, gg, c in zip(w, g2, C)), arb(0))
    bg = sum((wk*gg*i for wk, gg, i in zip(w, g2, I)), arb(0))
    ba = sum((wk*gg*a for wk, gg, a in zip(w, g2, Aq)), arb(0))
    print(f"      eps={float(eps):.2f}: B[P]={bp}  B_glatt={bg}  B[A]={ba}")

# explicit d/d(eps) bounds (requirement: state the bound in closed form)
e2lo = (arb(EPS_LO.numerator)/EPS_LO.denominator)**2
S4 = sum(((-e2lo*gg).exp()*gg*gg for gg in g2), arb(0))     # sum e^{-eps_lo^2 g^2} g^4
SP2 = sum(logp2, arb(0))
DP = 2*(arb(EPS_HI.numerator)/EPS_HI.denominator)*SP2*S4
absI = sum(((-e2lo*gg).exp()*gg*gg*(i.abs_upper()) for gg, i in zip(g2, I)), arb(0))
DG = 2*(arb(EPS_HI.numerator)/EPS_HI.denominator)*absI
print(f"    explicit derivative bounds on the window: |d B[P]/d eps| <= {float(DP.upper()):.3e}, "
      f"|d B_glatt/d eps| <= {float(DG.upper()):.3e}")
print(f"    (certificate itself uses the sharper per-subinterval Arb enclosure; "
      f"a grid+Lipschitz variant would need step h <~ {2*300/float(DP.upper()):.1e})")

# --------------------------------- 5. T-INV: explicit constants (Arb upper) --
L = u53
print("[5] T-INV explicit constants:")
# (Part 1) sup |E_pi|, E_pi(t) = pi(t) - (li(t)-li(2)) on [2,53]
# Candidates are kept as BALLS.  pi jumps +1 at each prime and decreases
# continuously in between, so the extrema of |E_pi| on [2,53] are exactly the
# 16 right values |E_pi(p^+)| and the 15 left limits |E_pi(p^-)|, p=3,...,53.
ARGMAX_KEY = "53^-"                      # claimed rank 1; certified below
RUNNER_KEY = "41^-"                      # claimed rank 2; certified below
cands = []
for i, p in enumerate(primes):           # on [p_i, p_{i+1}) : pi = i+1
    cands.append((f"{p}^+", abs(arb(i + 1) - (li(arb(p)) - li2))))
    if i < 15:
        cands.append((f"{primes[i+1]}^-",
                      abs(arb(i + 1) - (li(arb(primes[i+1])) - li2))))
# The FULL ranking is certified, not just the top: order the candidates by
# midpoint (a heuristic whose provenance is irrelevant) and then verify the
# chain of directed separations lower(v_j) > upper(v_{j+1}) in ball
# arithmetic.  A chain of 30 verified links is a total order on the 31
# candidates, so every rank statement the manuscript makes -- maximiser,
# runner-up, or any other -- is covered by one certificate.
ranked = sorted(cands, key=lambda kv: -float(kv[1].mid()))
for (ka, va), (kb, vb) in zip(ranked, ranked[1:]):
    if not bool(va.lower() > vb.upper()):
        fail(f"finite ordering of |E_pi| is not certified: "
             f"directed separation {ka} > {kb} failed")
if ranked[0][0] != ARGMAX_KEY or ranked[1][0] != RUNNER_KEY:
    fail(f"certified ranking disagrees with the claimed ranks: got "
         f"{ranked[0][0]}, {ranked[1][0]}; expected {ARGMAX_KEY}, {RUNNER_KEY}")
star, runner_v = ranked[0][1], ranked[1][1]
third_k, third_v = ranked[2]
supEpi = star
print(f"    sup_[2,53] |pi(t) - (li(t)-li(2))| = {supEpi.str(12)}   "
      f"(certified total ordering of {len(cands)} candidates, "
      f"{len(ranked)-1} directed links; rank 1 {ARGMAX_KEY}, "
      f"rank 2 {RUNNER_KEY} {runner_v.str(8)}, rank 3 {third_k}; "
      f"margins 1-2 {(star - runner_v).str(8)}, "
      f"2-3 {(runner_v - third_v).str(8)})")
print(f"    |16 - M| = {float(abs(16-M)):.6f}")
# (Part 2) midpoint-in-li constants for g_B: sup|g'|, sup|g''| -> sup|h''|
print("    quadrature constants (per eps):  E-bound = M^2/384 * sup|h''|,")
print("    sup|h''| <= L^2 sup|g''| + (L/2) sup|g'|,  L = log 53")
for eps in [Fraction(4,100), Fraction(5,100), Fraction(6,100), Fraction(7,100)]:
    e = arb(eps.numerator)/eps.denominator; e2 = e*e
    G1 = arb(0); G2 = arb(0)
    for g, gg in zip(gam, g2):
        w = (-e2*gg).exp()*gg
        G1 += w*(2*L + g*L*L)/2
        G2 += w*(2 + 2*L + 4*g*L + g*L*L + gg*L*L)/4
    H2 = L*L*G2 + L*G1/2
    Ebnd = M*M/384*H2
    epsf = float(eps)
    Eobs = B_disc_f(Af, epsf) - float((16/M).mid())*B_smooth_f(epsf)
    print(f"      eps={epsf:.2f}: sup|g'|<={float(G1.upper()):.3e}  sup|g''|<={float(G2.upper()):.3e}"
          f"  |E|<= {float(Ebnd.upper()):.3e}   (observed E = {Eobs:+.1f})")

# ------------------------------- 6. d_BL (flat metric) via LP + dual bound ---
print("[6] d_BL(mu_P, mu_A), flat/bounded-Lipschitz convention "
      "(||f||_inf <= 1 and Lip(f) <= 1):")
X  = np.concatenate([Pf, Af]); n = 32
cw = np.concatenate([np.log(Pf)**2, -np.log(Af)**2])     # objective: int f d(mu_P - mu_A)
print(f"    masses: mu_P(1) = {np.sum(np.log(Pf)**2):.4f}, "
      f"mu_A(1) = {np.sum(np.log(Af)**2):.4f}, "
      f"difference = {np.sum(np.log(Pf)**2)-np.sum(np.log(Af)**2):+.4f}")
rowsA = []; rhs = []
for i in range(n):
    for j in range(i+1, n):
        d = abs(X[i] - X[j])
        r = np.zeros(n); r[i], r[j] = 1, -1
        rowsA += [r, -r]; rhs += [d, d]
Aub = np.vstack(rowsA); bub = np.array(rhs)
from scipy.optimize import linprog
res = linprog(-cw, A_ub=Aub, b_ub=bub, bounds=[(-1, 1)]*n, method='highs')
assert res.status == 0
dbl = -res.fun
# rigorous LOWER bound: shrink optimal f slightly, verify feasibility in Arb
fstar = res.x * (1 - 1e-9)
lp2b  = logp2 + [(a.log())**2 for a in A_balls]
Xb    = [arb(p) for p in primes] + A_balls
LBv = arb(0); feas = True
for i in range(n):
    if not (abs(arb(f"{fstar[i]:.17g}")).upper() <= 1): feas = False
    for j in range(i+1, n):
        lhs = abs(arb(f"{fstar[i]:.17g}") - arb(f"{fstar[j]:.17g}"))
        if not (lhs.upper() <= abs(Xb[i]-Xb[j]).lower()): feas = False
if not feas:
    fail("d_BL primal witness is not certified feasible in Arb")
sgn = [1]*16 + [-1]*16
LBv = sum((s*w*arb(f"{f:.17g}") for s, w, f in zip(sgn, lp2b, fstar)), arb(0))
# rigorous UPPER bound: weak duality with residual absorption (||f||_inf <= 1)
y = np.clip(res.ineqlin.marginals * -1, 0, None)     # duals of A_ub f <= b_ub
zl = np.clip(res.lower.marginals, 0, None)           # duals of -f <= 1
zu = np.clip(res.upper.marginals * -1, 0, None)      # duals of  f <= 1
r_res = Aub.T @ y + zu - zl - cw                     # A^T y + z_u - z_l = c (float resid)
# absorb: c.f = (A^T y + zu - zl - r).f <= b.y + 1.(zu+zl) + ||r||_1
UB = arb(0)
k = 0
for i in range(n):
    for j in range(i+1, n):
        dij = abs(Xb[i]-Xb[j])
        UB += arb(f"{y[k]:.17g}")*dij + arb(f"{y[k+1]:.17g}")*dij
        k += 2
UB += sum((arb(f"{v:.17g}") for v in np.concatenate([zl, zu])), arb(0))
# residual recomputed in Arb (mass weights as balls)
cb = [s*w for s, w in zip(sgn, lp2b)]
r1 = arb(0)
for i in range(n):
    col = arb(0); k = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i: col += arb(f"{y[k]:.17g}") - arb(f"{y[k+1]:.17g}")
            if b == i: col += -arb(f"{y[k]:.17g}") + arb(f"{y[k+1]:.17g}")
            k += 2
    col += arb(f"{zu[i]:.17g}") - arb(f"{zl[i]:.17g}") - cb[i]
    r1 += abs(col)
UBtot = UB + r1
print(f"    LP value (float): {dbl:.6f}")
print(f"    rigorous: {float(LBv.lower()):.6f} <= d_BL <= {float(UBtot.upper()):.6f}"
      f"   (primal feasible in Arb: {feas}; dual residual ||r||_1 = {float(r1.upper()):.2e})")
if not (float(LBv.lower()) <= float(UBtot.upper())):
    fail("d_BL enclosure is inconsistent (lower bound exceeds upper bound)")
dbl_up = float(UBtot.upper())
Lstar = float(DISCRETE_MARGIN) / dbl_up
print(f"    T-FACT threshold: L* = Delta_disc / d_BL >= {Lstar:.2f}   "
      f"(from the certified discrete margin {DISCRETE_MARGIN.str(6, radius=False)})")
# scale of B's own BL-norm as a linear functional (for honest contrast), eps=0.05
Kinf = float(sum(((-arb('0.0025')*gg).exp()*gg for gg in g2), arb(0)).upper())
print(f"    contrast: ||K_eps||_inf <= sum w g^2 = {Kinf:.1f} at eps=0.05 "
      f"(B itself is BL with large constant — no contradiction)")

# --------------------------------------- 7. scatter / mass-wedge balance -----
print("[7] scatter band & mass wedge at eps=0.05 (float64 diagnostics):")
print(f"    B[A] = {B_disc_f(Af,0.05):.2f}, B[A'+] = {vApl:.2f}, B[A'-] = {vAmi:.2f}"
      f"  -> band width {vAmi - vApl:.1f} = O(10^3)")
print(f"    (16/M)*B_glatt = {float((16/M).mid())*B_smooth_f(0.05):.2f}  "
      f"(wedge factor M/16 = {float((M/16).mid()):.4f} explicit, not noise)")
print(f"    eps=0.06 'flip': B[A] = {B_disc_f(Af,0.06):+.2f} but B_glatt = "
      f"{B_smooth_f(0.06):.2f} < 0 — quadrature artifact, verdict carrier does not flip")

print(f"\n{'='*70}\nSMOOTH/DISCRETE MARGIN CERTIFICATE v2.0: {'PASS' if PASS else 'FAIL'}   "
      f"(total {time.time()-T0:.1f}s, prec={flint.ctx.prec} bits)")
sys.exit(0 if PASS else 1)
