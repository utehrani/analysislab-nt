#!/usr/bin/env python3
"""
cert_paper9_ratio.py v1.2 — ball certificate for the two-point separation
of the limiting spectral ratio (Paper 9, Theorem 3.1).

Changes in v1.2 (cycle 5)
  Series-harmonised ordinate location: the normative layer is found by the
  Papers-1-8 candidate cascade with header sniffing (an explicitly given
  --zeros path is authoritative and fail-closed; only the default cascades
  through script directory and CWD; a headerless file is detected by a
  numeric first cell).  Print and docstring notation now writes the
  certified two-point separation as R(u1)-R(u2) — the limiting ratio is
  what the certificate bounds; the finite-cutoff quantity r(x) of the
  manuscript converges to it.  The computational core is unchanged.

Changes in v1.1 (cycle 4)
  The ordinates are now obtained as CERTIFIED ENCLOSURES from Arb's
  zeta-zero isolation (acb.zeta_zeros) at working precision, instead of as
  an mpmath regeneration wrapped in an assumed input radius.  The former
  GAMMA_RADIUS hypothesis is gone: the input enclosures are themselves
  proved, so the certificate chain is closed end to end.  Retained: the
  mandatory cross-check of the ball midpoints against the normative
  double-precision layer, the Rayleigh lower bound, the even-power trace
  upper bound, outward interval division, and STOP semantics.

New in cycle 2. Before this script the separation was quoted from a
floating-point evaluation and could not carry the label "certified"; it now
does, because the bound below is produced in interval arithmetic with the
ordinates entering as balls.

Certifies, in Arb ball arithmetic via python-flint:

      R(u1) - R(u2)  >=  Delta_0 ,      R(u) := lambda_max(Q(u)) / Tr Q(u),

with u1 = 17.6170, u2 = 7.1815 (exact decimals).

Structure of the certificate
----------------------------
  lower bound on lambda_max(Q(u1))   Rayleigh quotient <v,Qv>/<v,v> of an
                                     explicit vector v.  Rigorous for ANY v,
                                     so the provenance of v is irrelevant to
                                     the validity of the bound; v is obtained
                                     deterministically in double precision.
  upper bound on lambda_max(Q(u2))   (Tr Q^m)^{1/m} for even m.  Since
                                     Tr Q^m = sum lambda_i^m >= lambda_max^m
                                     for even m, this majorises lambda_max
                                     with no positivity assumption.
  both traces                        summed directly as balls.
  ratio                              outward-rounded interval division.

Ordinate protocol
-----------------
Layer 1 (normative) is zeros_100.csv in double precision; it defines the
anchors.  The certificate layer consists of Arb's certified enclosures of
the first hundred non-trivial zeros; their ball midpoints must reproduce
layer 1 exactly in double precision (mandatory cross-check), and the
enclosures then enter the arithmetic as they are — the input uncertainty is
part of the certificate, not a hypothesis about it.  No enclosure list is
archived; this script is.

Usage:  python3 cert_paper9_ratio.py [--zeros zeros_100.csv] [--m 64]
"""

import argparse, csv, sys

PREC = 320                    # working precision in bits
EXPECTED_N = 100              # Theorem 3.1 is a statement at (N, eps) = (100, 0.05)
EPS = "0.05"
U1, U2 = "17.6170", "7.1815"
DELTA_0_TARGET = "0.0840112"


def fail(msg):
    print(f"\n  STOP: {msg}")
    sys.exit(1)


def load_layer1(path):
    """Series-standard cascade + header sniffing (fail-closed; an explicit
    path is authoritative, only the default cascades). Explicitness is
    decided by ``path is not None``, never by comparing against the default
    base name: a caller who passes --zeros zeros_100.csv means that file in
    the current directory and must not be silently served the copy next to
    the script (cycle-26 finding N1b)."""
    import os
    if path is not None:
        candidates = [path]
    else:
        sd = os.path.dirname(os.path.abspath(__file__))
        candidates = [os.path.join(sd, "..", "..", "data", "zeros_100.csv"),
                      os.path.join(sd, "data", "zeros_100.csv"),
                      os.path.join("data", "zeros_100.csv"),
                      os.path.join(sd, "zeros_100.csv"), "zeros_100.csv"]
    for c in candidates:
        if os.path.exists(c):
            with open(c) as fh:
                first = fh.readline().split(",")[0].strip()
            skip = 0 if first.lstrip("-").replace(".", "", 1).isdigit() else 1
            vals = []
            with open(c, newline="") as f:
                for i, line in enumerate(f):
                    if i < skip or not line.strip():
                        continue
                    vals.append(float(line.split(",")[-1]))
            print(f"  ordinates: {len(vals)} from {c} "
                  f"(header {'skipped' if skip else 'absent'})")
            return vals
    fail(f"ordinate list not found; tried {candidates}")


def certified_ordinates(n, layer1):
    """Obtain the first n ordinates as certified Arb enclosures and run the
    mandatory cross-check of the ball midpoints against the normative
    layer."""
    from flint import acb
    try:
        zs = acb.zeta_zeros(1, n)
    except AttributeError:
        fail("this python-flint build does not expose acb.zeta_zeros; "
             "the certificate cannot be built in this environment")
    if len(zs) != n:
        fail(f"zeta-zero isolation returned {len(zs)} of {n} enclosures")
    gam = [z.imag for z in zs]
    bad = sum(1 for k in range(n) if float(gam[k].mid()) != layer1[k])
    if bad:
        fail(f"midpoint cross-check failed for {bad} of {n} enclosures")
    print(f"  ordinate protocol: {n} certified enclosures from Arb zeta-zero "
          f"isolation; midpoint cross-check against the normative layer "
          f"exact in double precision for all {n}")
    return gam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zeros", default=None)
    ap.add_argument("--m", type=int, default=64)
    args = ap.parse_args()
    if args.m % 2:
        fail("the trace-power bound needs an even exponent")

    print("=" * 62)
    print("cert_paper9_ratio.py v1.2 — two-point separation, ball certificate")
    print("=" * 62)

    try:
        from flint import arb, arb_mat, ctx
    except ImportError:
        fail("python-flint is not available; the certificate cannot be built "
             "in this environment")
    ctx.prec = PREC
    print(f"  Arb working precision: {PREC} bits")

    layer1 = load_layer1(args.zeros)
    if len(layer1) != EXPECTED_N:
        fail(f"the normative ordinate layer must contain exactly {EXPECTED_N} "
             f"rows for the Theorem-3.1 certificate; got {len(layer1)}")
    n = EXPECTED_N
    gam = certified_ordinates(n, layer1)
    eps = arb(EPS)
    w = [(-eps ** 2 * g ** 2 / 2).exp() for g in gam]

    def h(theta, u):
        return ((theta * u).cos() + theta * (theta * u).sin()) / (1 + theta ** 2)

    def build_Q(u_str):
        u = arb(u_str)
        Q = arb_mat(n, n)
        for k in range(n):
            Q[k, k] = w[k] ** 2 / 2 * (1 - h(gam[k], u))
            for l in range(k + 1, n):
                v = w[k] * w[l] / 2 * (h((gam[k] - gam[l]) / 2, u)
                                       - h((gam[k] + gam[l]) / 2, u))
                Q[k, l] = v
                Q[l, k] = v
        return Q

    print("  building Q(u1), Q(u2) as ball matrices ...")
    Q1, Q2 = build_Q(U1), build_Q(U2)

    def trace(Q):
        t = arb(0)
        for i in range(n):
            t += Q[i, i]
        return t

    tr1, tr2 = trace(Q1), trace(Q2)
    print(f"  Tr Q(u1) in {tr1.str(12)}")
    print(f"  Tr Q(u2) in {tr2.str(12)}")
    if not (tr1 > 0) or not (tr2 > 0):
        fail("trace not certified positive")

    # ---- lower bound on lambda_max(Q(u1)) via a Rayleigh quotient ----------
    import numpy as np
    g1 = np.array([float(g.mid()) for g in gam])
    e = float(EPS)
    wf = np.exp(-e ** 2 * g1 ** 2 / 2)

    def hf(th, u):
        return (np.cos(th * u) + th * np.sin(th * u)) / (1 + th * th)

    def Qf(u):
        M = np.outer(wf, wf) / 2 * (hf((g1[:, None] - g1[None, :]) / 2, u)
                                    - hf((g1[:, None] + g1[None, :]) / 2, u))
        np.fill_diagonal(M, wf ** 2 / 2 * (1 - hf(g1, u)))
        return M

    vec = np.linalg.eigh(Qf(float(U1)))[1][:, -1]
    v = arb_mat(n, 1)
    for i in range(n):
        v[i, 0] = arb(float(vec[i]))          # exact doubles, no rounding claim
    Qv = Q1 * v
    num = arb(0); den = arb(0)
    for i in range(n):
        num += v[i, 0] * Qv[i, 0]
        den += v[i, 0] ** 2
    rayleigh = num / den
    lam1_low = arb(rayleigh.lower())
    print(f"  lambda_max(Q(u1)) >= {lam1_low.str(12)}   (Rayleigh quotient)")

    # ---- upper bound on lambda_max(Q(u2)) via (Tr Q^m)^{1/m} --------------
    m = args.m
    P, e_, base = None, m, Q2
    bits = bin(m)[2:]
    P = base
    for b in bits[1:]:
        P = P * P
        if b == "1":
            P = P * base
    trPm = arb(0)
    for i in range(n):
        trPm += P[i, i]
    if not (trPm > 0):
        fail("trace of the matrix power not certified positive")
    lam2_up = arb((trPm ** (arb(1) / m)).upper())
    print(f"  lambda_max(Q(u2)) <= {lam2_up.str(12)}   "
          f"(trace-power bound, m = {m})")

    # ---- ratios and separation -------------------------------------------
    r1_low = arb((lam1_low / arb(tr1.upper())).lower())
    r2_up = arb((lam2_up / arb(tr2.lower())).upper())
    sep = arb((r1_low - r2_up).lower())

    print(f"\n  R(u1) >= {r1_low.str(15)}")
    print(f"  R(u2) <= {r2_up.str(15)}")
    print(f"  R(u1) - R(u2) >= {sep.str(15)}")

    target = arb(DELTA_0_TARGET)
    ok = bool(sep >= target)
    print(f"\n  registered separation Delta_0 = {DELTA_0_TARGET}")
    print(f"  certificate {'CONFIRMS' if ok else 'DOES NOT CONFIRM'} "
          f"R(u1) - R(u2) >= Delta_0")
    print("=" * 62)
    if not ok:
        print("  RESULT: FAIL")
        return 1
    print("  RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
