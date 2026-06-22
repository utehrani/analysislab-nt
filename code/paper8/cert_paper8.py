#!/usr/bin/env python3
# cert_paper8.py
# Rigorous interval certificate for Paper 8, Theorem 6.2 (negativity at kappa=53):
#     B_GW^infty(53, eps) < 0   for all  0 < eps <= 0.05.
#
# This is the CERTIFICATE (not the consistency gate verify_paper8.py).
# It is the artifact the word "certified" in the paper refers to.
#
# Method (archimedean L^1 bound, fourth-order integration by parts):
#   Multiply the second-moment identity by 8 sqrt(pi) eps^3 (log p)^2 and sum:
#   the prime side is  S_3(53) + Cross(eps),  with
#       Cross(eps) = sum_{p<=53} (log p)^2 [ sum_{n=q^m, +-, (n,+)!=(p,+)}
#                        (Lambda(n)/sqrt n) h(x_{n,+-})  +  log(pi) h(x_{l_p}) ],
#       h(x) = (1 - 2x) e^{-x},   x_{n,+-} = (log n -+ l_p)^2 / (4 eps^2).
#   The (p,+) diagonal (x=0, value (log p)^3/sqrt p) is exactly S_3 and is pulled out.
#   h decreases on [0,3/2], increases on [3/2,inf), h(3/2) = -2 e^{-3/2}, h(inf)=0.
#   Conductor: A_p = h(x)/(4 sqrt(pi) eps^3); in Cross = -8 sqrt(pi) eps^3 B_GW it
#   appears as +(log pi) h(x_{l_p}) (positive coeff; see _derive_conductor, sympy).
#   On an interval eps in [a,b], each x ranges in [x(b), x(a)]; both the prime-power
#   terms and the positive log-pi conductor term are bounded below by min h, giving a
#   closed-form interval
#   lower bound LB[a,b] <= min_{[a,b]} Cross.  The interval criterion
#       S_3(53) + LB[a,b] > 8 sqrt(pi) b^3 * 2693   ==>   B_GW^infty(53,eps) < 0 on [a,b].
#
# The 7 certificate requirements:
#   (1) Arb ball arithmetic (FLINT 3, python-flint) with directed outward rounding.
#   (2) Every interval endpoint is an exact decimal rational.
#   (3) Every prime-power interval contribution is recomputed here.
#   (4) The discarded prime-power tail (|log n -+ l_p| >= 3) is bounded rigorously.
#   (5) For each covering interval: lower bound, right-hand side, and margin are printed.
#   (6) On overlap / gap / non-positive margin the script exits with code 1.
#   (7) Backend version, precision, source hash and parameter hash are logged.
#
# Usage:  python cert_paper8.py
# Requires: python-flint (FLINT 3 / Arb), sympy.

import sys
import hashlib
import math

try:
    import flint
    from flint import arb, acb
except ImportError:
    print("ERROR: python-flint (FLINT 3 / Arb) is required for the certificate.")
    sys.exit(2)

from sympy import primerange, factorint

PREC_BITS = 256            # ~77 decimal digits of working precision
flint.ctx.prec = PREC_BITS
KAPPA = 53
# GW_ARCH_BOUND is CERTIFIED below by certify_I4(): the weighted archimedean
# total satisfies  G_w = ceil( I_4(1/4)/(2 pi) * sum_{p<=53}(log p)^{-2} ) <= 2693,
# with I_4(1/4) <= 3561.1 produced as an Arb OUTPUT (not a numerical input) by
# recomputing the five-term Leibniz bounds of the archimedean residual lemma in Arb.
GW_ARCH_BOUND = 2693
TAIL_CUTOFF = 3            # prime powers with |log n -+ l_p| < 3 are summed exactly

# Covering intervals of (0, 0.05], endpoints as exact decimal rationals.
# a is taken as a tiny positive rational standing for the open left end (0, .].
COVER = [
    ("1/1000000000", "47/1000"),   # (0, 0.047]
    ("47/1000",      "1/20"),      # [0.047, 0.05]
]

PASS = True


def fail(msg):
    global PASS
    print(f"  CERT-FAIL  {msg}")
    PASS = False


def is_prime_power(n):
    f = factorint(n)
    if len(f) == 1:
        q = next(iter(f)); return (q, f[q])
    return None


primes = list(primerange(2, KAPPA + 1))
# prime powers up to kappa * e^3 (the |log n - l_p| < 3 window for the largest p)
PP = [(n,) + is_prime_power(n) for n in range(2, int(KAPPA * math.exp(TAIL_CUTOFF)) + 2)
      if is_prime_power(n)]
logpi = arb.pi().log()

# S_3(53) = sum (log p)^3 / sqrt p, as an Arb ball (lower bound used below).
S3 = arb(0)
for p in primes:
    S3 += arb(p).log() ** 3 / arb(p).sqrt()


def _derive_conductor():
    """Anchor the conductor term symbolically; this fixes the sign convention used by the interval certificate.
    A_p := (1/2pi) int_{-inf}^{inf} t^2 e^{-eps^2 t^2} cos(tL) dt  (paper def, line 704).
    sympy gives  A_p = (2 eps^2 - L^2) e^{-L^2/4eps^2} / (8 sqrt(pi) eps^5)
               = 1/(4 sqrt(pi) eps^3) * (1 - 2x) e^{-x},   x = L^2/4eps^2,
    and the certificate's  h(x) = (1-2x) e^{-x}, so  A_p = h(x)/(4 sqrt(pi) eps^3)  (NOT -h).
    Conductor in N_GW = 8 sqrt(pi) eps^3 * (-1/2 (log pi) A_p) = -(log pi) h(x);
    Cross carries -N_GW, so the conductor in Cross is +(log pi) h(x) (POSITIVE),
    hence the lower bound cross_LB uses +logpi * hmin (parallel to the +hmin prime
    cross terms).  Asserts the identity (1 - L^2/2eps^2) e^{-L^2/4eps^2} = h(x)."""
    import sympy as sp
    t, e, L = sp.symbols('t epsilon L', positive=True)
    I = sp.integrate(t**2 * sp.exp(-e**2 * t**2) * sp.cos(t * L), (t, -sp.oo, sp.oo))
    Ap = sp.simplify(I / (2 * sp.pi))
    x = L**2 / (4 * e**2)
    target = (1 - 2 * x) * sp.exp(-x) / (4 * sp.sqrt(sp.pi) * e**3)   # h(x)/(4 sqrt pi eps^3)
    assert sp.simplify(Ap - target) == 0, "conductor A_p != h(x)/(4 sqrt pi eps^3)"
    return True


_CONDUCTOR_OK = _derive_conductor()


def h(x):
    """h(x) = (1-2x)e^{-x} as an Arb ball; for x>600, |h|<=(2x-1)e^{-x} is below 1e-250."""
    if float(x.mid()) > 600:
        return arb(0, 1e-250)
    return (1 - 2 * x) * (-x).exp()


def _arb_lt(a, b):
    """True only if a < b is rigorously decided (a.upper() < b.lower())."""
    return float((b - a).lower()) > 0.0


def hmin(xb, xa):
    """Rigorous LOWER enclosure of min_{x in [xb,xa]} h(x).
    h is unimodal with global min h(3/2); use Arb bounds (not float mids) to
    decide the bracket, conservatively including h(3/2) when 1.5 may lie in
    [xb, xa]. Returns the candidate with the smallest lower bound."""
    cands = [h(xa), h(xb)]
    half = arb(3) / 2
    # conservative bracket test: 1.5 possibly in [xb, xa]
    if not (_arb_lt(half, xb) or _arb_lt(xa, half)):
        cands.append(h(half))
    return min(cands, key=lambda c: float(c.lower()))


def hmax(xb, xa):
    """Rigorous UPPER enclosure of max h over [xb,xa]; the max sits at an
    endpoint (h decreases then increases), so compare the two endpoints."""
    ha, hb = h(xa), h(xb)
    return ha if float(ha.upper()) > float(hb.upper()) else hb


def cross_LB(a, b):
    """Closed-form interval lower bound LB[a,b] <= min_{eps in [a,b]} Cross(eps)."""
    a, b = arb(a), arb(b)
    LB = arb(0)
    for p in primes:
        lp = arb(p).log()
        inner = arb(0)
        for (n, q, m) in PP:
            ln = arb(n).log()
            for sign in (1, -1):                    # x_{n,+} uses (ln - lp); x_{n,-} uses (ln + lp)
                if n == p and sign == 1:            # (p,+) is the x=0 diagonal = S_3 (pulled out)
                    continue
                d = ln - sign * lp
                if abs(float(d.mid())) >= TAIL_CUTOFF:
                    continue                        # tail, bounded separately
                xa = d ** 2 / (4 * a ** 2)
                xb = d ** 2 / (4 * b ** 2)
                inner += (arb(q).log() / arb(n).sqrt()) * hmin(xb, xa)   # positive weight -> min h
        xa = lp ** 2 / (4 * a ** 2)
        xb = lp ** 2 / (4 * b ** 2)
        inner += logpi * hmin(xb, xa)               # conductor +(log pi) h(x): A_p=h(x)/(4 sqrt pi eps^3),
        LB += lp ** 2 * inner
    return LB


def _tail_branch(lp, b, x_of_k, k_lo, k_hi):
    """Sum the Gauss/Chebyshev majorant over log-bins k in [k_lo, k_hi), PLUS an
    explicit remainder for k>=k_hi.  On bin [k,k+1): sum_{n} Lambda(n) =
    psi(e^{k+1}) - psi(e^k) <= psi(e^{k+1}) <= 1.04 e^{k+1} (Chebyshev);
    1/sqrt(n) <= e^{-k/2}; |h(x)| <= (2 x_k - 1) e^{-x_k} with x_k = x_of_k(k)
    the smallest x on the bin.  The successive ratio phi(k+1)/phi(k) =
    e^{1/2} e^{-(x_{k+1}-x_k)} < 1/2 once x grows (it does, x quadratic in k), so
    the tail k>=k_hi is bounded by phi(k_hi)/(1-1/2) = 2 phi(k_hi)."""
    s = arb(0)
    def phi(k):
        kk = arb(k); x = x_of_k(kk)
        return arb("1.04") * (kk + 1).exp() * (-kk / 2).exp() * (2 * x - 1) * (-x).exp()
    for k in range(k_lo, k_hi):
        s += phi(k)
    s += 2 * phi(k_hi)            # geometric remainder (ratio < 1/2)
    return s


def tail_bound(b):
    """Rigorous Gauss/Chebyshev majorant for the discarded prime powers.
    cross_LB keeps only |log n -+ l_p| < 3 and discards the rest; there are
    THREE discarded families:
      (+, high) sign=+1, log n - l_p >= 3   -> binned, x_k=(k - l_p)^2/(4 b^2)
      (+, low)  sign=+1, log n - l_p <= -3  -> log n in [0, l_p-3], distance >= 3
      (-)       sign=-1, log n + l_p >= 3   -> binned, x_k=(k + l_p)^2/(4 b^2)
    For the two BINNED families the nearest bin end is k (smallest distance,
    largest |h|), so x_k uses k.  The (+, low) family is NOT binned: across the
    whole region log n in [0, l_p-3] the distance |l_p - log n| >= 3, i.e.
    x >= x_min = 9/(4 b^2); since |h| decreases for x>3/2 and sum Lambda/sqrt(n)
    <= psi(e^{l_p-3}) <= 1.04 e^{l_p-3} (Chebyshev, 1/sqrt(n)<=1), one CLOSED
    term suffices.  (Binning this family with the bin end k+1 would reach into
    the kept |.|<3 zone for primes whose boundary bin straddles l_p-3 and badly
    over-count -- e.g. p=23 has no discarded low prime power at all.)
    On each binned bin sum_{n} Lambda(n) <= psi(e^{k+1}) <= 1.04 e^{k+1}, and
    |h| decays like e^{-x_k} with x_k quadratic in k -- a Gaussian in k that
    dominates the single-exponential bin growth.  The window W is chosen so the
    first omitted term has x > (W)^2/(4 b^2) > 6e5, i.e. < 10^{-400}; the whole
    omitted remainder (a super-exponentially decreasing series) is then far
    below 10^{-380} and is neglected rigorously."""
    import math
    b = arb(b)
    tot = arb(0)
    W = 80
    x_min = arb(9) / (4 * b ** 2)
    h_min = (2 * x_min - 1) * (-x_min).exp()          # |h| bound at distance exactly 3
    for p in primes:
        lp = arb(p).log()
        wp = lp ** 2
        flp = float(lp)
        # (+, high): partial boundary bin [l_p+3, ceil(l_p+3)) (distance>=3 -> x_min), then full bins
        k0 = int(math.ceil(flp + 3.0))
        tot += wp * arb("1.04") * arb(k0).exp() * h_min      # boundary bin: count<=psi(e^{k0}), 1/sqrt(n)<=1
        tot += wp * _tail_branch(lp, b, lambda kk: (kk - lp) ** 2 / (4 * b ** 2), k0, k0 + W)
        # (+, low): closed-form, whole region has distance >= 3
        if flp - 3.0 > 0.0:
            tot += wp * arb("1.04") * (lp - 3).exp() * h_min
        # (-): log n + l_p >= 3, i.e. k >= 3 - l_p (>=0)
        km0 = max(0, int(math.ceil(3.0 - flp)))
        if km0 > 0:                                          # partial boundary bin only when 3-l_p>0
            tot += wp * arb("1.04") * arb(km0).exp() * h_min
        tot += wp * _tail_branch(lp, b, lambda kk: (kk + lp) ** 2 / (4 * b ** 2), km0, km0 + W)
    return tot



_S_TRUNC = 400             # series truncation for the digamma sums S_1, S

def _S1_maj(tau):
    """RIGOROUS upper bound for S_1(tau)=sum_{n>=0} a_n/(a_n^2+tau^2)^2, a_n=n+1/4.
    Option B: the TRUE series truncated at N=_S_TRUNC plus a rigorous constant
    tail, sum_{n>N} a_n/(a_n^2+tau^2)^2 <= sum_{n>N} 1/a_n^3 <= 1/(2(N+3/4)^2).
    Verified maj >= true for tau in {0.5,1,2,3,5,10,20} (the old closed-form
    'integral' majorant UNDER-shot from tau~2 because the summand peaks at
    a=tau/sqrt 3 inside the sum; this truncation avoids any monotonicity claim)."""
    N = _S_TRUNC
    s = arb(0)
    for n in range(N + 1):
        an = arb(n) + arb(1) / 4
        s += an / (an ** 2 + tau ** 2) ** 2
    return s + 1 / (2 * (arb(N) + arb(3) / 4) ** 2)


def _S_maj(tau):
    """RIGOROUS upper bound for S(tau)=sum a_n/(a_n^2+tau^2)^3 (truncation + tail
    sum_{n>N} 1/a_n^5 <= 1/(4(N+3/4)^4)); same Option-B construction as _S1_maj."""
    N = _S_TRUNC
    s = arb(0)
    for n in range(N + 1):
        an = arb(n) + arb(1) / 4
        s += an / (an ** 2 + tau ** 2) ** 3
    return s + 1 / (4 * (arb(N) + arb(3) / 4) ** 4)


def _self_test_majorants():
    """Self-test gate: confirm maj >= true at sample tau before certify_I4 runs.
    Compares the Arb majorant lower edge against a high-term reference sum."""
    for ft in ("0.5", "1", "2", "3", "5", "10", "20"):
        tau = arb(ft)
        ref1 = sum((arb(n) + arb(1) / 4) / ((arb(n) + arb(1) / 4) ** 2 + tau ** 2) ** 2
                   for n in range(4000))
        refS = sum((arb(n) + arb(1) / 4) / ((arb(n) + arb(1) / 4) ** 2 + tau ** 2) ** 3
                   for n in range(4000))
        if not (float(_S1_maj(tau).lower()) >= float(arb(ref1).upper())):
            fail(f"_S1_maj not an upper bound at tau={ft}")
        if not (float(_S_maj(tau).lower()) >= float(arb(refS).upper())):
            fail(f"_S_maj not an upper bound at tau={ft}")



def certify_I4():
    """Certify I_4(1/4) <= 3561.1 and hence G_w <= GW_ARCH_BOUND in Arb, by
    recomputing the five Leibniz terms of the archimedean residual bound as rigorous balls.
    Returns (I4_upper, Gw_upper).  Main terms are exact (pi * psi'(1/4) multiples);
    T_2 is the exact closed form 1920 sqrt(pi) eps; the residual majorants
    A_1 (log weight), A_3, A_4 (digamma series, bounded by closed-form majorants)
    are evaluated by Arb rigorous integration (acb.integral)."""
    eps = arb(1) / 4
    psi1 = arb.pi() ** 2 + 8 * arb.const_catalan()       # psi'(1/4) = pi^2 + 8 G
    sqrtpi = arb.pi().sqrt()
    # exact closed forms
    T2 = 1920 * sqrtpi * eps                              # = 4*int|phi'''|*20/(1+t), Gaussian moments
    main = (9 * arb.pi() + 6 * arb.pi() + arb("3.75") * arb.pi()) * psi1   # T3+T4+T5 main terms
    # residual majorants by rigorous integration (upper edge of each ball)
    def fA1(z, _):
        t = z
        poly = 24 * eps ** 2 + 156 * eps ** 4 * t ** 2 + 112 * eps ** 6 * t ** 4 + 16 * eps ** 8 * t ** 6
        return poly * (-eps ** 2 * t ** 2).exp() * ((2 + t).log() + 8)
    def fA3(z, _):
        t = z
        return (10 * eps ** 2 * t ** 2 + 4 * eps ** 4 * t ** 4) * (-eps ** 2 * t ** 2).exp() * _S1_maj(t / 2)
    def fA4(z, _):
        t = z
        return t ** 4 * (-eps ** 2 * t ** 2).exp() * _S_maj(t / 2)
    A1 = acb.integral(fA1, 0, 200).real
    A3 = 9 * acb.integral(fA3, 0, 200).real
    A4 = 12 * eps ** 2 * acb.integral(fA4, 0, 200).real
    # explicit tail for the truncated upper limit: at eps=1/4 every
    # integrand is <= 1e6 * t^6 * e^{-t^2/16}; int_200^inf 1e6 t^6 e^{-t^2/16} dt
    # < 1e-300 (the integrand at t=200 is ~ e^{-2500}); add it rigorously.
    tail_I4 = arb(10) ** (-300)
    I4 = T2 + main + A1 + A3 + A4 + tail_I4
    sinv = sum(1 / arb(p).log() ** 2 for p in primes)    # sum_{p<=53} (log p)^{-2}
    Gw = I4 / (2 * arb.pi()) * sinv
    return I4, Gw


def certify_Gw_uniform(eps_hi="1/20"):
    """Certify  G_w(eps) <= GW_ARCH_BOUND  uniformly for all eps in
    (0, eps_hi]  WITHOUT any monotonicity claim (the old 'I_4(eps) <= I_4(1/4)'
    argument is removed).  After the substitution u = eps t the residual splits as

       A_1(eps) = eps * int_0^inf (24+156u^2+112u^4+16u^6) e^{-u^2}(log(2+u/eps)+8) du,
       A_3(eps) = 18  * int_0^inf (10u+4u^3) e^{-u^2} [tau S_1(tau)]  du,   tau=u/2eps,
       A_4(eps) = 96  * int_0^inf  u      e^{-u^2} [tau^3 S(tau)]   du.

    (i) A_1 is INCREASING in eps on (0,eps_hi]:
        d/deps[ eps(log(2+u/eps)+8) ] = log(2+u/eps)+8 - 1/(1+2eps/u) >= 8-1 > 0,
        so  A_1(eps) <= A_1(eps_hi)  (one rigorous point evaluation).
    (ii) tau S_1(tau) <= M_1 := (9/(16 sqrt 3)) psi'(1/4)  for all tau>=0
         (sum of the per-term maxima: max_tau tau a_n/(a_n^2+tau^2)^2 = 9/(16 sqrt3 a_n^2)
          at tau=a_n/sqrt3, summed gives (9/(16 sqrt3)) sum a_n^{-2} = (9/(16 sqrt3)) psi'(1/4)),
         hence  A_3(eps) <= 18 * 7 * M_1 = 126 M_1  (int(10u+4u^3)e^{-u^2}=7), all eps.
    (iii) tau^3 S(tau) <= M_3 := (1/8) psi'(1/4)  for all tau>=0
         (max_tau tau^3 a_n/(a_n^2+tau^2)^3 = 1/(8 a_n^2) at tau=a_n), hence
         A_4(eps) <= 96 * (1/2) * M_3 = 48 M_3  (int u e^{-u^2}=1/2), all eps.

    Therefore I_4(eps) <= T_2(eps_hi)+main+A_1(eps_hi)+126 M_1+48 M_3 for all eps<=eps_hi,
    giving a single eps-uniform G_w bound.  Returns (I4_uniform_upper, Gw_uniform_upper)."""
    num, den = eps_hi.split("/"); eh = arb(int(num)) / arb(int(den))
    psi1 = arb.pi() ** 2 + 8 * arb.const_catalan()
    sqrtpi = arb.pi().sqrt()
    M1 = arb(9) / (16 * arb(3).sqrt()) * psi1            # sup_tau tau S_1(tau)
    M3 = psi1 / 8                                        # sup_tau tau^3 S(tau)
    T2 = 1920 * sqrtpi * eh
    main = (9 * arb.pi() + 6 * arb.pi() + arb("3.75") * arb.pi()) * psi1
    def fA1u(z, _):
        u = z
        return (24 + 156 * u ** 2 + 112 * u ** 4 + 16 * u ** 6) * (-u ** 2).exp() * ((2 + u / eh).log() + 8)
    A1 = eh * acb.integral(fA1u, 0, 18).real + arb(10) ** (-100)   # u-tail at u=18 is ~e^{-324}
    A3 = 126 * M1
    A4 = 48 * M3
    I4u = T2 + main + A1 + A3 + A4
    sinv = sum(1 / arb(p).log() ** 2 for p in primes)
    Gwu = I4u / (2 * arb.pi()) * sinv
    return I4u, Gwu



print("=" * 64)
print("  cert_paper8.py  --  rigorous Arb interval certificate")
print(f"  backend: python-flint {flint.__version__} (FLINT 3 / Arb ball arithmetic)")
print(f"  working precision: {PREC_BITS} bits (~{int(PREC_BITS*0.301)} digits), directed rounding")
print(f"  kappa={KAPPA}, archimedean bound G_w<= {GW_ARCH_BOUND}, tail cutoff |.|<{TAIL_CUTOFF}")
print("  covers [1e-9, 0.05]; eps in (0, 1e-9) is handled by the closed-form theorem")
print("=" * 64)
print(f"  S_3(53) lower bound = {float(arb(S3).lower()):.6f}")

# Certify the archimedean constant: I_4(1/4) <= 3561.1 ==> G_w <= 2693,
# as Arb OUTPUTS (2693 is no longer a numerical input).
_self_test_majorants()        # confirm maj >= true before trusting certify_I4
I4_cert, Gw_cert = certify_I4()
print(f"  I_4(1/4) upper bound (Arb)   = {arb(I4_cert).upper().str(7)}   (<= 3561.1 required)")
print(f"  G_w(1/4) upper bound (Arb)   = {arb(Gw_cert).upper().str(7)}   (<= {GW_ARCH_BOUND} required)")
if not (float(arb(I4_cert).upper()) <= 3561.1):
    fail("I_4(1/4) exceeds 3561.1")
if not (float(arb(Gw_cert).upper()) <= GW_ARCH_BOUND):
    fail(f"certified G_w exceeds GW_ARCH_BOUND={GW_ARCH_BOUND}")

# Certify G_w(eps) <= GW_ARCH_BOUND uniformly on (0,0.05] -- no monotonicity.
I4_uni, Gw_uni = certify_Gw_uniform("1/20")
print(f"  I_4(eps) uniform bound (Arb) = {arb(I4_uni).upper().str(7)}   (all eps in (0,0.05])")
print(f"  G_w(eps) uniform bound (Arb) = {arb(Gw_uni).upper().str(7)}   (<= {GW_ARCH_BOUND} required, no monotonicity)")
if not (float(arb(Gw_uni).upper()) <= GW_ARCH_BOUND):
    fail(f"uniform G_w(eps) exceeds GW_ARCH_BOUND={GW_ARCH_BOUND}")

# (6) coverage check: intervals chain without gap/overlap (Arb-rigorous equality of endpoints).
def rat(s):
    num, den = s.split("/"); return arb(int(num)) / arb(int(den))

prev_sb = None
for (sa, sb) in COVER:
    if prev_sb is not None and sa != prev_sb:
        fail(f"coverage gap/overlap at {sa}: previous right end {prev_sb} != this left end")
    prev_sb = sb
if COVER[-1][1] != "1/20":
    fail("coverage does not reach eps=0.05")

# (5) per-interval certificate
THR = arb(10) ** (-380)
for (sa, sb) in COVER:
    a, b = rat(sa), rat(sb)
    LB = cross_LB(sa, sb)
    tail = tail_bound(sb)
    rhs = 8 * arb.pi().sqrt() * b ** 3 * GW_ARCH_BOUND
    margin = S3 + LB - tail - rhs
    lo_LB   = float(arb(LB).lower())
    up_rhs  = float(arb(rhs).upper())
    lo_marg = float(arb(margin).lower())
    tail_ok = float((tail / THR).upper()) < 1.0          # rigorous: tail < 1e-380 (ratio avoids underflow)
    print(f"  interval ({sa}, {sb}] :  LB>= {lo_LB:10.4f}   "
          f"RHS<= {up_rhs:8.4f}   tail< 1e-380[{tail_ok}]   margin>= {lo_marg:8.4f}")
    if not (lo_marg > 0):
        fail(f"non-positive margin on ({sa}, {sb}]")
    if not tail_ok:
        fail(f"tail bound not below 1e-380 on ({sa}, {sb}]")

# (4) explicit tail bound report (Arb comparison; the value underflows float)
tb = tail_bound("1/20")
print(f"  discarded-tail bound at eps=0.05: {arb(tb).upper().str(4)}  "
      f"(rigorously < 1e-380: {float((tb / THR).upper()) < 1.0})")

# (7) hashes
src = open(__file__, "rb").read()
print(f"  source-hash (sha256, self):   {hashlib.sha256(src).hexdigest()[:32]}...")
cert_payload = f"{KAPPA}|{GW_ARCH_BOUND}|{COVER}|{PREC_BITS}".encode()
print(f"  parameter-hash  (sha256):    {hashlib.sha256(cert_payload).hexdigest()[:32]}...")

print("=" * 64)
if PASS:
    print("  CERTIFICATE VALID:  B_GW^infty(53, eps) < 0  on [1e-9, 0.05]")
    print("  (the region 0 < eps < 1e-9 is covered by the closed-form theorem,")
    print("   the closed-form theorem, which handles eps -> 0 directly).")
    print("=" * 64)
    sys.exit(0)
else:
    print("  CERTIFICATE FAILED")
    print("=" * 64)
    sys.exit(1)
