#!/usr/bin/env python3
r"""
verify_paper9.py v0.31 — numerical and textual checks for Paper 9.

Paper 9: Scaling Laws, Class Invariance, and the Limits of Audibility
         of an Explicit Prime-Zero Coupling Operator
All normative parameters: kappa=53, eps=0.05, N=100, sigma=0.5

Paper-specific gate for Paper 9 v0.31 (August 2026): 277 checks, of which the
script verifies its own count and exits non-zero on any shortfall, so a partial
run cannot end in success.  It is a drift gate over the quoted values — it
recomputes the anchors marked for recomputation and does not itself reproduce
the experiments.

What is checked (grouped):
  A. Anchors and identities — B, O''(1/2) = -2B, O'(1/2), D_SEL, the two
     energy-asymmetry values, W_eps_N, and the kernel identity of Theorem 2.2,
     each recomputed from the normative ordinate list.
  B. Leading profile — <f>, the half-width, the outward-rounded band and the
     exact positivity envelope of Lemma 3.3.
  C. Certified quantities — the two-point separation Delta_0, the curvature
     margins, the bounded-Lipschitz distance and the finite ordering of
     |E_pi|: registered as text anchors, produced by the named certificate
     artefacts.
  D. Data contracts — the normative layer has exactly 100 rows; the first 100
     rows of the extended list are binary64-identical to it.
  E. Source-level contracts of the named artefacts — each certificate script
     must carry its own input contract and its fail-closed guards.
  F. Textual anchors and prohibitions — status headings, scope sentences,
     measurement wording, rounding direction, and the forbidden-pattern list.

Two kinds of check are performed and they are reported separately:

  [R]  RECOMPUTED — the value is recomputed from the normative ordinate list
       and first principles, then compared against the value printed in the
       manuscript.  A pass means the manuscript agrees with the arithmetic.

  [T]  TEXT-ANCHOR — the value is not recomputable inside this script (it is
       a certified interval-arithmetic result, or it belongs to a construction
       archived elsewhere).  The check verifies that the manuscript carries
       exactly the registered value, so that silent drift is impossible.

Conventions used here are those of the manuscript:
    Phi(sigma)_{k,p} = e^{-eps^2 g_k^2/2} sin(sigma * g_k * log p)
    O(sigma)         = 1/2 sum_{k,p} e^{-eps^2 g_k^2} cos(2 sigma g_k log p)
    B                = sum_{k,p} e^{-eps^2 g_k^2} (g_k log p)^2 cos(g_k log p)
CAUTION: the two energy-asymmetry functionals eta live in the earlier convention with
phasor argument g_k*log p and NO trace parameter.  Mixing the conventions
changes eta by roughly 0.15 and is the single most likely source of a false
alarm; the check below therefore pins the convention explicitly.

Scope of the prohibitions.  The forbidden-pattern list is LOADED from
paper9_contract.json rather than duplicated here, so that the gate enforces the
normative contract and not a private opinion.  That file is a generated export
of the working notation contract: it carries the pattern list and nothing else,
so no internal working record is published in order to run this gate.  Three groups are added on top,
each with its source: (a) prohibitions required by the sprint brief; (b)
prohibitions arising from the first external review; (c) one local
typographic convention, marked as such.  No check in this script mandates a
replacement phrase: a gate that forces wording produces overclaims, so it may
only forbid what the contract forbids.

Figure: the block that writes figures/paper9/fig_paper9_main.png emits a status
line, never a recorded check, so the check count is identical with and without
matplotlib.

Usage:  python3 verify_paper9.py [--tex paper9_v0_31.tex] [--zeros650 zeros_650.csv] [--zeros zeros_100.csv]
                                 [--contract paper9_contract.json]
"""

import argparse, csv, json, math, re, sys
import struct
import numpy as np

# ---------------------------------------------------------------- parameters
KAPPA, EPS, N_ZEROS, SIGMA0 = 53, 0.05, 100, 0.5

# Registered anchors.  Source of record: the project's single source of truth;
# every one of them was re-pulled at the start of this sprint, not copied.
ANCHORS = {
    "B":            -19342.5476,
    "O2_half":      +38685.0952,          # = -2B
    "O1_half":      +2.4751,
    "D_SEL":         10.9845,
    "eta_orig":       0.69078176,         # weight sqrt(V_p(1/2))
    "eta_ren":        0.66926873,         # weight sqrt(f_p)
    "Delta_0":        0.0840112,
    "R_u1":           0.51690911320,
    "R_u2":           0.43289790300,
    "u1":            17.6170,
    "u2":             7.1815,
}
# Values that are certified elsewhere or belong to archived constructions.
TEXT_ANCHORS = [
    ("5343.90",                 "certified margin vs smooth control world (downward, E-6)"),
    ("5416.28",                 "certified margin vs discrete control world"),
    ("152.529455",              "bounded-Lipschitz distance d_BL"),
    ("32.78",                   "Lipschitz threshold, per-window form"),
    ("35.51",                   "Lipschitz threshold, sharp uniform form"),
    ("-241.150637407723",       "smooth-world curvature at eps=0.05"),
    ("-458.627462650173",       "discrete control curvature at eps=0.05"),
    ("-19342.5476062333843",    "certified ball value of B[P](0.05)"),
    ("0.0021613227",            "class constant of the mean law"),
    ("1.04749",                 "sideband ratio, predicted"),
    ("0.630",                   "max D/S at unfolded scale 0.25"),
    ("0.799",                   "max D/S at unfolded scale 0.10"),
    ("7.1",                     "effective sample at unfolded scale 0.25"),
    ("17.6",                    "effective sample at unfolded scale 0.10"),
    ("146521",                  "last cutoff with positive energy asymmetry"),
    ("2.184690",               "count defect of the smooth control world"),
    ("3.184690",               "sup of the counting error term"),
    ("1.136543",               "exact mass wedge M/16"),
    ("0.065",                   "resolution threshold in log distance"),
    ("0.977",                   "fusion correlation 31 with 32"),
    ("0.908",                   "fusion correlation 16 with 17"),
    ("0.907",                   "fusion correlation 47 with 49"),
    ("0.0840112102",            "value returned by the ratio certificate"),
    ("0.0020938741",            "lower end of the leading-profile band"),
    ("0.0020849161",            "exact positivity envelope of the profile"),
    ("0.0022287713",            "upper end of the leading-profile band"),
    ("6.7448532",               "half-width of the leading-profile band"),
    ("0.0021613227",            "class constant of the mean law"),
    ("18.1846895",              "mass of the smooth measure"),
    ("141.6150",                "mass of the prime measure"),
    ("139.0471",                "mass of the discrete control measure"),
    ("123.8",                   "pointwise threshold at eps = 0.05"),
    ("283.8",                   "pointwise threshold at eps = 0.04"),
    ("146527",                  "first prime cutoff with negative quotient"),
    ("0.6380807",               "uniform lower bound on the limiting trace"),
    ("640.995",                 "scatter scale of the first observable"),
    ("10.8",                    "curvature shift, weight (log n)^2, per cent"),
    ("3.8",                     "curvature shift, weight Lambda(n) log n"),
    ("1.3",                     "curvature shift, weight Lambda(n)^2"),
]
def load_contract_patterns(path):
    """E.2: the contract is the source of truth for prohibitions."""
    try:
        contract = json.load(open(path, encoding="utf-8"))
    except FileNotFoundError:
        return None
    out = []
    for e in contract.get("forbidden_patterns", []):
        if isinstance(e, dict) and e.get("pattern"):
            out.append((e["pattern"],
                        f"contract [{e.get('severity','?')}] "
                        f"{str(e.get('reason',''))[:60]}"))
    return out


# Additive to the contract: (a) sprint brief, (b) external review,
# (c) one local typographic convention.
FORBIDDEN = [
    (r"\bW[012]\b",                     "bare weight label (formulas only)"),
    (r"\b(INVARIANT|KANON|RENORM|SMOOTH-A|G2|EF-LENS)\b",
                                        "internal programme name"),
    (r"UNI-\u039b|UNI-Lambda",          "internal sprint name"),
    (r"\bSSOT\b|\bMPI\b|AnalysisLab|Sprint\s+[A-Z0-9]", "internal document label"),
    (r"P3-lite",                        "internal convention label"),
    (r"Schlussstein|Drei-Achsen|FG-Schranke", "internal German label"),
    (r"\\epsilon(?![a-z])",             "(c) local convention: \\varepsilon"),
    (r"\bK1[89]\b|\bK20\b",             "internal codename"),
    (r"blind to zero locations",
     "(b) undirected barrier phrasing; the barrier concerns orientation and "
     "the vanishing first-order response"),
    (r"independent of \$N\$|N-independent|exactly \$N\$-independent",
     "(b) exact independence of the ordinate count is false"),
]

CONTRACT_PATH = ["paper9_contract.json"]
EXPECTED_CHECKS = 277  # full-run count for v0.31; update on any deliberate change
results = []
def _load_650_prefix(n=100, cli_path=None):
    """Locate zeros_650.csv through the same candidate cascade as the
    normative list (script directory, then CWD) and return its first n
    ordinates as float, or None if the file is absent (fail-closed)."""
    import os
    _sd = os.path.dirname(os.path.abspath(__file__))
    cands = [cli_path] if cli_path else [
        os.path.join(_sd, "..", "..", "data", "zeros_650.csv"),
        os.path.join(_sd, "data", "zeros_650.csv"),
        os.path.join("data", "zeros_650.csv"),
        os.path.join(_sd, "zeros_650.csv"), "zeros_650.csv"]
    for c in cands:
        if c and os.path.exists(c):
            with open(c) as fh:
                first = fh.readline().split(",")[0].strip()
            skip = 0 if first.lstrip("-").replace(".", "", 1).isdigit() else 1
            data = np.loadtxt(c, delimiter=",", skiprows=skip)
            g = data[:, -1] if data.ndim == 2 else data.flatten()
            return np.asarray(g[:n], float), c
    return None, cands


def record(kind, name, ok, detail=""):
    results.append((kind, name, ok, detail))
    mark = "PASS" if ok else "FAIL"
    print(f"  [{kind}] {mark}  {name}" + (f"   — {detail}" if detail else ""))

def close(a, b, tol):
    return abs(a - b) <= tol

# ------------------------------------------------------------------ ordinates

def locate_and_load_zeros(cli_path, n_expected=100):
    """Series-standard ordinate loader (harmonised with Papers 1-8).

    Candidate cascade: an explicitly given path is authoritative — if it is
    missing the run stops (no silent fallback to a different file); only the
    default cascades through script directory and working directory. Header
    sniffing as in verify_paper8: if the first cell of the first line is
    numeric the file is headerless, otherwise one line is skipped — the
    normative zeros_100.csv historically shipped without a header, and
    zeros_650.csv ships with one."""
    import os
    explicit = cli_path is not None and cli_path != ""
    if explicit:
        candidates = [cli_path]
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
            data = np.loadtxt(c, delimiter=",", skiprows=skip)
            g = data[:, -1] if data.ndim == 2 else data.flatten()
            print(f"  ordinates: {len(g)} from {c} "
                  f"(header {'skipped' if skip else 'absent'})")
            return np.asarray(g[:n_expected], float)
    print(f"\n  FATAL: ordinate list not found; tried {candidates}")
    sys.exit(2)


def load_zeros(path):
    return locate_and_load_zeros(path)


def primes_upto(n):
    return [p for p in range(2, n + 1)
            if all(p % d for d in range(2, int(p ** 0.5) + 1))]

# ------------------------------------------------------------------- checks
def check_numeric(gam, zeros650=None):
    print("\n=== recomputed anchors ===")
    P = primes_upto(KAPPA)
    record("R", f"prime count at cutoff {KAPPA} is 16", len(P) == 16, f"{len(P)}")
    record("R", f"ordinate list has {N_ZEROS} entries", len(gam) == N_ZEROS,
           f"{len(gam)}")

    # cycle-24: the prose promise "its first 100 entries reproduce the
    # normative layer exactly in double precision" is now gated, not just
    # asserted. Missing file = FAIL, never a silent skip.
    _g650, _src650 = _load_650_prefix(100, zeros650)
    if _g650 is None:
        record("R", "zeros_650 first 100 rows are binary64-identical to "
               "the normative layer", False, f"file not found; tried {_src650}")
    else:
        _g100 = np.asarray(gam[:100], float)
        _bad = [k for k in range(min(100, len(_g650), len(_g100)))
                if struct.pack("<d", _g650[k]) != struct.pack("<d", _g100[k])]
        record("R", "zeros_650 first 100 rows are binary64-identical to "
               "the normative layer",
               len(_g650) == 100 and len(_bad) == 0,
               f"{100-len(_bad)}/100 identical" if not _bad
               else f"differs at k={[b+1 for b in _bad][:5]}")

    lp = np.log(np.array(P, float))
    w = np.exp(-EPS ** 2 * gam ** 2)
    u = np.outer(gam, lp)

    B = float(np.sum(w[:, None] * u ** 2 * np.cos(u)))
    record("R", "B at reference parameters", close(B, ANCHORS["B"], 5e-4),
           f"{B:.4f}")

    O2 = -2.0 * B
    record("R", "curvature identity O''(1/2) = -2B",
           close(O2, ANCHORS["O2_half"], 1e-3), f"{O2:+.4f}")

    O1 = -float(np.sum(w[:, None] * u * np.sin(u)))
    record("R", "O'(1/2) sign and value", close(O1, ANCHORS["O1_half"], 5e-4),
           f"{O1:+.4f}")

    # trace identity  Tr T~ = D_SEL - O(1/2)
    wk = np.exp(-EPS ** 2 * gam ** 2)
    Phi = np.exp(-EPS ** 2 * gam ** 2 / 2)[:, None] * np.sin(SIGMA0 * u)
    Tr = float(np.sum(Phi ** 2))
    D_SEL = len(P) / 2 * float(np.sum(wk))
    O_half = 0.5 * float(np.sum(wk[:, None] * np.cos(u)))
    record("R", "D_SEL constant", close(D_SEL, ANCHORS["D_SEL"], 5e-4),
           f"{D_SEL:.6f}")
    record("R", "trace identity Tr = D_SEL - O(1/2)",
           close(Tr, D_SEL - O_half, 1e-9), f"residual {Tr-(D_SEL-O_half):.2e}")

    # kernel identity  T_pq = 1/2 [ G(log p/q) - G(log pq) ]
    G = lambda v: float(np.sum(wk * np.cos(SIGMA0 * gam * v)))
    T = Phi.T @ Phi
    err = max(abs(T[i, j] - 0.5 * (G(lp[i] - lp[j]) - G(lp[i] + lp[j])))
              for i in range(len(P)) for j in range(len(P)))
    record("R", "kernel identity is exact", err < 1e-12, f"max dev {err:.2e}")

    # eta pair — EARLIER convention: phasor argument g_k log p, no sigma
    Phi_eta = (np.exp(-EPS ** 2 * gam ** 2 / 2)[:, None]) * np.sin(u)
    def eta(cf):
        c = np.array([math.sqrt(cf(p)) for p in P])
        E_str = float(np.sum(c ** 2 * np.sum(Phi_eta ** 2, axis=0)))
        v = Phi_eta @ c
        return 1.0 - float(v @ v) / E_str
    V = lambda p: 4 * math.log(p) ** 2 * p / (p - 1) ** 2
    f = lambda p: 4 * math.log(p) ** 2 * (2 * p - 1) / (p * (p - 1) ** 2)
    eo, er = eta(V), eta(f)
    record("R", "eta_orig with weight sqrt(V_p(1/2))",
           close(eo, ANCHORS["eta_orig"], 5e-8), f"{eo:.8f}")
    record("R", "eta_ren with weight sqrt(f_p)",
           close(er, ANCHORS["eta_ren"], 5e-8), f"{er:.8f}")
    record("R", "the two eta values are distinct", abs(eo - er) > 0.02,
           f"difference {eo-er:.8f}")

    # two-point separation of the limiting ratio
    wk2 = np.exp(-EPS ** 2 * gam ** 2 / 2)
    def h(th, uu): return (np.cos(th * uu) + th * np.sin(th * uu)) / (1 + th * th)
    def R_of(uu):
        Q = np.outer(wk2, wk2) / 2 * (
            h((gam[:, None] - gam[None, :]) / 2, uu)
            - h((gam[:, None] + gam[None, :]) / 2, uu))
        np.fill_diagonal(Q, wk2 ** 2 / 2 * (1 - h(gam, uu)))
        return float(np.linalg.eigvalsh(Q)[-1] / np.trace(Q))
    r1, r2 = R_of(ANCHORS["u1"]), R_of(ANCHORS["u2"])
    record("R", "limiting ratio at the upper certification point",
           close(r1, ANCHORS["R_u1"], 1e-9), f"{r1:.11f}")
    record("R", "limiting ratio at the lower certification point",
           close(r2, ANCHORS["R_u2"], 1e-9), f"{r2:.11f}")
    record("R", "two-point separation Delta_0",
           r1 - r2 >= ANCHORS["Delta_0"] - 1e-7, f"{r1-r2:.7f}")

    # class constants of the mean law and of the sideband formula
    fm = 0.5*float(np.sum(w/(0.25 + gam**2)))
    record("R", "mean of the ordinate profile", close(fm, 0.0021613227, 5e-10),
           f"{fm:.10f}")
    aa, om = 0.3, math.pi*2/3
    band = (aa**2/8)*float(np.sum(w*(1/(0.25+(gam+om)**2)+1/(0.25+(gam-om)**2))))
    record("R", "sideband ratio from the closed formula",
           close((fm+band)/fm, 1.04749, 5e-6), f"{(fm+band)/fm:.6f}")
    W = float(np.sum(np.exp(-EPS**2*gam**2)))
    record("R", "W, the sum of squared damping weights",
           close(W, 1.37306, 5e-5), f"{W:.5f}")
    tau0 = W/2*(1 - 1/math.sqrt(1 + gam[0]**2))
    record("R", "uniform lower bound tau_0 on the limiting trace",
           close(tau0, 0.6380807, 5e-7), f"{tau0:.7f}")
    from mpmath import mp, quad, log as mlog
    mp.dps = 25
    Msm = float(quad(lambda t: 1/mlog(t), [2, KAPPA]))
    record("R", "mass M of the smooth measure", close(Msm, 18.1846895, 5e-7),
           f"{Msm:.7f}")
    record("R", "count defect |16 - M|", close(abs(16-Msm), 2.184690, 5e-6),
           f"{abs(16-Msm):.6f}")
    record("R", "mass wedge M/16", close(Msm/16, 1.136543, 5e-6),
           f"{Msm/16:.6f}")
    tail = float(np.sum(np.exp(-EPS**2*gam[50:]**2)))
    Phi50 = np.exp(-EPS**2*gam[:50]**2/2)[:, None]*np.sin(
        SIGMA0*np.outer(gam[:50], lp))
    bnd = 2*len(P)*tail/float(np.sum(Phi50**2))
    record("R", "ordinate tail bound at N = 50 is of order 1e-23",
           1e-24 < bnd < 1e-22, f"{bnd:.3e}")

    # near-coincidence fusions on the completed mode set
    modes = sorted(n for n in range(2, KAPPA + 1)
                   if any(n == p ** m for p in P for m in range(1, 7)))
    record("R", "completed mode set has 24 members", len(modes) == 24,
           f"{len(modes)}")
    ln = np.log(np.array(modes, float))
    Ph = np.exp(-EPS ** 2 * gam ** 2 / 2)[:, None] * np.sin(
        SIGMA0 * np.outer(gam, ln))
    Gm = Ph.T @ Ph
    d = np.sqrt(np.diag(Gm)); rho = Gm / np.outer(d, d)
    ix = {n: i for i, n in enumerate(modes)}
    for a, b, target in [(31, 32, 0.977), (16, 17, 0.908), (47, 49, 0.907)]:
        v = rho[ix[a], ix[b]]
        record("R", f"fusion correlation {a} with {b}", close(v, target, 5e-4),
               f"{v:+.3f}")
    for a, b in [(31, 32), (16, 17), (47, 49)]:
        dl = abs(math.log(a) - math.log(b))
        record("R", f"fused pair {a},{b} lies below the 0.065 threshold",
               dl < 0.065, f"dlog {dl:.4f}")
    for a, b, target in [(25, 27, 0.540), (8, 9, 0.224)]:
        v = rho[ix[a], ix[b]]
        record("R", f"non-fusing pair {a},{b} stays below 0.9",
               v < 0.9 and close(v, target, 5e-3), f"{v:+.3f}")

def check_text(tex, gam):
    print("\n=== manuscript anchors and prohibitions ===")
    body = tex
    # strip comments so that editorial placeholders are not scanned as body,
    # then normalise whitespace: LaTeX line breaks must not hide a phrase.
    body_nc = re.sub(r"(?m)^\s*%.*$", "", body)
    body_nc = re.sub(r"\s+", " ", body_nc)

    for key, val in [("B", "-19342.5476"), ("O2_half", "38685.0952"),
                     ("O1_half", "2.4751"), ("D_SEL", "10.9845"),
                     ("eta_orig", "0.69078176"), ("eta_ren", "0.66926873"),
                     ("Delta_0", "0.0840112")]:
        record("T", f"manuscript carries anchor {key}", val in body_nc, val)

    for val, what in TEXT_ANCHORS:
        record("T", f"manuscript carries {what}", val in body_nc, val)

    contract = load_contract_patterns(CONTRACT_PATH[0])
    if contract is None:
        record("T", "pattern contract loaded", False, "paper9_contract.json missing")
    else:
        record("T", f"pattern contract loaded ({len(contract)} patterns)", True)
        bad = []
        for pat, why in contract:
            try:
                if re.search(pat, body_nc):
                    bad.append(why)
            except re.error:
                pass
        record("T", "no contract-forbidden pattern", not bad,
               "; ".join(bad[:3]) if bad else "")

    for pat, why in FORBIDDEN:
        hits = re.findall(pat, body_nc)
        record("T", f"no {why}", len(hits) == 0,
               f"found {sorted(set(hits))[:4]}" if hits else "")

    # E.3 — measurement wording and artefact naming
    record("T", "measurement outcome avoids the word 'confirmed'",
           not re.search(r"expectation[^.]{0,60}\bconfirmed\b", body_nc))
    record("T", "measurement framing avoids test vocabulary",
           "not rejected" not in body_nc and "null expectation" not in body_nc)
    record("T", "measurement framing carries delimitation and registration",
           "does not test an $n$-level" in body_nc and
           "none of the five" in body_nc and
           "separates the declared worlds" in body_nc)
    record("T", "the five observables are given by formula",
           body_nc.count("\\mathcal{O}_") >= 5)
    record("T", "the canonicity proposition is qualified by its hypotheses",
           "under (H1)--(H2)" in body_nc)
    check_cycle3(body, body_nc, gam)
    record("T", "the Beurling-space delineation is cited",
           "Kouroupis23" in body_nc)
    m = re.search(r"\\section\{Methods and Verification\}(.*?)"
                  r"\\section\{Open Problems\}", body, re.DOTALL)
    meth = re.sub(r"\s+", " ", m.group(1)) if m else ""
    for art in ["verify\\_paper9.py", "cert\\_paper9\\_ratio.py",
                "certify\\_smooth\\_controls.py", "gen\\_zeros.py", "python-flint", "mpmath"]:
        record("T", f"methods section names {art.replace(chr(92), '')}",
               art in meth)
    n_cert = len(re.findall(r"; certified\]", body_nc))
    n_prot = body_nc.count("\\emph{Certification.}")
    record("T", "each certified heading has a certification paragraph",
           n_cert <= n_prot, f"{n_cert} headings, {n_prot} protocols")

    # structural checks
    record("T", "Open Problems is the last numbered section",
           body.rfind(r"\section{Open Problems}") >
           max(body.rfind(r"\section{" + s) for s in
               ["Introduction", "The Coupling Operator", "Scaling Laws",
                "Class Invariance", "The Certified", "Limits of",
                "The Audibility", "Methods and"]))
    record("T", "Non-Circularity follows Open Problems",
           body.find(r"\section*{Non-Circularity}") >
           body.find(r"\section{Open Problems}"))
    record("T", "symmetry barrier is stated before the measurement",
           body.find("thm:barrier") < body.find("Outcome (numerical)"))
    record("T", "no theorem status appears in upper case brackets",
           not re.search(r"\[(PROVED|NUMERICAL|CERTIFIED|CONDITIONAL|OPEN)\]",
                         body_nc))
    record("T", "every certified claim names a protocol",
           body_nc.count("certified") > 0 and "Certification." in body_nc)
    record("T", "no forward citation beyond Paper 8",
           not re.search(r"\\cite\{Paper(9|1[0-9])\}", body_nc))
    record("T", "the five required 2026 delineation citations are present",
           all(k in body_nc for k in
               ["Connes26", "ConnesVS25", "Groskin26a", "Groskin26b",
                "Suzuki26"]))
    record("T", "measurement carries the explicit no-go disclaimer",
           "No no-go theorem is derived" in body_nc)
    record("T", "measurement scope names battery, window and coupling",
           "this battery, this window and this coupling" in body_nc)
    record("T", "pre-registration is framed as a methods import",
           "methodological import" in body_nc)

def _figure_dir():
    """Repo layout first (figures/paperN next to code/ and data/), then the
    script directory; mirrors the ordinate cascade."""
    import os
    sd = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.join(sd, "..", "..")
    if os.path.isdir(os.path.join(repo, "data")):
        return os.path.normpath(os.path.join(repo, "figures", "paper9"))
    if os.path.isdir("data"):
        return os.path.join("figures", "paper9")
    return os.path.join(sd, "figures", "paper9")


def make_figure(gam):
    """Four panels, all from quantities the manuscript already states.

    Wrapped in try/except so that a missing matplotlib is a printed notice,
    never a gate failure: the check count is invariant to this block.
    """
    import os
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpmath import mp, mpf, li as _li
        mp.dps = 30
        outdir = _figure_dir()
        os.makedirs(outdir, exist_ok=True)
        fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))

        # (a) finite ordering of |E_pi| — Thm 5.2 (statement and proof)
        Pr = primes_upto(KAPPA)
        li2 = _li(mpf(2))
        cands = []
        for i, q in enumerate(Pr):
            cands.append((f"{q}+", abs(mpf(i + 1) - (_li(mpf(q)) - li2))))
            if i < len(Pr) - 1:
                cands.append((f"{Pr[i+1]}-",
                              abs(mpf(i + 1) - (_li(mpf(Pr[i + 1])) - li2))))
        cands.sort(key=lambda kv: -float(kv[1]))
        vals = [float(v) for _, v in cands]
        ax[0, 0].bar(range(len(vals)), vals, color="#c7c7c7", width=0.8)
        for j, col in zip((0, 1, 2), ("#d62728", "#ff7f0e", "#1f77b4")):
            ax[0, 0].bar([j], [vals[j]], color=col, width=0.8)
        for j, col in zip((0, 1, 2), ("#d62728", "#ff7f0e", "#1f77b4")):
            ax[0, 0].hlines(vals[j], j, len(vals) - 1, color=col, lw=0.7,
                            ls=":", alpha=0.8)
        ax[0, 0].text(0.97, 0.95,
                      "certified margins\n"
                      r"rank 1$\to$2  $\geq 0.1201$" "\n"
                      r"rank 2$\to$3  $\geq 0.0920$",
                      transform=ax[0, 0].transAxes, ha="right", va="top",
                      fontsize=7.5,
                      bbox=dict(boxstyle="round,pad=0.35", fc="white",
                                ec="#999999", lw=0.6))
        ax[0, 0].set_ylim(0, vals[0] * 1.30)
        ax[0, 0].axhline(2.184690, color="k", ls=":", lw=0.8)
        ax[0, 0].text(len(vals) * 0.55, 2.23, r"right endpoint $|16-M|=2.184690$",
                      fontsize=7)
        ax[0, 0].set_title(r"(a) certified ordering of $|E_\pi|$: "
                           r"31 candidates, 30 directed links", fontsize=9)
        ax[0, 0].set_xlabel("rank"); ax[0, 0].set_ylabel(r"$|E_\pi|$")

        # (b) leading profile: band and exact envelope — Lemma 3.3
        ax[0, 1].axhspan(0.0020938741, 0.0022287713, color="#1f77b4", alpha=0.18)
        for y, lab, st in ((0.0022287713, "upper band edge  0.0022287713", "-"),
                           (0.0021613227, r"$\langle f\rangle$  0.0021613227", "--"),
                           (0.0020938741, "lower band edge  0.0020938741", "-"),
                           (0.0020849161, "exact envelope  0.0020849161", "-.")):
            ax[0, 1].axhline(y, color="#1f77b4" if st != "-." else "#d62728",
                             ls=st, lw=1.1)
            ax[0, 1].text(0.02, y, lab, fontsize=7, va="bottom")
        ax[0, 1].set_ylim(0.002080, 0.002233); ax[0, 1].set_xticks([])
        ax[0, 1].set_title("(b) leading profile: band and exact positivity "
                           "envelope", fontsize=9)

        # (c) curvature comparison and certified margins — Thm 5.4 / 5.5
        names = [r"$B[P]$", r"$B_{\rm sm}$ margin", r"$B[A]$ margin"]
        heights = [abs(ANCHORS["B"]), 5343.90, 5416.28]
        bars = ax[1, 0].bar(names, heights,
                            color=["#333333", "#2ca02c", "#9467bd"])
        for b, h in zip(bars, heights):
            ax[1, 0].text(b.get_x() + b.get_width() / 2, h * 1.02,
                          f"{h:.2f}", ha="center", fontsize=8)
        ax[1, 0].set_yscale("log")
        ax[1, 0].set_title(r"(c) $|B[P]|$ and the two certified margins "
                           r"(downward-rounded)", fontsize=9)

        # (d) five-world battery, D_i/S_i — section 6 table
        M = [[0.299, 0.069, 0.055, 0.580, 0.135],
             [0.342, 0.031, 0.064, 0.324, 0.073],
             [0.154, 0.630, 0.586, 0.270, 0.270],
             [0.082, 0.060, 0.018, 0.168, 0.012]]
        im = ax[1, 1].imshow(M, cmap="viridis", vmin=0, vmax=3.0, aspect="auto")
        ax[1, 1].set_xticks(range(5))
        ax[1, 1].set_xticklabels([r"$\mathcal{O}_%d$" % (i + 1) for i in range(5)])
        ax[1, 1].set_yticks(range(4))
        ax[1, 1].set_yticklabels([r"$W_\zeta$", r"$W_{\chi_4}$",
                                  r"$W_{\rm Ded}$", r"$W_{\rm DH}$"])
        for i in range(4):
            for j in range(5):
                ax[1, 1].text(j, i, f"{M[i][j]:.3f}", ha="center", va="center",
                              color="w", fontsize=7)
        fig.colorbar(im, ax=ax[1, 1], label="separation threshold 3")
        ax[1, 1].set_title(r"(d) $D_i/S_i$: largest entry 0.630, "
                           r"a factor 4.8 below 3", fontsize=9)

        fig.tight_layout()
        path = os.path.join(outdir, "fig_paper9_main.png")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"\n  figure written to {path}")
    except ImportError:
        print("\n  figure skipped (matplotlib not installed) "
              "— gate result unaffected")
    except Exception as exc:                     # pragma: no cover
        print(f"\n  figure skipped ({exc}) — gate result unaffected")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", default=None)
    ap.add_argument("--zeros", default=None)
    ap.add_argument("--zeros650", default=None)
    ap.add_argument("--contract", default=None)
    a = ap.parse_args()

    print("=" * 62)
    print("verify_paper9.py v0.31")
    print("=" * 62)
    import os as _os
    _sd = _os.path.dirname(_os.path.abspath(__file__))
    _explicit = a.contract is not None
    _cands = [a.contract] if _explicit else [
        _os.path.join(_sd, "..", "..", "data", "paper9_contract.json"),
        _os.path.join(_sd, "data", "paper9_contract.json"),
        _os.path.join("data", "paper9_contract.json"),
        _os.path.join(_sd, "paper9_contract.json"), "paper9_contract.json"]
    for _c in _cands:
        if _os.path.exists(_c):
            CONTRACT_PATH[0] = _c
            break
    else:
        print(f"\n  FATAL: pattern contract not found; tried {_cands}")
        return 2
    gam = load_zeros(a.zeros)          # exits 2 itself if missing
    # manuscript: an explicit --tex is authoritative (single path, fail-
    # closed); only the default cascades — repo layout (papers/paper9/,
    # per README "Setup and Run"), then script dir, then CWD.
    if a.tex is not None:
        _tex_cands = [a.tex]
    else:
        _tex_cands = [
            _os.path.join(_sd, "..", "..", "papers", "paper9",
                          "paper9_v0_31.tex"),
            _os.path.join(_sd, "..", "..", "papers", "paper9_v0_31.tex"),
            _os.path.join(_sd, "paper9_v0_31.tex"),
            _os.path.join("papers", "paper9", "paper9_v0_31.tex"),
            "paper9_v0_31.tex"]
    for _t in _tex_cands:
        if _os.path.exists(_t):
            a.tex = _t
            break
    else:
        print(f"\n  FATAL: manuscript not found; tried {_tex_cands}")
        return 2
    tex = open(a.tex, encoding="utf-8").read()
    print(f"  manuscript: {a.tex}")

    check_numeric(gam, a.zeros650)
    check_text(tex, gam)
    make_figure(gam)

    npass = sum(1 for *_, ok, _ in results if ok)
    nfail = len(results) - npass
    print("\n" + "=" * 62)
    print(f"  {npass} PASS, {nfail} FAIL   ({len(results)} of "
          f"{EXPECTED_CHECKS} expected checks recorded)")
    for kind, name, ok, detail in results:
        if not ok:
            print(f"    FAIL [{kind}] {name}   {detail}")
    if len(results) != EXPECTED_CHECKS:
        print(f"  FATAL: check count mismatch — a block was skipped or the "
              f"expectation is stale; update EXPECTED_CHECKS on any "
              f"deliberate change")
        print("=" * 62)
        return 2
    print("=" * 62)
    return 1 if nfail else 0


def check_cycle3(body, body_nc, gam):
    """One machine watcher per HIGH finding of the second external review."""
    import numpy as np
    print("\n=== cycle-3 watchers ===")

    # D.1 — B[A](0.05) recomputed from the node definition, both readings
    from mpmath import mp, quad, log as mlog, mpf
    mp.dps = 25
    w2 = np.exp(-EPS ** 2 * gam ** 2)
    Kk = lambda t: float(np.sum(w2 * (gam * np.log(t)) ** 2 * np.cos(gam * np.log(t))))
    li = lambda a: quad(lambda t: 1 / mlog(t), [2, a])
    M = float(li(KAPPA))

    def node(target):
        lo, hi = mpf("2.0000001"), mpf("1e6")
        for _ in range(120):
            mid = (lo + hi) / 2
            if li(mid) < target:
                lo = mid
            else:
                hi = mid
        return float((lo + hi) / 2)

    B_ok = sum(Kk(node((j - 0.5) * M / 16)) for j in range(1, 17))
    B_bad = sum(Kk(node(j - 0.5)) for j in range(1, 17))
    record("R", "B[A](0.05) from the equal-mass node definition",
           close(B_ok, -458.627462650173, 1e-5), f"{B_ok:.6f}")
    record("R", "the equal-mass and the naive node definitions differ grossly",
           abs(B_bad - B_ok) > 5000, f"naive gives {B_bad:.2f}")
    record("T", "manuscript prints the equal-mass node definition",
           "\\Bigl(j-\\tfrac12\\Bigr)\\frac{M}{16}" in body
           or "(j-\\tfrac12)\\frac{M}{16}" in body)

    # D.2 — leading-profile band, recomputed, and the profile form in the text
    b = w2 / (2 * (0.25 + gam ** 2))
    half = float(np.sum(b / np.sqrt(1 + 4 * gam ** 2)))
    fm = float(np.sum(b))
    record("R", "half-width of the leading-profile band",
           close(half, 6.7448532e-5, 5e-12), f"{half:.10e}")
    record("R", "lower end of the leading-profile band (one-sided: the "
           "printed decimal must not exceed the exact edge)",
           (fm - half) >= 0.0020938741, f"{fm-half:.12f}")
    _env = fm * (1.0 - 1.0 / math.sqrt(1.0 + 4.0 * gam[0] ** 2))
    record("R", "the exact envelope is positive and below the band edge",
           _env > 0 and _env < fm - half and _env >= 0.0020849161,
           f"{_env:.12f}")
    record("R", "the band is bounded away from zero", fm - half > 0)
    record("T", "the divergence theorem uses the profile, not a constant rate",
           "F(\\log X)" in body)
    record("T", "no constant-rate leading term",
           not re.search(r"4\\langle f\\rangle\}\{W\}\\cdot\\frac\{X\}", body_nc))

    # D.3 — class transfer speaks of the profile
    m = re.search(r"\\section\{Class Invariance over Beurling Systems\}(.*?)"
                  r"\\section\{The Certified", body, re.DOTALL)
    cls = re.sub(r"\s+", " ", m.group(1)) if m else ""
    record("T", "class section avoids 'identical constants'",
           "identical constants" not in cls)
    record("T", "class section states profile invariance",
           "same leading profile" in cls)
    record("T", "class section carries the exact envelope",
           "exact positivity envelope" in cls)

    # D.4 — conjugate profile with the record formula; smooth measure is not
    #        called a Beurling system
    record("T", "conjugate profile defined by the record formula",
           "\\sin\\varphi u-\\varphi\\cos\\varphi u" in body)
    m = re.search(r"\\begin\{definition\}\[Smooth generalised-prime measure(.*?)"
                  r"\\end\{definition\}", body, re.DOTALL)
    sm = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "the smooth measure is not called a Beurling system",
           sm != "MISSING" and "Beurling system" not in sm)
    record("T", "the sine convention is declared", "sine convention" in body_nc)

    # D.5 — tail theorem and proof live on the prime space
    m = re.search(r"\\label\{thm:sat\}(.*?)\\end\{proof\}", body, re.DOTALL)
    tail = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "tail theorem is stated and proved on the loop operator",
           "T_\\infty-T_N" in tail and "\\Hstr" in tail)
    record("T", "tail proof does not operate on the ordinate-space operator",
           "\\Tt_\\infty-\\Tt_N" not in tail)
    record("T", "the bridge to the ordinate space is stated",
           "non-zero spectra" in tail)

    # D.6 — qualified quadratic; 'matched exactly' only near the rescaling
    for m in re.finditer(r"[^.]*quadratic[^.]*\.", body_nc):
        seg = m.group(0)
        if "response" in seg or "surviving term" in seg or "first term" in seg:
            record("T", "quadratic claim carries its amplitude qualification",
                   "orbit average" in seg or "orbit\naverage" in seg
                   or "cancel" in seg or "differentiab" in seg
                   or "real-analytic" in seg, seg[:70] + "...")
    for m in re.finditer(r"[^.]*matched exactly[^.]*\.|[^.]*matches the prime "
                         r"count exactly[^.]*\.", body_nc):
        seg = m.group(0)
        record("T", "exact-match claim is attached to the rescaling",
               "norm" in seg or "rescaling" in seg, seg[:70] + "...")

    # ---- cycle-4 watchers (third review) --------------------------------
    for m in re.finditer(r"[^.]*position of (?:a|any) zero[^.]*\.", body_nc):
        seg = m.group(0)
        record("T", "position-of-a-zero sentence carries the C3 form",
               "available" not in seg and
               ("orientation" in seg or "drawn" in seg or "made" in seg),
               seg[:70] + "...")
    m = re.search(r"Axis & Finding & Status(.*?)\\end\{tabular\}",
                  body, re.DOTALL)
    tab = m.group(1) if m else "MISSING"
    record("T", "classification row carries the point-configuration form",
           "point configuration" in tab and "stated controls" in tab)
    record("T", "no Euler attribution inside the classification table",
           "Euler" not in tab)
    m = re.search(r"\\label\{thm:sideband\}(.*?)\\end\{theorem\}",
                  body, re.DOTALL)
    sb = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "sideband theorem states the non-resonance hypothesis",
           "non-resonance" in sb and "pairwise distinct" in sb)
    record("T", "the exact-resonance remark is present",
           "Exact resonance" in body)
    for m in re.finditer(r"[^.]*146521[^.]*\.", body_nc):
        seg = m.group(0)
        record("T", "146521 appears only with its search range",
               ("2\\cdot10^{5}" in seg) or
               ("scan" in seg and "prime cutoff" in seg),
               seg[:70] + "...")
    record("T", "no unqualified sharpness claim for the Beurling exponent",
           "sharp in a precise sense" not in body_nc)
    record("T", "the majorant qualification is present",
           "majorant threshold of the present argument" in body_nc)
    from mpmath import quad as _q, log as _l
    li = lambda a: float(_q(lambda t: 1 / _l(t), [2, a]))
    Pr = primes_upto(KAPPA)
    vals = [abs(j - li(p)) for j, p in enumerate(Pr)]
    Mv = li(KAPPA)
    record("R", "finite 16-point comparison: max |E_pi| = M - 15 at 53^-",
           abs(max(vals) - (Mv - 15)) < 1e-9 and
           vals.index(max(vals)) == len(Pr) - 1,
           f"max {max(vals):.7f}, runner-up {sorted(vals)[-2]:.4f}")
    record("R", "quadrature constants in both normalisations",
           close(Mv * Mv / 384, 0.861153, 5e-7) and
           close(Mv ** 3 / 6144, 0.978738, 5e-7),
           f"{Mv*Mv/384:.6f} / {Mv**3/6144:.6f}")

    # ---- cycle-5 watchers (fourth review) -------------------------------
    m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}",
                  body, re.DOTALL)
    abst = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "abstract does not attribute audibility to multiplicative "
           "or Euler structure",
           "multiplicative structure is audible" not in abst and
           "Euler structure is audible" not in abst)
    record("T", "abstract carries the point-configuration form",
           "point configuration" in abst or "point\nconfiguration" in abst)
    record("T", "abstract carries the non-identification sentence",
           "not multiplicativity as its cause" in abst)
    m = re.search(r"\\label\{thm:eta\}(.*?)\\end\{proof\}", body, re.DOTALL)
    eta = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "Step 1 does not argue from continuity of the Dirichlet "
           "series", "extends continuously" not in eta)
    record("T", "Step 1 carries the PNT error term",
           "e^{-c\\sqrt{\\log t}}" in eta)
    record("T", "structural energy carries W(log x)^2 + O(log x)",
           "\\WeN(\\log x)^2+O(\\log x)" in eta)
    m = re.search(r"\\begin\{remark\}\[The smooth measure carries no zeros\]"
                  r"(.*?)\\end\{remark\}", body, re.DOTALL)
    zrem = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "E1 formula is attached to the infinite measure",
           "\\emph{infinite} measure" in zrem and
           "\\mu_{\\mathrm{sm}}^{\\infty}" in zrem)
    record("T", "the restriction carries no singularity claim",
           "entire function" in zrem and "no singularity claim" in zrem)
    m = re.search(r"part of\s+the\s+certificate rather than(.*?)"
                  r"\\texttt\{cert", body, re.DOTALL)
    cpar = m.group(1) if m else "MISSING"
    record("T", "certification paragraph uses R(u1)-R(u2)",
           "R(u_1)-R(u_2)" in cpar and "r(u_1)-r(u_2)" not in cpar)

    # ---- cycle-6 watchers (fifth review) --------------------------------
    record("T", "measurement Gram matrix is G^un",
           "G^{\\mathrm{un}}=\\Phi" in body)
    record("T", "no bare G = Phi^T Phi anywhere",
           not re.search(r"\$G=\\Phi", body_nc))
    m = re.search(r"\\label\{rem:lambda\}(.*?)\\end\{remark\}", body, re.DOTALL)
    rl = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "Remark 2.6 separates curvature and spectral-ratio "
           "robustness",
           "between $1.3\\%$ and $10.8\\%$" in rl and
           "separately, the two quoted spectral ratios" in rl)
    record("T", "the contradictory sub-4%-second-moment sentence is gone",
           "sub-$4\\%$ effect on the second moment" not in rl)
    record("T", "canonicity frame uses character rigidity",
           "Character rigidity under unimodularity" in body and
           "character part is rigid" in body_nc)
    record("T", "forced/modelling-choice overclaims are gone",
           "is in fact forced" not in body_nc and
           "not a modelling choice" not in body_nc)
    m = re.search(r"\\label\{op:pointwise\}(.*?)\\end\{openproblem\}",
                  body, re.DOTALL)
    op1 = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "OP 9.1 carries the two-task form",
           "splits into two tasks" in op1 and "no rate" in op1)
    record("T", "OP 9.1 makes no equivalence claim",
           "equivalent" not in op1)
    record("T", "audibility definition precedes the classification table",
           body.find("audible if some observable of the stated class")
           < body.find("Axis & Finding & Status") and
           "no absolute notion of audibility" in body_nc)
    record("T", "factorisation barrier is proved from certified inputs",
           re.search(r"Factorisation barrier: bounded-Lipschitz resolution lower\s+bound; proved from certified inputs\]", body) is not None)
    record("T", "the certified d_BL lemma is present and referenced",
           "Bounded-Lipschitz distance of the weighted measures" in body_nc
           and "as certified in Lemma~\\ref{lem:dbl}" in body_nc)
    record("T", "partial summation carries the plus sign",
           "+\\int_2^x E(t)\\,\\frac{\\log t-1}{t^2}" in body and
           "-\\int_2^x E(t)\\,\\frac{\\log t-1}{t^2}" not in body)
    record("T", "structure sentence uses the labels-as-stated form",
           "carries its status in its heading" in body_nc)
    record("T", "identification design appears in both open lists",
           body_nc.count("identification design") >= 2)
    record("T", "resolution diagnostics and fusion paragraphs are labelled",
           body_nc.count("Numerical observation") >= 4)

    # ---- cycle-7 watchers (sixth review) --------------------------------
    record("T", "the ordinate-density-audible row is gone",
           "Ordinates, density class & audible" not in body and
           "only the density class is audible" not in body_nc)
    m = re.search(r"Axis & Finding & Status(.*?)\\end\{tabular\}",
                  body, re.DOTALL)
    tab7 = m.group(1) if m else "MISSING"
    record("T", "table carries the Beurling within-class invariance row",
           "Prime-side Beurling within-class variation" in tab7 and
           "not distinguished" in tab7)
    record("T", "table carries the not-separated ordinate row",
           "beyond unfolded one-point density" in tab7 and
           "not separated on this battery" in tab7)
    rows_aud = [r for r in tab7.split("\\\\")
                if "& audible" in r and "not audible" not in r]
    record("T", "every audible table row carries certified status",
           all("certified" in r for r in rows_aud) and len(rows_aud) >= 1,
           f"{len(rows_aud)} audible row(s)")
    record("T", "abstract carries the within-class invariance form",
           "invariant across the Beurling" in body_nc)
    record("T", "Theorem 4.2 carries the M0 class",
           "\\Dclass(C,x_0,M_0)" in body and
           "(C,x_0,\\delta,M_0)" in body)
    record("T", "the unbounded uniformity clause is gone",
           "only on $(C,x_0,\\delta)$;" not in body_nc)
    m = re.search(r"\\label\{thm:class\}(.*?)\\end\{proof\}", body, re.DOTALL)
    cls14 = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "the necessity clause is gone from Theorem 4.2 (INVERTED)",
           cls14 != "MISSING" and "cannot be dropped" not in cls14)
    m = re.search(r"\\label\{op:pointwise\}(.*?)\\end\{openproblem\}",
                  body, re.DOTALL)
    op7 = m.group(1) if m else "MISSING"
    record("T", "OP 9.1 uses the Step-3 pointwise remainder",
           "(\\log X)^4" not in op7 and
           ("e^u/u^3" in op7 or "\\frac{e^u}{u^3}" in op7))
    record("T", "Methods carries the differentiated reproducibility form",
           "drift anchors" in body_nc and
           "recomputed from" in body_nc and
           "certificate artefacts" in body_nc)

    # ---- cycle-8 watchers (seventh review) ------------------------------
    record("T", "the causal identification phrase is gone",
           "collapse of multiplicativity, not of counting density"
           not in body_nc and
           "visible in both directions" not in body_nc)
    m = re.search(r"Numerical observation\} \(function-side(.*?)\\subsection",
                  body, re.DOTALL)
    mir = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "the mirror paragraph carries the non-identification clause",
           "multiplicative" in mir and "point geometry together" in mir and
           "no common causal feature is established" in mir)
    record("T", "the mirror paragraph declares itself a separate experiment",
           "separate" in mir and "not part of the" in mir)
    record("T", "the record figures are quoted from the artefact",
           "98.1" in mir and "34.6" in mir and "85.9" in mir and
           "n=6=2\\cdot3" in mir)
    record("T", "function_side_reconstruction.py is named in the methods section",
           "function\\_side\\_reconstruction.py" in body)
    m = re.search(r"\\textbf\{What is numerical\.\}(.*?)\\textbf",
                  body, re.DOTALL)
    numl = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "the abstract numerical list carries the experiment",
           "function-side reconstruction" in numl)
    m = re.search(r"establishes numerically(.*?)labelled", body, re.DOTALL)
    ctr = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "the reader contract carries the experiment",
           "function-side reconstruction" in ctr)
    record("T", "the Liu record is pinned to v2",
           "17645v2" in body and
           "arXiv:2605.17645," not in body.replace("\\allowbreak", ""))

    # ---- cycle-9 watchers (eighth review) -------------------------------
    record("T", "the persistence overclaim is gone",
           "scaling laws persist" not in body_nc)
    record("T", "Theorem 4.4 carries the transfer title form",
           "Sideband transfer of the leading representation" in body_nc)
    record("T", "the delimitation sentence is present",
           "No general assertion" in body_nc and
           "stands as stated" in body_nc)
    record("T", "sideband law is retired in favour of transfer",
           "sideband law" not in body_nc.lower() and
           body_nc.count("sideband transfer") >= 5)
    record("T", "all three metadata places carry August 2026",
           "% Date:     August 2026" in body and
           body.count("August~2026") == 2 and
           "July~2026" not in body)

    # ---- cycle-10 watchers (ninth review) -------------------------------
    body_ws = re.sub(r"\s+", " ", body_nc)
    record("T", "Theorem 3.1 states the reference-data hypothesis",
           "$(N,\\varepsilon)=(100,0.05)$" in body and
           "reference ordinate and regularisation data" in body_ws)
    record("T", "the scope-free held-fixed clause is gone",
           "the ordinate set and the regularisation width being held fixed"
           not in body_ws)
    record("T", "Theorem 4.2 carries the reference-data anchor sentence",
           "Keep the same reference ordinate" in body_ws and
           "vary only the prime-side Beurling system" in body_ws)
    record("T", "the synopsis carries the reference-data half-sentence",
           body_ws.count("reference data $(100,0.05)$") >= 3)
    record("T", "the notation table prints the series form of B_int",
           body.count("B_{\\mathrm{int}}") ==
           body.count("B_{\\mathrm{int}}^{+}") and
           "B_{\\mathrm{int}}^{+}(\\varepsilon,N)" in body)
    record("T", "zeta_A is retired in favour of zeta_sm",
           re.search(r"\\zeta_\{?A\}?", body) is None and
           body.count("\\zeta_{\\mathrm{sm}}") >= 2)

    # ---- cycle-11 watchers (tenth review, release) ----------------------
    m = re.search(r"Step 3 \(uniformity\)(.*?)\\end\{proof\}", body, re.DOTALL)
    st3 = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "the invalid construction is deleted from Step 3 (RE-KEYED)",
           st3 != "MISSING" and "adjoin" not in st3 and
           "accumulation just above" not in st3)
    record("T", "finite deletion stands as a separate statement (RE-KEYED)",
           "separate statement about the union" in st3 and
           "finite deletion" in st3)
    record("T", "rational independence is marked as unproved",
           ("rational independence" not in body_ws) and
           ("unproved hypothesis" in body_ws) and
           ("linearly independent over $\\mathbb{Q}$" in body_ws))
    m = re.search(r"No inference from a discriminator to a location\.\}(.*?)"
                  r"Theorem~\\ref\{thm:barrier\}", body, re.DOTALL)
    ncirc = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "the non-circularity paragraph carries the two-level form",
           "two levels must be kept apart" in ncirc and
           "never transferred" in ncirc and
           "one of them has no zeros at all" not in ncirc)
    record("T", "the 0.065 threshold carries its scope",
           "24$-mode reference set" in body_ws and
           "one-sided empirical rule" in body_ws and
           "not monotone" in body_ws)

    # ---- cycle-12 watchers (eleventh review, final fix) -----------------
    record("T", "the Tt-is-a-kernel prose form is gone",
           "is} a kernel, evaluated at arithmetic points" not in body_ws)
    record("T", "the identity is attributed to T, not Tt",
           "entrywise formula for $T$, not for $\\Tt$" in body_ws)
    record("T", "abstract and 1.2 carry the companion-Gram wording",
           body_ws.count("companion prime-space Gram operator") >= 2)
    record("T", "the channel sentence is present",
           "first channel" in body_ws and
           "logarithmic differences" in body_ws and
           "logarithmic sums" in body_ws)
    record("T", "the alpha claim carries the Stieltjes derivation",
           "Stieltjes partial summation against $E$" in body_ws and
           "\\frac{\\sqrt{x}}{(\\log x)^{\\alpha}}" in body)
    record("T", "the narrow alpha claim is gone, every alpha>0 stated",
           "with $0<\\alpha\\le 2$ preserves the laws" not in body_ws and
           "for \\emph{every} $\\alpha>0$" in body_ws)

    # ---- cycle-13 watchers (twelfth review, release fixes) --------------
    record("T", "the abstract carries the unless-clause",
           "first surviving term is quadratic. This is a theorem"
           not in body_ws and
           body_ws.count("unless the aggregate quadratic coefficient") >= 2)
    record("T", "the abstract carries the weak two-margin form",
           "separates from three" not in body_ws and
           "Certified separation margins are established against the smooth "
           "measure and the discrete world" in body_ws)
    record("T", "the localisation is labelled with its parameters",
           "with an overlap of $0.94$ against" not in body_ws and
           "Numerical observation} (leading-mode localisation" in body_ws and
           "\\kappa\\approx10^{4}" in body and
           "squared overlap" in body_ws)
    record("T", "the 312 zeros belong to the extended list",
           "zeros of the normative" not in body_ws and
           body.count("zeros\\_650.csv") >= 3)
    record("T", "the methods sentence is the normative-qualified form",
           "No further parallel ordinate file exists" not in body_ws and
           "No further \\emph{normative} ordinate file exists" in body_ws and
           "defines no anchor and no tolerance" in body_ws)

    # ---- cycle-14 watchers (thirteenth review) --------------------------
    record("T", "the M0 openness is recorded in OP 9.2",
           "can be dispensed with at fixed $(C,x_0)$" in body_ws and
           "not decided here" in body_ws)

    # ---- cycle-15 watchers (fourteenth review) --------------------------
    m = re.search(r"\\label\{thm:barrier\}(.*?)\\end\{proof\}", body, re.DOTALL)
    bar = re.sub(r"\s+", " ", m.group(1)) if m else "MISSING"
    record("T", "G does not carry the orbit average in Thm 6.1",
           bar != "MISSING" and "G(\\gamma_k+i\\delta" not in bar and
           "G(\\gamma_k-i\\delta" not in bar)
    record("T", "the K kernel is defined with the DISTINCT sentence",
           "\\mathcal K_{\\kappa,\\varepsilon}(z)" in body and
           "must not be interchanged" in body_ws and
           "sum over \\emph{primes} in the \\emph{ordinate} variable"
           in body_ws)
    record("T", "the anchor identity sentence is present",
           "\\mathcal K_{\\kappa,\\varepsilon}(\\gamma_k)=B=-19342.5476"
           in body)
    record("T", "the damped-inverse-spectrum orbit phrase is retired",
           "orbit averages of the damped inverse spectrum" not in body_ws and
           "entire weighted curvature test kernel" in body_ws)
    # [R] attribution watcher: recompute sum_k K(gamma_k) == B
    _P15 = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53]
    _lp15 = np.log(np.array(_P15, dtype=float))
    _S15 = float(np.sum((gam**2*np.exp(-EPS**2*gam**2))[:,None]
                        *(_lp15**2)[None,:]*np.cos(np.outer(gam,_lp15))))
    record("R", "sum_k K(gamma_k) recomputes the series anchor B",
           abs(_S15-(-19342.5476)) < 1e-3, f"{_S15:.4f}")

    # ---- cycle-16 watchers (fifteenth review, record E-6) ---------------
    record("T", "the nearest-rounded margin literal is gone",
           "5343.91" not in body)
    record("T", "the smooth margin is the downward bound at all sites",
           body.count("5343.90") >= 3 and
           "downward-rounded\ncertified lower bounds" in body or
           "downward-rounded certified lower bounds" in body_ws)
    record("T", "the even-order close carries its qualifiers",
           "unconditionally" in body_ws and
           "differentiable in $\\delta$" in body_ws and
           "only even powers" in body_ws and
           "begins at even order.\\end{theorem}" not in body_ws.replace(" ",""))
    record("T", "the classification row is the qualified form",
           "first derivative zero when" in body_ws and
           "no first-order response" not in body_ws)
    record("T", "the percent claim is the point statement",
           "1.247" in body and "$1.2\\%$ and $1.7\\%$" not in body)
    # ---- cycle-17 watchers (sixteenth review) ---------------------------
    record("T", "the bare-W display form is gone",
           "\\frac{4}{W}" not in body.replace("\\frac{4}{\\WeN}",""))
    record("T", "the naked W definition is gone",
           "with $W=\\sum_k w_k^2" not in body_ws)
    record("T", "the W_{eps,N} carrier is defined and used in Thm 3.4",
           "\\WeN:=\\sum_{k=1}^{N}e^{-\\varepsilon^2\\gamma_k^2}" in body
           and "\\frac{4}{\\WeN}" in body
           and "\\newcommand{\\WeN}{W_{\\varepsilon,N}}" in body)
    m = re.search(r"\\begin\{openproblem\}(.*?)\\end\{openproblem\}",
                  body, re.DOTALL)
    op1 = m.group(1) if m else "MISSING"
    m=re.search(r"Background and Motivation(.*?)Main Results", body, re.DOTALL)
    p1 = m.group(1) if m else "MISSING"
    record("T", "Paper-4 status is the record wording (section 1)",
           "established that it is not a Hilbert" not in re.sub(r"\s+"," ",p1)
           and "no numerical support" in p1 and "finite-cutoff grid" in p1)
    record("T", "the categorical HP form is gone from Non-Circularity",
           "is not a Hilbert--P\'olya operator\n(\\cite{Paper4})" not in body
           and "numerical\ndiagnostic failure" in body or
           "numerical diagnostic failure" in body_ws)
    record("T", "K_eps is defined before Theorem 5.3",
           "K_\\varepsilon(t):=(\\log t)^2" in body)
    record("T", "the Paper-6 citation targets section 6",
           "(\\cite{Paper6}, \\S\\,6)" in body and
           "(\\cite{Paper6}, \\S\\,7)" not in body)
    # ---- cycle-29 watchers (28th review) --------------------------------
    record("T", "the completed-mode Gram identity is stated as algebraic",
           body_nc.count("never the primality of the label") == 1 and
           "Theorem~\\ref{thm:kernel} holds verbatim" in body_nc and
           "[Restriction to prime modes as a design choice; numerical]" in body_nc)

    # ---- cycle-28 watchers (27th review) --------------------------------
    record("T", "the Thm-5.2 candidate count is 31, not sixteen",
           "sixteen one-sided limits" not in body_nc and
           body_nc.count("$31$ one-sided values") == 1)
    record("T", "Thm 5.2 carries the two-part status of Thm 3.1",
           body_nc.count("proved, finite ordering certified") == 1 and
           "proved, separation certified" in body_nc and
           "\\ge0.1201" in body_nc and "\\ge0.0920" in body_nc and
           "the runner-up being" not in body_nc)
    import os as _os3
    _csc = ""
    for _c in (_os3.path.join(_os3.path.dirname(_os3.path.abspath(__file__)),
                              "certify_smooth_controls.py"),
               "certify_smooth_controls.py"):
        if _os3.path.exists(_c):
            _csc = open(_c, encoding="utf-8").read(); break
    record("T", "the smooth-control certificate certifies the full ranking",
           'ARGMAX_KEY = "53^-"' in _csc and 'RUNNER_KEY = "41^-"' in _csc and
           "ranked = sorted(cands" in _csc and "directed separation" in _csc and
           "certified ranking disagrees" in _csc and "directed links" in _csc and
           "max(vals)" not in _csc and "max(others, key=" not in _csc)

    # ---- cycle-27 watchers (26th review) --------------------------------
    _cmp = "\\bigl(1-(1+4\\gamma_1^2)^{-1/2}\\bigr)"
    _dsp = "\\Bigl(1-\\frac{1}{\\sqrt{1+4\\gamma_1^2}}\\Bigr)"
    record("T", "the profile positivity rests on the exact envelope",
           body_nc.count(_cmp) == 2 and body_nc.count(_dsp) == 1 and
           body_nc.find(_dsp) < body_nc.find(_cmp) < body_nc.rfind(_cmp))
    record("T", "the superseded rounding forms are gone",
           "0.0020938742" not in body_nc and "0.861153" not in body_nc
           and "0.978738" not in body_nc)
    record("T", "bound decimals are written with an inequality sign",
           "\\le 0.861154" in body_nc and "\\le 0.978739" in body_nc and
           "\\ldots\\approx3.184690" in body_nc)
    record("T", "eta is named as the energy-asymmetry functional",
           body_nc.count("energy-asymmetry") >= 5 and
           len(re.findall(r"Rayleigh", body_nc)) == 3)

    # ---- cycle-26 watchers (25th review) --------------------------------
    record("T", "the Paper-7 status sentence claims no RH implication",
           "strictly weaker" not in body_nc and
           body_nc.count("not known to follow from the Riemann Hypothesis") == 1)
    import os as _os2
    def _artefact(_name):
        for _c in (_os2.path.join(_os2.path.dirname(_os2.path.abspath(__file__)),
                                  _name), _name):
            if _os2.path.exists(_c):
                return open(_c, encoding="utf-8").read()
        return ""
    _cr = _artefact("cert_paper9_ratio.py")
    _fs = _artefact("function_side_reconstruction.py")
    record("T", "the two-point certificate enforces its N = 100 contract",
           "EXPECTED_N = 100" in _cr and
           "must contain exactly {EXPECTED_N}" in _cr and
           'default="zeros_100.csv"' not in _cr)
    record("T", "the function-side experiment enforces its 650-row contract",
           "EXPECTED_ROWS = 650" in _fs and
           "must contain exactly {EXPECTED_ROWS} rows" in _fs)

    # ---- cycle-25 watchers (24th review) --------------------------------
    record("T", "no Beurling index remnants anywhere in the document",
           re.search(r"\\vartheta_B\b", body) is None and
           re.search(r"\\pi_B\b", body) is None and
           re.search(r"\\Dclass\s*=\s*\\\{\\vartheta", body) is None)
    record("T", "the global notation block binds the class over frak-B",
           "\\Dclass=\\{\\mathfrak{B}:\\ \\vartheta_{\\mathfrak{B}}(x)=x+O(" in body)
    import os as _os
    _cert = ""
    for _c in (_os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                            "certify_smooth_controls.py"),
               "certify_smooth_controls.py"):
        if _os.path.exists(_c):
            _cert = open(_c, encoding="utf-8").read(); break
    record("T", "the d_BL certificate is fail-closed at both guards",
           "d_BL primal witness is not certified feasible" in _cert and
           "d_BL enclosure is inconsistent" in _cert)

    # ---- cycle-24 watchers (23rd review) --------------------------------
    _i0=body.index("\\label{thm:sideband}")
    _pr=body.index("\\begin{proof}", _i0)
    _blk=body[_pr:body.index("\\end{proof}", _pr)]
    _ic=_blk.find("\\textit{Count-side transfer.}")
    _isp=_blk.find("\\textit{Spectral-energy transfer.}")
    _inc=_blk.find("Theorem~\\ref{thm:nonconv}, Step~1")
    _ie=_blk.find("Theorem~\\ref{thm:eta}, Step~2")
    record("T", "the section-8 battery names its own zeta regeneration",
           "for the $\\zeta$-layer" not in body and
           "\\texttt{mpmath.zetazero}" in body and
           re.search(r"first \$31\$\s*\n?ordinates", body) is not None and
           body.count("not an input to either Arb certificate") == 1)
    record("T", "the Thm-4.4 cross-references are correct AND ordered",
           _blk.count("Theorem~\\ref{thm:nonconv}, Step~1") == 1 and
           _blk.count("Theorem~\\ref{thm:eta}, Step~2") == 1 and
           -1 not in (_ic, _isp, _inc, _ie) and _ic < _inc < _isp < _ie)

    # ---- cycle-23 watchers (22nd review) --------------------------------
    record("T", "the Thm-4.4 proof carries the count-measure display",
           "d\\Pi_{a,\\omega}(t):=\\frac{d\\vartheta_{a,\\omega}(t)}{\\log t}" in body
           and "Count-side transfer" in body)
    record("T", "the Thm-4.4 proof carries the spectral-side display",
           "t^{-1/2+i\\gamma}\\,d\\vartheta_{a,\\omega}(t)" in body
           and "Spectral-energy transfer" in body
           and re.search(r"giving the printed \$a\^\{2\}/8\$ shift", body) is not None)
    record("T", "the certificate provenance names acb, not gen_zeros",
           re.search(r"performed by\s+\\texttt\{gen\\_zeros", body) is None and
           "not an input to either Arb certificate" in body)
    record("T", "Thm 5.5 carries the delimitation bracket and sentence",
           re.search(r"Factorisation barrier: bounded-Lipschitz resolution lower\s+bound; proved from certified inputs", body) is not None
           and "no stronger obstruction to factorisation is claimed" in body)

    # ---- cycle-22 watchers (21st review; minimum only, see note) --------
    _s0=body.index("\\section{Class Invariance")
    _s1=body.index("\\section{The Certified Factorisation")
    _sec4=body[_s0:_s1]
    _r0=body.index("\\begin{remark}[Finite insertion")
    _rem=body[_r0:body.index("\\end{remark}",_r0)]
    record("T", "section-4 segment is free of bare A-index forms",
           not any(re.search(p,_sec4) for p in (r"\br_A\b", r"\bO_A\b",
           r"E\^\{A\}", r"E\^A\b", r"\\pi_A", r"\\vartheta_A", r"\\Phi_A",
           r"\\eta_A", r"\\\{A\\in", r"Let \$A\\in", r"index \$A\$")))
    record("T", "the frak-B set binder defines the M0 class",
           "\\bigl\\{\\mathfrak{B}\\in\\Dclass(C,x_0):" in _sec4)
    record("T", "the cannot-distinguish and fingerprint overclaims are gone",
           "cannot distinguish" not in body and
           "fingerprint of a density class" not in body and
           "no fingerprint" in body)
    record("T", "both cycle-22 replacement texts are anchored",
           "asks what the stated" in body and "within-class comparison" in body
           and re.search(r"does not\s+establish uniqueness", body) is not None)
    record("T", "the Remark hypothesis is the uniform-slack form",
           re.search(r"fixed\s+positive uniform slack", body) is not None and
           "strict" not in _rem and
           re.search(r"C'\\,b\(x\)\+\\varepsilon_0\\le C\\,b\(x\)", body) is not None)

    # ---- cycle-21 watchers (twentieth review) ---------------------------
    _s0=body.index("\\section{Class Invariance")
    _s1=body.index("\\section{The Certified Factorisation")
    _sec4=body[_s0:_s1]
    record("T", "a_q^eta is defined before Theorem 4.2",
           "a_q^{\\eta}:=\\bigl(e^{-\\varepsilon^2\\gamma_k^2/2}" in body and
           body.index("a_q^{\\eta}:=") < body.index("\\label{thm:class}"))
    record("T", "the exhausted/all-they-retain overclaims are gone",
           re.search(r"exhausted by the\s+density class", body) is None and
           "all they retain" not in body)
    record("T", "the leading-law qualifier forms are present",
           "lower-order discrimination not excluded" in body and
           re.search(r"contributions are not\s+excluded", body) is not None)
    record("T", "OP 9.2 carries the sharpened uniformity clause",
           re.search(r"dispensed with at fixed \$\(C,x_0\)\$ is not", body) is None
           and "$(C,x_0)$: whether" in body and "bypass uniform" in body
           and "Finite insertion near $1$" in body)
    record("T", "the section-4 system index is frak-B (A-forms gone there)",
           "\\mathfrak{B}" in _sec4 and not any(re.search(p,_sec4) for p in
           (r"\\vartheta_A", r"\\Phi_A", r"r_A\(", r"\\eta_A",
            r"E\^A_", r"\\Tt_A", r"\\pi_A", r"\\mathcal\{P\}_A")))
    record("T", "artefact names are the public ones",
           all(x not in body for x in ("smooth\\_a.py","invariant\\_2\\_g2.py",
           "area\\_bridge.py")) and
           all(x in body for x in ("certify\\_smooth\\_controls.py",
           "unfolded\\_discrimination.py","function\\_side\\_reconstruction.py")))
    record("T", "the initiative phrase replaced the RH phrase",
           re.search(r"toward the\s+Riemann Hypothesis", body) is None and
           "open research initiative on prime--zero" in body)

    _i1=body_ws.find("Step 1 (positive semidefiniteness")
    _i2=body_ws.find("Step 2 (ratio perturbation")
    _i3=body_ws.find("Step 3 (Riemann--von-Mangoldt")
    _i4=body_ws.find("Step 4 (transfer to the companion")
    record("T", "the four Thm-3.5 step markers are present AND ordered",
           -1 not in (_i1,_i2,_i3,_i4) and _i1<_i2<_i3<_i4)
    record("T", "the post-theorem reprise carries the unless-qualifier",
           re.search(r"quadratic\s+term,\s+unless the aggregate quadratic coefficient", body) is not None
           and re.search(r"quadratic term\.\s+What the theorem bars", body) is None)
    record("T", "the Thm-3.5 rate sentence is the precise form",
           "for fixed $\\kappa$ and $\\varepsilon>0$" in body_ws and
           "Gaussian ordinate tail $e^{-\\varepsilon^2\\gamma_N^2}$" in body_ws
           and "Riemann--von-Mangoldt factor" in body_ws)

    record("T", "OP 9.1 carries the WeN carrier form",
           "\\frac{4e^uf(e^u)}{\\WeN u^2}" in op1)
    record("T", "the bare-W OP form is gone",
           "\\frac{4e^uf(e^u)}{Wu^2}" not in op1 and
           "{W u^2}" not in op1)

    record("T", "B_sm matches the v2.8.56 contract symbol",
           "\\newcommand{\\Bsm}{B_{\\mathrm{sm}}}" in body)

    record("T", "section 8 names both acb certificate paths and v2.0",
           "certify\\_smooth\\_controls.py} (v2.0)" in body_ws and
           "acb.zeta\\_zeros" in body and
           "ess\\_gate.py" in body)
    record("T", "the cosh attribution is gone",
           "cosh(2\\delta" not in body)
    record("T", "the even-orbit-average attribution is present (RE-KEYED)",
           body.count("\\mathcal K_{\\kappa,\\varepsilon}(\\gamma_k+i\\delta)") >= 2
           and "even orbit average" in body_ws)
    record("T", "the two-line expansion is in the Thm-6.1 proof (RE-KEYED)",
           "-\\tfrac{\\delta^2}{2}\\,\\mathcal K_{\\kappa,\\varepsilon}''" in body
           and "O(\\delta^4)" in body)
    record("T", "the function-side paragraph is active and labelled",
           "Numerical observation} (function-side reconstruction" in body_nc
           and "directionally compatible" in body_nc)


if __name__ == "__main__":
    sys.exit(main())
