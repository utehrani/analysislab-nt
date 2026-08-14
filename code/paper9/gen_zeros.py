import mpmath as mp, csv
import os as _os
def _z100():
    _sd = _os.path.dirname(_os.path.abspath(__file__))
    for _p in (_os.path.join(_sd, '..', '..', 'data', 'zeros_100.csv'),
               _os.path.join(_sd, 'data', 'zeros_100.csv'),
               _os.path.join('data', 'zeros_100.csv'),
               _os.path.join(_sd, 'zeros_100.csv'),
               'zeros_100.csv'):
        if _os.path.exists(_p): return _p
    raise SystemExit('zeros_100.csv not found (repo data/, script dir, CWD)')
mp.mp.dps = 30
rows=[]
for k in range(1,101):
    g = mp.zetazero(k).imag
    rows.append((k, mp.nstr(g, 28)))
_sd = _os.path.dirname(_os.path.abspath(__file__))
for _r in (_os.path.join(_sd, '..', '..', 'data', 'results'),
           _os.path.join(_sd, 'data', 'results'),
           _os.path.join('data', 'results')):
    if _os.path.isdir(_os.path.dirname(_r)):
        _out = _r; break
else:
    _out = _sd
_os.makedirs(_out, exist_ok=True)
with open(_os.path.join(_out, 'zeros_100_dps30.csv'), 'w') as f:
    f.write('k,gamma\n')
    for k,g in rows: f.write(f'{k},{g}\n')
# cross-check vs normative CSV
# v-bump (revision, Aug 2026): path cascade for zeros_100.csv —
#   script dir, then CWD; no fixed mount path (review C6).
# v-fix (later revision): the cascade FUNCTION itself was missing here —
#   the earlier bump copied only the call (regression caught in review).
import math
ref = {}
with open(_z100()) as f:
    next(f)
    for line in f:
        k,g = line.split(','); ref[int(k)] = float(g)
if len(ref) < 100:
    raise SystemExit(f'zeros_100.csv must contain at least 100 rows for the '
                     f'cross-check; got {len(ref)}')
maxd = max(abs(float(mp.mpf(g)) - ref[k]) for k,g in rows)
print('zeros computed, max |zetazero - csv| =', maxd)
