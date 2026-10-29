"""recompute selected T8 edges (TL=100, TR=60) and merge into T8_cert json.  usage: t8_patch.py N J G tol edge_lo edge_hi"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, math
import cert_T8 as C
from multiprocessing import Pool
from flint import arb
N, J, G = sys.argv[1:4]; tol = float(sys.argv[4]); lo, hi = float(sys.argv[5]), float(sys.argv[6])
f = _wout('g0_constants', 'T8_cert_%s_%s_%s.json' % (N, J, G)); d = json.load(open(f))
Na, Ja, Ga = arb(N), arb(J), arb(G); g = Ga/Ja; nu1 = 1 + Na; ds = 4*Na/Ja*(2*Na/(1 + Na - g) - 1)
nu1s, gs, dss = nu1.mid().str(20, radius=False), g.mid().str(20, radius=False), ds.mid().str(20, radius=False)
ED = d['edges']; idx = [i for i in range(len(ED) - 1) if ED[i] >= lo - 1e-12 and ED[i + 1] <= hi + 1e-12]
jobs = [(ED[i], ED[i + 1], nu1s, gs, dss, str(d['sigma_min']), 100.0, 60.0, tol) for i in idx]
with Pool(2) as P: res = P.map(C.cell_sup, jobs[::-1], chunksize=1)
res = res[::-1]
te = {(a, b): th for a, b, th in d['Theta_edge']}
for i, r in zip(idx, res):
    d['cells'].append(dict(r, edge=i, patch=True))
    new = r['Theta'] if r['ok'] else math.inf
    old = te.get((ED[i], ED[i + 1]), math.inf)
    te[(ED[i], ED[i + 1])] = min(old, new)       # both are valid upper bounds
    print(ED[i], ED[i + 1], r.get('Theta'), r.get('why'), 'old', old, flush=True)
d['Theta_edge'] = [[a, b, te[(a, b)]] for (a, b) in sorted(te)]
json.dump(d, open(f, 'w'), indent=1)
