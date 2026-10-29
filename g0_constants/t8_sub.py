"""T8 for a lower sigma_min, only beta cells inside [blo, bhi] -> T8_cert_<N>_<J>_<G>_s<smin>.json (same format)."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, math, time
import cert_T8 as C
from multiprocessing import Pool
from flint import arb
N, J, G, smin = sys.argv[1:5]; blo, bhi = float(sys.argv[5]), float(sys.argv[6])
Na, Ja, Ga = arb(N), arb(J), arb(G); g = Ga/Ja; nu1 = 1 + Na; ds = 4*Na/Ja*(2*Na/(1 + Na - g) - 1)
nu1s, gs, dss = nu1.mid().str(20, radius=False), g.mid().str(20, radius=False), ds.mid().str(20, radius=False)
ED = [0, 0.002, 0.005, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.8, 1, 1.3, 1.7, 2.2, 3, 4, 5, 6.5]
cells = [(ED[i], ED[i + 1], i) for i in range(len(ED) - 1) if ED[i] >= blo - 1e-12 and ED[i + 1] <= bhi + 1e-12]
TOLX = float(sys.argv[7]) if len(sys.argv) > 7 else None
tol = lambda b: TOLX if TOLX else (0.04 if b <= 1.3 else 0.1)
jobs = [(a, b, nu1s, gs, dss, smin, 100.0, 60.0, tol(b)) for (a, b, i) in cells][::-1]
T0 = time.time()
with Pool(2) as P: res = P.map(C.cell_sup, jobs, chunksize=1)
res = res[::-1]
out = dict(NJG=[N, J, G], sigma_min=float(smin), edges=ED, cells=[dict(r, edge=i) for (a, b, i), r in zip(cells, res)],
           Theta_edge=[[a, b, (r['Theta'] if r['ok'] else math.inf)] for (a, b, i), r in zip(cells, res)], seconds=time.time() - T0)
fn = _wout('g0_constants', 'T8_cert_%s_%s_%s_s%s.json' % (N, J, G, smin))
import os
if os.path.exists(fn):        # merge: keep the smaller valid bound per cell
    old = json.load(open(fn)); te = {(a, b): th for a, b, th in old['Theta_edge']}
    for a, b, th in out['Theta_edge']: te[(a, b)] = min(te.get((a, b), math.inf), th)
    out['cells'] = old['cells'] + out['cells']; out['Theta_edge'] = [[a, b, te[(a, b)]] for (a, b) in sorted(te)]
json.dump(out, open(fn, 'w'), indent=1)
for r in out['Theta_edge']: print('beta in [%g,%g]: Theta <= %.4f' % tuple(r))
print('seconds', out['seconds'])
