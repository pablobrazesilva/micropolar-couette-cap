"""Fine certified table of P(S') upper bounds (eta=0.05), using pbound.C_bound(0, S') (rigorous, see G0-constantes 4.4).
P is nonincreasing, so P(S') <= P(S'_k) for S' >= S'_k.  Writes ptable_0.05.json (values rounded UP)."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, math
sys.argv = ['pbound.py', '0.05']
sys.path.insert(0, _os.path.join(_R, 'g0_constants'))
import pbound as PB
grid = [round(1 + 0.25*k, 4) for k in range(0, 101)]         # 1 .. 26
grid = sorted(set(grid + [round(8 + 0.05*k, 4) for k in range(0, 61)]))   # finer on [8, 11]
tab = []
for Sp in grid:
    v = PB.C_bound(0.0, Sp)
    e = math.floor(math.log10(v)); m = math.ceil(v/10**e*1e4)/1e4       # round up, 5 significant digits
    tab.append((Sp, '%.4fe%d' % (m, e)))
json.dump(tab, open(_wout('g0_constants', 'ptable_0.05.json'), 'w'))
for Sp, v in tab[::8]: print(Sp, v)
print(len(tab), 'entries')
