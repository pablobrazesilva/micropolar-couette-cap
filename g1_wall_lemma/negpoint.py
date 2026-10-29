"""Negative control: at (sigma, beta in [0,1e-9]) all 16 parameter choices of the criterion must give margin <= 0."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys; sys.path.insert(0, _os.path.join(_R, 'g1_wall_lemma')); sys.path.insert(0, _os.path.join(_R, 'g1_wall_lemma/heuristic'))
from arbcrit import crit_box, params
from arbprim import arb, iv
N, J, G, s = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]); cp = int(sys.argv[5]) if len(sys.argv) > 5 else 1
P = params(N, J, G); best = -9
for lam in (0, 1):
    for th in (0.35, 0.5, 0.65, 0.8):
        for cone in (True, False):
            o = crit_box(P, arb(s), iv(0, 1e-9), lam, th, cone, coupling=cp)
            best = max(best, o['margin'])
print('(%s,%s,%s) sigma=%g coupling=%d, beta in [0,1e-9]: best rigorous margin over 16 choices = %.4f  -> %s'
      % (N, J, G, s, cp, best, 'NOT certified (as expected)' if best <= 0 else 'certified'))
