"""Tighten the certified lower bounds of a(eta) (key 'u') or at(eta) (key 'ut'): every box of line_<eta>.json whose
   lower bound is <= target is bisected (prec 96) until its bound exceeds target.  usage: refine.py eta key target"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json; sys.path.insert(0, _os.path.join(_R, 'g0_constants'))
from cl import *
ctx.prec = 96
eta, key, target = sys.argv[1], sys.argv[2], float(sys.argv[3])
d = json.load(open(_os.path.join(_CERT, 'g0_constants', 'line_%s.json' % eta))); col = 4 if key == 'u' else 5
tails = d['tail_right'][col-4], d['tail_left'][col-4]
assert min(tails) > target
bad = [(b[0], b[1]) for b in d['boxes'] if b[col] <= target]
n = 0; worst = 9; minw = 9
st = list(bad)
while st:
    a1, a2 = st.pop(); n += 1
    r = box_line(a1, a2, eta)
    v = None if r is None else r['ulo' if key == 'u' else 'utlo']
    if v is not None and v > target:
        worst = min(worst, v); minw = min(minw, a2 - a1); continue
    if a2 - a1 < 1e-7: print('FAIL', a1, a2, v); sys.exit(1)
    m = 0.5*(a1 + a2); st += [(a1, m), (m, a2)]
print('eta=%s %s: CERTIFIED  %s > %.9f on the whole line  (%d coarse boxes refined, %d evaluations, min width %.2e)' % (eta, key, 'a' if key == 'u' else 'at', target, len(bad), n, minw))
