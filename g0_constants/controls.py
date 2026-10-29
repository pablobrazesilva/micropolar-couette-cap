import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
# controls for the G1 constants: (i) negative control for the cone at eta=0.05/0.1, (ii) rigorous LOWER bounds
# for P(S') (point values) to show the certified upper bounds are sharp.
import sys; sys.path.insert(0, _os.path.join(_R, 'g0_constants')); sys.path.insert(0, _os.path.join(_R, 'g0_sector_lemma/audit'))
from cl import *
from negctl import attempt
print('cone eta=0.05, th2=0.54935 (< max 0.5493617):', attempt(-3.58, -3.48, 0, 0.05, '0.05', '-0.0836', '0.54935', None))
print('cone eta=0.10, th2=0.56135 (< max 0.5613603):', attempt(-3.56, -3.46, 0, 0.05, '0.1', '-0.0837', '0.56135', None))
import numpy as np
from scipy.optimize import minimize_scalar
ctx.prec = 96
def ratio(t, Sp, eta):
    t = float(t)
    F, u, ut = point(repr(t), eta, 96)
    z0 = acb(arb(repr(t)))*EM - arb(eta)*E3; z1 = acb(arb(repr(t)) + arb(repr(Sp)))*EM - arb(eta)*E3
    return abs(F)*abs(z1.airy_ai()/z0.airy_ai())
for eta in ['0.05']:
    for Sp in [2, 6, 10, 14, 20]:
        ts = np.linspace(-12, 4, 161); v = [float(ratio(t, Sp, eta).mid()) for t in ts]; t0 = ts[int(np.argmax(v))]
        r = minimize_scalar(lambda t: -float(ratio(t, Sp, eta).mid()), bounds=(t0-0.1, t0+0.1), method='bounded')
        print("eta=%s S'=%d: P(S') >= %s  (t=%.4f)" % (eta, Sp, ratio(r.x, Sp, eta).str(6), r.x))
