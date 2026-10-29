"""Rigorous certification of the crit3 two-point tables OmA(l), OmPsi(l) (eta = 0.05) on the grid l = 0:0.025:16:
   OmA(l)   = sup_{Z in H}|Ai(Z+l e)/Ai(Z)|   = sup_t exp(-int_t^{t+l} ut),   ut = Re(-e^{-i pi/6} Ai'/Ai) >= at > 0
   OmPsi(l) = sup_{Z in H}|psi(Z+l e)/psi(Z)| = sup_t exp(-int_t^{t+l} u),    u  = Re( e^{-i pi/6} F)     >= a  > 0
   (e = e^{-i pi/6}).  Both functions are holomorphic and bounded (<= 1) in H_eta, so the sup is on the line (Phragmen-Lindelof);
   on the line z(t)+l e = z(t+l).  Uses the certified per-box lower bounds of u, ut (line_0.05.json) + explicit tails (Lemmas B, D).
   Writes tp_cert_0.05.json and compares with the numerical tables of g1/coupling/tp.py."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, math, bisect
sys.path.insert(0, _os.path.join(_R, 'g0_constants'))
from cl import *
import numpy as np
eta = '0.05'
d = json.load(open(_os.path.join(_CERT, 'g0_constants', 'line_%s.json' % eta))); boxes = d['boxes']; TL, TR = d['TL'], d['TR']
etaA = arb(eta)
def tail_rate(Tabs, side, which):
    """rigorous lower bound of the rate (u: which=0, ut: which=1) for |t| >= Tabs on the given side"""
    return tail_side(Tabs, eta, side)[which]
def table(which, Lgrid, TF=40.0, hc=1/16):
    col = 4 if which == 0 else 5
    cells = []
    n = int(round((TF - TL)/hc))
    for i in range(n):                                        # left tail cells [-TF, -TL]
        a1 = -TF + i*hc; a2 = a1 + hc
        cells.append((a1, a2, tail_rate(-a2, 'left', which)))
    cells += [(b[0], b[1], b[col]) for b in boxes]
    rR = tail_rate(TR, 'right', which); rL = tail_rate(TF, 'left', which)
    lo_rate = min(min(c[2] for c in cells), rR)
    assert lo_rate > 0
    # exact minimisation of the piecewise-linear window integral f(t) = Fc(t+l) - Fc(t), Fc = cumulative lower rate:
    # minimum over t is attained at breakpoints t in {edges} U {edges - l}; outside [-TF, TR] slopes rL, rR.
    edges = np.array([cells[0][0]] + [c[1] for c in cells])
    rt = np.array([c[2] for c in cells])
    Fe = np.concatenate([[0.0], np.cumsum(rt*np.diff(edges))])
    def Fc(x):
        x = np.asarray(x, float)
        v = np.interp(x, edges, Fe)
        v = np.where(x < edges[0], Fe[0] - rL*(edges[0] - x), v)
        v = np.where(x > edges[-1], Fe[-1] + rR*(x - edges[-1]), v)
        return v
    out = []
    for l in Lgrid:
        if l == 0: out.append(1.0); continue
        cand = np.concatenate([edges, edges - l])
        f = Fc(cand + l) - Fc(cand)
        fmin = min(f.min(), l*min(rL, rR))
        out.append(math.exp(-fmin*(1 - 1e-9))*(1 + 1e-9))
    return np.array(out)
if __name__ == '__main__':
    L = np.arange(0, 16.0001, 0.025)
    OmPsi = table(0, L); OmA = table(1, L)
    # monotone majorant (nonincreasing) as in tp.py
    OmPsi = np.maximum.accumulate(OmPsi[::-1])[::-1]; OmA = np.maximum.accumulate(OmA[::-1])[::-1]
    json.dump(dict(L=L.tolist(), OmA=OmA.tolist(), OmPsi=OmPsi.tolist()), open(_wout('g0_constants', 'tp_cert_0.05.json'), 'w'))
    try:
        tb = dict(np.load(_os.path.join(_R, 'g1_wall_lemma', 'heuristic', 'tp_tables2.npz')))
        for nm, v in (('OmA', OmA), ('OmPsi', OmPsi)):
            r = v/tb[nm]
            print('%-6s certified/numeric(tp.py, +0.2%%): max %.5f at l=%.3f, min %.5f ; numeric table is a valid majorant where ratio<=1: %s'
                  % (nm, r.max(), L[r.argmax()], r.min(), bool((r <= 1).all())))
    except Exception as ex:
        print('no tp table to compare:', ex)
    for l in [0.5, 1, 2, 4, 8, 16]:
        i = int(round(l/0.025)); print('l=%5.2f  OmA<=%.4e  OmPsi<=%.4e' % (l, OmA[i], OmPsi[i]))
