"""Certification (Arb + rigorous window minimisation) of the crit3 two-point tables T3-T6 (eta=0.05), l = 0:0.025:16:
  T3 K_Pi(l)  = sup|Ai(Z+le) Ai(wb Z)|              = sup_t W_A(t,l) * |Ai Ai(wb.)|(z(t))
  T4 K_Pi1(l) = sup max(|Ai'(Z+le)Ai(wbZ)|, |Ai(Z+le)Ai'(wbZ)|)
              <= sup_t W_A(t,l) * max( rA(t+l) |AiAi(wb.)|(t),  |Ai Ai'(wb.)|(t) ),   rA = |Ai'/Ai|
  T5 Q(l)     = sup|Ai'(Z+le) psi(wb Z)|             = sup_t rA(t+l) W_A(t,l) |Ai(z(t))| |Ai(wb z(t))| / |F(z(-t))|
  T6 Q*(l)    = sup|Ai'(wb Z) psi(Z+le)|             = sup_t W_psi(t,l) |Ai'(wb z(t))| |Ai(z(t))| / |F(z(t))|
  with W_A(t,l) = exp(-int_t^{t+l} ut), W_psi = exp(-int_t^{t+l} uF)  (exact identities on the line; z(t)+le = z(t+l)),
  |F(wb z(t))| = |F(z(-t))| (reflection), sup over H_eta = sup over the line (Phragmen-Lindelof, Lemma F1).
  Cells: G0 boxes on [-20,30] (h=1/256) with Arb Airy enclosures; tail cells (h=1/16) on [-40,-20] and [30,40] with the
  Wronskian identity + Lemmas B, D; |t|>40 analytically (rates >= tail rates).  Output: T36_cert_0.05.json"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, math
sys.path.insert(0, _os.path.join(_R, 'g0_constants'))
from cl import *
from flint import acb, arb, ctx
import numpy as np
ctx.prec = 64
eta = '0.05'; ETA = 0.05
d = json.load(open(_os.path.join(_CERT, 'g0_constants', 'line_%s.json' % eta))); boxes = d['boxes']
TL, TR, TF = 20.0, 30.0, 40.0
EMc = acb(arb(3).sqrt()/2, -arb(1)/2); E3c = acb(arb(1)/2, arb(3).sqrt()/2); WB = acb(-arb(1)/2, -arb(3).sqrt()/2)
def ball(t1, t2):
    tc = (arb(t1) + arb(t2))/2; r = (arb(t2) - arb(t1))/2
    zc = acb(tc)*EMc - arb(eta)*E3c
    rad = ub(r) + 1e-15 + ub(zc.real.rad()) + ub(zc.imag.rad())
    return acb(arb(zc.real.mid(), rad), arb(zc.imag.mid(), rad))
def U(x): return ub(x)
# ---------------- compact cells: Airy enclosures ----------------
cells = []      # (t1, t2, P0, P1b, rA, Flo, ut, uF)
Flo_of = {}
for (t1, t2, Flo, Fup, ulo, utlo, *_r) in boxes:
    P0 = P1b = rA = 0.0; NSUB = 6
    for q in range(NSUB):
        a1 = t1 + (t2 - t1)*q/NSUB; a2 = t1 + (t2 - t1)*(q + 1)/NSUB
        Z = ball(a1, a2); ai, aip, _, _ = Z.airy(); bi, bip, _, _ = (WB*Z).airy()
        P0 = max(P0, U(abs(ai*bi))); P1b = max(P1b, U(abs(ai*bip))); rA = max(rA, U(abs(aip/ai)))
    cells.append([t1, t2, P0, P1b, rA, Flo, utlo, ulo])
    Flo_of[round(t1*256)] = Flo
# ---------------- tail cells ----------------
def tail_vals(t1, t2, side):
    """profiles on a tail cell (|t| >= T = min|t| on the cell)"""
    T = min(abs(t1), abs(t2)); Tm = max(abs(t1), abs(t2))
    tr = tail_side(T, eta, side); tre = tail_side(T, eta, 'left' if side == 'right' else 'right')
    eB, eD = tr[3], max(tr[4], tre[4])
    x = math.sqrt(T)
    P0 = 1/(2*math.pi*(2*x - 2*eD))*(1 + 1e-9)              # Wronskian identity (cert_pi.py)
    rmax = (Tm**2 + ETA**2)**0.25
    rA = (rmax + eD)*(1 + 1e-9)
    P1b = rA*P0                                               # |Ai'/Ai|(wb z) <= |z|^{1/2}+eps_D as well
    Flo = x/(1 + eB)*(1 - 1e-9)
    return [t1, t2, P0, P1b, rA, Flo, tr[1], tr[0]]
left = [tail_vals(-TF + k/16, -TF + (k + 1)/16, 'left') for k in range(int((TF - TL)*16))]
right = [tail_vals(TR + k/16, TR + (k + 1)/16, 'right') for k in range(int((TF - TR)*16))]
allc = left + cells + right
A = np.array(allc)
t1, t2, P0, P1b, rA, Flo, ut, uF = [A[:, k] for k in range(8)]
edges = np.concatenate([[t1[0]], t2])
# reflected |F| lower bound: F(z(-t)) for t in cell [t1,t2]  -> cell [-t2,-t1]
def Flo_reflect(a1, a2):
    k = round(-a2*256)
    if -a2 >= -TL - 1e-12 and -a1 <= TR + 1e-12 and k in Flo_of: return Flo_of[k]
    T = min(abs(a1), abs(a2)); side = 'left' if -a1 <= 0 else 'right'
    tr = tail_side(max(T, 1e-9), eta, side)
    return math.sqrt(T)/(1 + tr[3])*(1 - 1e-9)
FloR = np.array([Flo_reflect(a, b) for a, b in zip(t1, t2)])
rate_far = {'ut': min(tail_side(TF, eta, 'left')[1], tail_side(TF, eta, 'right')[1]),
            'uF': min(tail_side(TF, eta, 'left')[0], tail_side(TF, eta, 'right')[0])}
eDF = max(tail_side(TF, eta, 'left')[4], tail_side(TF, eta, 'right')[4])
P0F = 1/(2*math.pi*(2*math.sqrt(TF) - 2*eDF))
def window_min(rate, l, rL, rR):
    """for each cell: lower bound of min_{t in cell} int_t^{t+l} rate (rate: step lower bound; slopes rL, rR outside)"""
    Fe = np.concatenate([[0.0], np.cumsum(rate*np.diff(edges))])
    def Fc(x):
        v = np.interp(x, edges, Fe)
        v = np.where(x < edges[0], Fe[0] - rL*(edges[0] - x), v)
        return np.where(x > edges[-1], Fe[-1] + rR*(x - edges[-1]), v)
    cand = [t1, t2]
    # interior breakpoints: t = e - l in (t1,t2) for edges e
    j = np.searchsorted(edges, t1 + l, side='right')
    e = edges[np.clip(j, 0, len(edges) - 1)] - l
    cand.append(np.where((e > t1) & (e < t2), e, t1))
    vals = np.min([Fc(c + l) - Fc(c) for c in cand], axis=0)
    return vals*(1 - 1e-9)
def range_max(arr, a, b):
    """max of arr over cells intersecting [a,b] (vectorised); outside the cell range: analytic bound"""
    i0 = np.clip(np.searchsorted(edges, a, side='right') - 1, 0, len(arr) - 1)
    i1 = np.clip(np.searchsorted(edges, b, side='left') - 1, 0, len(arr) - 1)
    out = arr[i0].copy()
    for k in range(1, 20):
        out = np.maximum(out, arr[np.minimum(i0 + k, i1)])
    assert (i1 - i0 < 20).all()
    return out
def rA_point_far(x):   # x > TF
    return ((x**2 + ETA**2)**0.25 + eDF)*(1 + 1e-9)
L = np.arange(0, 16.0001, 0.025)
T3, T4, T5, T6 = [], [], [], []
rLut, rRut = rate_far['ut'], tail_side(TR, eta, 'right')[1]
rLuF, rRuF = rate_far['uF'], tail_side(TR, eta, 'right')[0]
for l in L:
    WA = np.exp(-window_min(ut, l, rLut, rRut))*(1 + 1e-9)
    WP = np.exp(-window_min(uF, l, rLuF, rRuF))*(1 + 1e-9)
    a, b = t1 + l, t2 + l
    rAs = np.where(b <= edges[-1], range_max(rA, a, np.minimum(b, edges[-1])), rA_point_far(b))
    k3 = (P0*WA).max(); k4 = np.maximum(rAs*WA*P0, WA*P1b).max()
    q5 = (rAs*WA*P0/FloR).max(); q6 = (WP*P1b/Flo).max()
    # |t| > TF: profiles <= values at TF (decreasing in |t|), rates >= far rates
    far = math.exp(-l*min(rLut, rRut)*(1 - 1e-9)); farp = math.exp(-l*min(rLuF, rRuF)*(1 - 1e-9))
    rAfar = rA_point_far(TF + 16 + l + 1)    # crude: |z(t+l)| bounded for... (only multiplies tiny far terms)
    T3.append(max(k3, P0F*far)); T4.append(max(k4, rAfar*P0F*far)); T5.append(max(q5, rAfar*P0F*far/(math.sqrt(TF)/2)))
    T6.append(max(q6, rAfar*P0F*farp/(math.sqrt(TF)/2)))
# min with the uniform certified bounds of Lemma F1(c): K_Pi <= Pi e^{-at l}, K_Pi1 <= Pi1 e^{-at l}, Q,Q* <= Pi1psi e^{-a l}
cp = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pi_0.05.json'))); cq = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pipsi.json')))
at_, a_ = 0.35346270, 0.44452699
T3 = np.minimum(T3, cp['Pi']*np.exp(-at_*L)*(1 + 1e-9)); T4 = np.minimum(T4, cp['Pi1']*np.exp(-at_*L)*(1 + 1e-9))
T5 = np.minimum(T5, cq['Pi1_psi']*np.exp(-a_*L)*(1 + 1e-9)); T6 = np.minimum(T6, cq['Pi1_psi']*np.exp(-a_*L)*(1 + 1e-9))
mono = lambda v: np.maximum.accumulate(np.array(v)[::-1])[::-1]
T3, T4, T5, T6 = map(mono, (T3, T4, T5, T6))
json.dump(dict(L=L.tolist(), KPi=T3.tolist(), KPi1=T4.tolist(), Q=T5.tolist(), Qs=T6.tolist()),
          open(_wout('g0_constants', 'T36_cert_0.05.json'), 'w'))
tb = dict(np.load(_os.path.join(_R, 'g1_wall_lemma', 'heuristic', 'tp_tables2.npz')))
for nm, v in (('KPi', T3), ('KPi1', T4), ('Q', T5), ('Qs', T6)):
    r = v/tb[nm]
    print('%-5s cert: l=0 %.5f l=1 %.5f l=4 %.3e l=16 %.3e | cert/num: max %.4f (l=%.3f) min %.4f ; int_0^16 cert %.4f'
          % (nm, v[0], v[40], v[160], v[-1], r.max(), L[r.argmax()], r.min(), np.sum(v[:-1])*0.025))
