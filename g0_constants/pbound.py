"""Rigorous upper bounds for P(S') = sup_{z in H_eta} |Ai(z + S' e^{-i pi/6}) / psi(z)|, psi(z) = int_z^inf Ai.
On the boundary line (Phragmen-Lindelof reduction, see report):
   |Ai(z(t+S'))/psi(z(t))| = |F(z(t))| exp(-int_t^{t+S'} ut(tau) dtau),   ut = Re(-e^{-i pi/6} Ai'/Ai) >= at > 0.
Uses the certified per-box data of line_<eta>.json (|F| upper, ut lower) + explicit tails.
   C(c, S'min) := sup_t |F(z(t))| exp(-int_t^{t+S'min} (ut - c))   =>   P(S') <= C e^{-c S'} for all S' >= S'min  (needs c <= at_lo).
usage: pbound.py eta"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, math; sys.path.insert(0, _os.path.join(_R, 'g0_constants'))
from cl import *
eta = sys.argv[1]
d = json.load(open(_os.path.join(_CERT, 'g0_constants', 'line_%s.json' % eta))); boxes = d['boxes']; TL = d['TL']; TR = d['TR']
etaA = arb(eta)
# right tail data (constant lower bound of ut on [TR, inf), monotone)
tr = tail_side(TR, eta, 'right'); utR, coefR = tr[1], tr[2]
tlL = tail_side(TL, eta, 'left'); coefL = tlL[2]
alL = (etaA/arb(TL)).atan()/2; kL = (PI/4 + alL).cos(); eDL = arb(tlL[4])
def ut_left_lo(tau_abs_min):   # lower bound of ut at tau <= -TL with |tau| >= tau_abs_min (>= TL)
    return lb(arb(tau_abs_min).sqrt()*kL - eDL)
def Fup_tail(t_abs_max, coef):
    return ub(arb(coef)*((arb(t_abs_max)**2 + etaA**2)**(arb(1)/4)))
def build_cells(Sp):
    """cells (t1,t2,Fup,utlo) covering [-(TL+Sp+1), TR]: left-tail cells of width 1/64 + certified boxes."""
    cells = []
    TF = TL + math.ceil(Sp) + 1
    hc = 1/64; n = int(round((TF - TL)/hc))
    for i in range(n):
        a1 = -TF + i*hc; a2 = a1 + hc
        cells.append((a1, a2, Fup_tail(-a1, coefL), ut_left_lo(-a2)))
    cells += [(b[0], b[1], b[3], b[5]) for b in boxes]
    return cells, TF
def C_bound(c, Sp):
    """sup over t of |F| exp(-int_t^{t+Sp}(ut - c)), rigorous upper bound (float, rounded up generously)."""
    assert c <= min(d['at_lo'], utR)
    cells, TF = build_cells(Sp)
    g = [max(cl[3] - c, 0.0) for cl in cells]   # (ut - c) >= cl[3] - c ; and ut - c >= 0 everywhere (c <= at_lo)
    # prefix sums of g*width (lower bounds); floats rounded down via a 1e-12 relative safety factor
    pref = [0.0]
    for cl, gg in zip(cells, g): pref.append(pref[-1] + gg*(cl[1] - cl[0]))
    import bisect
    starts = [cl[0] for cl in cells]
    best = 0.0
    for i, cl in enumerate(cells):
        lo_end = cl[0] + Sp            # integrate over [cl[1], cl[0]+Sp] (subset of [t, t+Sp] for every t in the cell)
        if lo_end <= cl[1]: integ = 0.0
        else:
            # cells fully inside [cl[1], lo_end]
            j = i + 1; k = bisect.bisect_right(starts, lo_end) - 1
            while k >= j and cells[k][1] > lo_end: k -= 1
            integ = pref[k+1] - pref[j] if k >= j else 0.0
            if lo_end > TR:            # part beyond the compact: ut >= utR on [TR, lo_end]
                integ += max(utR - c, 0.0)*(lo_end - max(TR, cl[1]))
        integ *= (1 - 1e-9)
        best = max(best, cl[2]*math.exp(-integ))
    # far left: t <= -TF :  |F| exp(-Sp (kL sqrt(|t|-Sp) - eDL - c)), decreasing in |t| (Sp kL sqrt|t| > 1) -> value at TF
    fl = ub(arb(Fup_tail(TF, coefL))*(-(arb(Sp)*(kL*(arb(TF) - Sp).sqrt() - eDL - c))).exp())
    assert lb(arb(Sp)*kL*arb(TF).sqrt()) > 1
    # right tail t >= TR: coefR (t^2+eta^2)^{1/4} exp(-Sp (sqrt t k_R - eD_R - c)) decreasing -> value at TR
    alR = (etaA/arb(TR)).atan()/2; kR = (PI/4 + alR).cos(); eDR = arb(tr[4])
    assert lb(arb(Sp)*kR*arb(TR).sqrt()) > 1
    fr = ub(arb(Fup_tail(TR, coefR))*(-(arb(Sp)*(kR*arb(TR).sqrt() - eDR - c))).exp())
    return max(best, fl, fr)*(1 + 1e-9)
if __name__ == '__main__':
    at = d['at_lo']
    print('eta=%s  at_lo=%.6f (box data), utR=%.4f' % (eta, at, utR))
    for Sp in [2, 6, 10, 11.5, 14, 18.2, 20, 23.6]:
        print("  S'=%5.1f  P(S') <= %.3e" % (Sp, C_bound(0.0, Sp)))
    for c in [0.30, 0.35, round(math.floor(at*1e4)/1e4, 4)]:
        if c > at: continue
        for Smin in [1, 2, 5, 10]:
            print("  c=%.4f S'min=%4.1f:  P(S') <= %.4f e^{-c S'}  for all S' >= S'min" % (c, Smin, C_bound(c, Smin)))
