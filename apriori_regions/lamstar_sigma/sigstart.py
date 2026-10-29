"""para cada Lambda* inteiro, o menor sigma (bissecao geometrica) com pior rho pontual <= THR (amostra de pontos de scan_pts)"""
import sys, math, json, ap_arb as AA, fa_table as FT
from flint import arb
from scan_pts import pts
THR = float(sys.argv[1])
def worst(sg, Ls):
    S = 2*sg**(1/3); bmax = FT.alpha_max(sg*1.0001)*sg**(-1/3); w = 0
    for b in [0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.6, 1, bmax]:
        if b > bmax: continue
        for L in pts(Ls, S):
            o = AA.box_eval(arb(sg), arb(b), arb(L.real), arb(L.imag))
            v = o.get('rho', 9) if o.get('why') is None else 9
            w = max(w, v)
            if w > THR: return w
    return w
res = {}
for Ls in map(float, sys.argv[2:]):
    lo, hi = 3.7498, 288.0
    if worst(lo, Ls) <= THR: res[Ls] = lo; print(Ls, 'ok ja em', lo, flush=True); continue
    if worst(hi, Ls) > THR: res[Ls] = None; print(Ls, 'nunca', flush=True); continue
    for _ in range(9):
        m = math.sqrt(lo*hi)
        if worst(m, Ls) <= THR: hi = m
        else: lo = m
    res[Ls] = hi; print(Ls, 'sigma_start=%.4f' % hi, flush=True)
json.dump(res, open('sigstart_%s.json' % sys.argv[2], 'w'))
