"""pior rho pontual na reta Re = -delta' com as linhas (m_L, r_4') Airy exatas (e LG / Airy LR como alternativas)"""
import sys, numpy as np, ap_arb as AA, ap_airy3 as A3
from flint import arb
D = -0.0396
def rho_pt(sg, b, y):
    best = 9
    r1 = AA.box_eval(arb(sg), arb(b), arb(D), arb(y)); best = r1['rho'] if r1.get('ok') else 9
    if best < 0.3: return best
    for rows in (('mL', 'r4'), ('mR', 'mL')):
        if rows == ('mR', 'mL') and b <= 0: continue
        r = A3.airy3_safe(arb(sg), arb(b), arb(D), arb(y), rows=rows)
        if 'rho' in r and r.get('rA', 9) < 1: best = min(best, r['rho'])
    return best
def worst(sg, b, stop=None):
    S = 2*sg**(1/3); ys = list(np.arange(S/2, -3, -0.25)) + [-4, -6, -9, -13, -20]; w = (0, None)
    for y in ys:
        v = rho_pt(sg, b, y)
        if v > w[0]: w = (v, y)
        if stop and w[0] > stop: break
    return w
if __name__ == '__main__':
    sg = float(sys.argv[1])
    for b in map(float, sys.argv[2:]):
        w = worst(sg, b); print('sigma=%g beta=%g pior rho na reta=%.3f em Im=%s' % (sg, b, w[0], w[1]), flush=True)
