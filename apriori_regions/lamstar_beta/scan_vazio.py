"""pior rho pontual na regiao inteira {Re >= -delta', Im <= S/2, |Lambda_1| < 100} para (sigma, beta)"""
import sys, math, ap_arb as AA
from flint import arb
D = 0.0396
def full_pts(S):
    xs = [-D, 0, 0.1, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 20]
    ys = [S/2 - k*S/16 for k in range(9)] + [-0.25, -0.5, -1, -2, -3, -4, -6, -8, -12, -20]
    return [complex(x, y) for x in xs for y in ys]
def worst(sg, b, stop=None):
    S = 2*sg**(1/3); w = (0, None)
    for L in full_pts(S):
        r = AA.box_eval(arb(sg), arb(b), arb(L.real), arb(L.imag))
        v = r.get('rho', 9) if r.get('why') is None else 9
        if v > w[0]: w = (v, L)
        if stop and w[0] > stop: break
    return w
if __name__ == '__main__':
    sg = float(sys.argv[1])
    for b in map(float, sys.argv[2:]):
        w = worst(sg, b); print('sigma=%g beta=%g pior rho=%.3f em %s' % (sg, b, w[0], w[1]), flush=True)
