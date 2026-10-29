import sys, math, ap_arb as AA, fa_table as FT
from flint import arb
dlt = 0.0397
def pts(Ls, S):
    out = []
    for re in [Ls, 1.1*Ls, 1.3*Ls, 2*Ls]:
        for k in range(17): out.append(complex(re, S/2 - k*(S/2 + 2*Ls)/16))
    for re in [-dlt, 0, Ls/4, Ls/2, 3*Ls/4, Ls]:
        for im in [-Ls, -1.3*Ls, -2*Ls]: out.append(complex(re, im))
    return out
def _main():
  sg = float(sys.argv[1])
  for Ls in map(float, sys.argv[2:]):
      S = 2*sg**(1/3); bmax = FT.alpha_max(sg*1.0001)*sg**(-1/3); w = (0, None)
      for b in [0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.6, 1, bmax]:
          if b > bmax: continue
          for L in pts(Ls, S):
              o = AA.box_eval(arb(sg), arb(b), arb(L.real), arb(L.imag))
              v = o.get('rho', 9) if o.get('why') is None else 9
              if v > w[0]: w = (v, (b, L))
      print('sigma=%g Lam*=%g worst rho=%.3f at %s' % (sg, Ls, w[0], w[1]), flush=True)

if __name__ == '__main__': _main()