"""Rigorous line tables for route 4.  Line z(tau) = tau e - eta e^{i pi/3}, e = e^{-i pi/6}, eta = 0.05, grid tau_i = i/128.
Point values a_i = Ai(z(tau_i)), a1_i = Ai'(z(tau_i)) as complex float centres with radii (Arb, prec 80, outward), and
cell suprema over [tau_i, tau_{i+1}]: Amax_i >= sup|Ai|, A1max_i >= sup|Ai'|, Zmax_i >= sup|z| (Arb on the cell ball)."""
import math, numpy as np, json, sys, time
from flint import arb, acb, ctx
ctx.prec = 128
PI = arb.pi(); I = acb(0, 1)
EMa = (-I*PI/6).exp(); E3a = (I*PI/3).exp(); ETA = arb('0.05')
DT = 1/128.0
def up(x):   # float upper bound of a real ball / nonneg
    v = float(x.upper()) if hasattr(x, 'upper') else float(x)
    return math.nextafter(v, math.inf)
def build(imax):
    n = 2*imax + 1
    a = np.zeros(n, complex); a1 = np.zeros(n, complex); ra = np.zeros(n); ra1 = np.zeros(n)
    Am = np.zeros(n); A1m = np.zeros(n); Zm = np.zeros(n)
    half = arb(DT)/2
    for k in range(n):
        i = k - imax
        t = arb(i)/128
        z = acb(t)*EMa - ETA*E3a
        A, Ap = z.airy()[:2]
        a[k] = complex(float(A.real.mid()), float(A.imag.mid())); a1[k] = complex(float(Ap.real.mid()), float(Ap.imag.mid()))
        ra[k] = up(abs(A - acb(a[k].real, a[k].imag)))*(1 + 1e-15) + 1e-300
        ra1[k] = up(abs(Ap - acb(a1[k].real, a1[k].imag)))*(1 + 1e-15) + 1e-300
        # cell [t, t + DT]: ball centre t + DT/2, radius DT/2 in t -> z ball (|e| = 1)
        tc = t + half
        zc = acb(tc)*EMa - ETA*E3a
        zb = acb(arb(zc.real.mid(), up(half) + up(abs(zc.real.rad())) if False else zc.real), arb(zc.imag))
        # enclose the segment: add radius DT/2 to both real and imaginary parts (covers |dz| <= DT/2)
        zb = acb(arb(zc.real.mid(), float(half.upper()) + float(zc.real.rad())), arb(zc.imag.mid(), float(half.upper()) + float(zc.imag.rad())))
        Ab, Apb = zb.airy()[:2]
        Am[k] = up(abs(Ab)); A1m[k] = up(abs(Apb)); Zm[k] = up(abs(zb))
    return dict(imax=imax, a=a, a1=a1, ra=ra, ra1=ra1, Am=Am, A1m=A1m, Zm=Zm)
if __name__ == '__main__':
    imax = int(sys.argv[1]) if len(sys.argv) > 1 else 100*128
    t = time.time(); d = build(imax)
    np.savez_compressed(sys.argv[2] if len(sys.argv) > 2 else 'line_tab.npz', **d)
    print('built', len(d['a']), 'points in %.1f s; max rel radius %.2e; max cell/point ratio %.3f' % (time.time() - t,
          np.max(d['ra']/np.maximum(abs(d['a']), 1e-300)), np.max(d['Am'][:-1]/np.maximum(abs(d['a'][:-1]), 1e-300))))
