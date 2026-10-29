"""F = Ai/psi at the t_f grid t = k/128, k in [k0, k1] (eta = 0.05): rigorous Arb point values (g0_constants/cl.point,
integral representation of Appendix A), stored as complex centres and radii.  usage: ftab.py k0 k1 out.npz"""
import sys, os, math, numpy as np, time
from multiprocessing import Pool
_R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, os.path.join(_R, 'g0_constants'))
from cl import point
def one(k):
    F, u, ut = point(repr(k/128), '0.05', prec=80)
    c = complex(float(F.real.mid()), float(F.imag.mid()))
    r = math.nextafter(float(F.real.rad()) + float(F.imag.rad()) + 1e-16*abs(c), math.inf)
    return k, c, r
if __name__ == '__main__':
    k0, k1 = int(sys.argv[1]), int(sys.argv[2]); out = sys.argv[3]
    t = time.time()
    with Pool(2) as P: res = P.map(one, range(k0, k1 + 1), chunksize=64)
    res.sort(); np.savez(out, k=np.array([r[0] for r in res]), F=np.array([r[1] for r in res]), rF=np.array([r[2] for r in res]))
    print('F table', len(res), 'points, %.1f s, max rel radius %.2e' % (time.time() - t, max(r[2]/abs(r[1]) for r in res)))
