"""T8 (Lemma F8, relative cone), rigorous (Arb), per triple, on beta-CELLS (no Lipschitz constants T9/T10 needed):
   ThetaBar[cell] >= sup_{beta in cell} sup_{t in R} | arg G(z_g(t), b_g) - arg G(z_1(t), b_1) |,
   z_1(t) = t e - eta_1 e^{i pi/3},           eta_1 = delta/kr - nu1^{2/3} beta^2,      b_1 = beta nu1^{1/3} e^{i pi/6}
   z_g(t) = (kr t) e - eta_g e^{i pi/3},      eta_g = delta - g^{2/3} beta^2 + g^{-1/3}(s1 - s_g) sigma^{-2/3},  b_g = beta g^{1/3} e^{i pi/6}
   kr = (nu1/g)^{1/3}.  The sigma-shift (s1-s_g) sigma^{-2/3} (Lemma F8 step 4) is put EXACTLY into eta_g, for all
   sigma >= sigma_min (interval in eta), so the theta_rel of F8 is ThetaBar + asin(1-bl1) + asin(1-blg) [+ asin eps_x].
   (Theta = arg G(kr z - g^{-1/3}e^{i pi/3} d, b_g) - arg G(z, b_1), d = (nu1-g)beta^2 + (s1-s_g)sigma^{-2/3}, is harmonic
   and bounded on H_{eta_1} since both lines have eta <= 0.1 < 0.8895 (G0 W1: L != 0, arg G in a fixed strip); PL.)
   Compact part t in [-TL, TR]: second-order (affine) Arb boxes of the G0 integral representation (g0lib.box_affine,
   generalised to eta-intervals), branch and bound on t-cells.  Tails |t| >= T: Lemma B (G = c/(1+eps), c = sqrt z + b)
   with Arb balls on geometric cells up to 1e6; beyond: explicit bound.
usage: cert_T8.py N J G sigma_min  -> T8_cert_<N>_<J>_<G>.json"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, math
sys.path.insert(0, _os.path.join(_R, 'g0_sector_lemma'))
from g0lib import *
from multiprocessing import Pool
ctx.prec = 64
DELTA = arb('0.05')

def box_aff_eta(t1, t2, eta_c, h_eta, r1, r2, amin=0.6):
    """g0lib.box_affine with z-centre on the line eta_c and |eta - eta_c| <= h_eta (the eta-variation enters the
    remainder through |dlogG/dz| h_eta and the second-order radius rZ).  Returns arg interval (lo, hi) or None."""
    tc = 0.5*(t1 + t2); ht = max(ub(arb(t2) - arb(tc)), ub(arb(tc) - arb(t1)))
    rc = 0.5*(r1 + r2); hr = max(ub(arb(r2) - arb(rc)), ub(arb(rc) - arb(r1)))
    zt = line_point(tc, eta_c); zc_c = complex(float(zt.real.mid()), float(zt.imag.mid()))
    zc = acb(zc_c.real, zc_c.imag); ez = arb(ub(abs(zt - zc))) + arb(h_eta)
    bt = acb(rc)*EP; bc_c = complex(float(bt.real.mid()), float(bt.imag.mid()))
    bc = acb(bc_c.real, bc_c.imag); eb = arb(ub(abs(bt - bc)))
    rZ = arb(ht) + ez; rB = arb(hr) + eb
    w0 = choose_w0(zc_c, amin); a = w0.real
    gam = arb(ub(abs((w0*w0 - zc).imag))) + rZ
    reb = max(min(lb(bc.real), lb(arb(r1)*SQ3H)), 0.0)
    ab = a + arb(reb)
    m0, m1, m2 = moments(a, gam)
    N0, N1, D0, Db, c = centre_data(zc, bc, w0)
    Dz = c*D0 - N0
    G0 = N0/D0
    aN0 = abs(N0); aD0 = abs(D0)
    if lb(aN0) <= 0 or lb(aD0) <= 0: return None
    RN = rZ*rZ*m2/2
    RD = (rZ*rZ*m2/ab + 2*rZ*rB*m1/(ab*ab) + 2*rB*rB*m0/(ab**3))/2
    X = (abs(N1)*rZ + RN)/aN0
    Y = (abs(Dz)*rZ + abs(Db)*rB + RD)/aD0
    if ub(Y) >= 1: return None
    lz = N1/N0 - Dz/D0
    k1 = lz*EM; k2 = -(Db/D0)*EP
    rem = abs(lz)*ez + abs(Db/D0)*eb + RN/aN0 + RD/aD0 + (Y*Y + X*Y)/(1 - Y)
    # G = G0 (1+u), u = k1 dt + k2 drho + e (|dt|<=ht, |drho|<=hr real, |e|<=rem):  |arg(1+u)| <= atan(Im/(1-|Re|))
    Aim = abs(k1.imag)*ht + abs(k2.imag)*hr + rem
    Bre = abs(k1.real)*ht + abs(k2.real)*hr + rem
    if ub(Bre) >= 0.9: return None
    phi = (Aim/(1 - Bre)).atan()
    ag = G0.arg()
    return lb(ag - phi), ub(ag + phi), ub(abs(k1)*ht), ub(abs(k2)*hr), ub(abs(lz))

def tail_arg(t1, t2, eta_c, h_eta, r1, r2):
    """Lemma B on a tail cell: arg G in arg(sqrt z + b) +- asin(eps).  None if not applicable."""
    tc = 0.5*(t1 + t2); ht = max(ub(arb(t2) - arb(tc)), ub(arb(tc) - arb(t1)))
    zt = line_point(tc, eta_c)
    rad = arb(ub(arb(ht) + arb(h_eta)))
    zb = zt + acb(arb(0, ub(rad)), arb(0, ub(rad)))
    w = zb.sqrt()
    if not lb(w.real) > 0: return None
    kap = lb(w.real/abs(w)); r = lb(abs(w))
    if kap <= 0: return None
    den = arb(r)**3 - arb(5)/48*arb(kap)**(-arb(7)/2)
    if not lb(den) > 0: return None
    eps = tailB_eps(r, kap)
    if not ub(eps) < 0.5: return None
    b = acb(arb(0.5*(r1 + r2), 0.5*(r2 - r1) + 1e-15))*EP
    cc = w + b
    ac = cc.arg(); ae = eps.asin()
    return lb(ac - ae), ub(ac + ae), 0.0, 0.0, 0.0

def far_bound(T, eta1, etag, kr, rho1, rhog):
    """|t| >= T (T huge): |Theta| <= [atan(|eta1|/T) + atan(|etag|/(kr T))]/2 + sum_nu [asin(rho_nu/r_nu) + asin eps_nu]."""
    T = arb(T); tot = arb(0)
    for (eta, tt, rho) in ((eta1, T, rho1), (etag, kr*T, rhog)):
        tot += (arb(abs(eta))/tt).atan()/2
        r = tt.sqrt()                          # |z| >= |t|
        kap = arb('0.25')                      # Re sqrt z/|sqrt z| >= cos(5pi/12 + small) > 0.25 for |t| >= 1e6
        tot += (arb(rho)/r).asin() + tailB_eps(r, kap).asin()
    return ub(tot)

def cell_sup(args):
    (be_a, be_b, nu1s, gs, dss, smin, TL, TR, tol) = args
    nu1 = arb(nu1s); g = arb(gs); ds = arb(dss); kr = (nu1/g)**(ONE3)
    n23, g23, g13, n13 = nu1**(2*ONE3), g**(2*ONE3), g**ONE3, nu1**ONE3
    shm = g**(-ONE3)*ds*arb(smin)**(-2*ONE3)                  # eta_g shift in [0, shm] (or [shm, 0])
    krl, kru = lb(kr), ub(kr)
    def ranges(b1, b2, u1, u2):
        B = arb(b1).union(arb(b2)); U = arb(u1).union(arb(u2))
        e1 = DELTA/kr - n23*B*B; eg = DELTA - g23*B*B + U*shm
        out = []
        for e in (e1, eg):
            lo, hi = lb(e), ub(e); c = 0.5*(lo + hi); h = max(ub(arb(hi) - arb(c)), ub(arb(c) - arb(lo)))
            out.append((c, h))
        r1 = (lb(B*n13), ub(B*n13)); rg = (lb(B*g13), ub(B*g13))
        return out[0], out[1], r1, rg
    def tg(t1, t2):
        a = t1*(krl if t1 >= 0 else kru); b = t2*(kru if t2 >= 0 else krl)
        return math.nextafter(a, -math.inf), math.nextafter(b, math.inf)
    def evalcell(c, tail=False):
        t1, t2, b1, b2, u1, u2 = c
        (c1, h1), (cg, hg), r1, rg = ranges(b1, b2, u1, u2)
        ga, gb = tg(t1, t2)
        f = tail_arg if tail else box_aff_eta
        A = f(t1, t2, repr(c1), h1, *r1); B = f(ga, gb, repr(cg), hg, *rg)
        if A is None or B is None: return None
        U = max(B[1] - A[0], A[1] - B[0]); L = abs(0.5*(B[0] + B[1]) - 0.5*(A[0] + A[1]))
        # split heuristics: t, beta (rho + eta_beta), shift
        hb1 = h1; hbg = ub(g23*(arb(b2)**2 - arb(b1)**2))/2
        ct = max(A[2], B[2]); cb = max(A[3] + A[4]*h1, B[3] + B[4]*hbg); cs = B[4]*ub(abs(shm))*(u2 - u1)/2
        return U, L, ct, cb, cs
    w0 = 0.5
    pend = []; t = -TL
    while t < TR - 1e-12:
        c = (t, min(t + w0, TR), be_a, be_b, 0.0, 1.0); pend.append((c, evalcell(c))); t += w0
    best_L = max([r[1] for c, r in pend if r] + [0.0])
    final = 0.0; worst = None; nev = len(pend)
    while pend:
        c, r = pend.pop()
        t1, t2, b1, b2, u1, u2 = c
        if r is not None and r[0] <= best_L + tol:
            if r[0] > final: final, worst = r[0], c
            continue
        if nev > 25000 or (t2 - t1 < 1/1024 and b2 - b1 < 1e-5 and u2 - u1 < 1e-3):
            return dict(cell=(be_a, be_b), ok=False, why='box', at=c, nev=nev, L=best_L)
        if r is None:
            (c1_, h1_), (cg_, hg_), _, _ = ranges(b1, b2, u1, u2)
            hbg_ = ub(g23*(arb(b2)**2 - arb(b1)**2))/2; hs_ = ub(abs(shm))*(u2 - u1)/2
            ht_ = (t2 - t1)/2*max(1.0, krl)
            cand = [(ht_ if t2 - t1 > 1/1024 else -1, 't'), (max(h1_, hbg_, (b2 - b1)*1.3) if b2 - b1 > 1e-5 else -1, 'b'),
                    (hs_ if u2 - u1 > 1e-3 else -1, 's')]
            dim = max(cand)[1]
        else:
            ct, cb, cs = r[2], r[3], r[4]
            if t2 - t1 < 1/1024: ct = -1
            if b2 - b1 < 1e-5: cb = -1
            if u2 - u1 < 1e-3: cs = -1
            dim = 't' if ct >= max(cb, cs) else ('b' if cb >= cs else 's')
        if dim == 't': m = 0.5*(t1 + t2); kids = [(t1, m, b1, b2, u1, u2), (m, t2, b1, b2, u1, u2)]
        elif dim == 'b': m = 0.5*(b1 + b2); kids = [(t1, t2, b1, m, u1, u2), (t1, t2, m, b2, u1, u2)]
        else: m = 0.5*(u1 + u2); kids = [(t1, t2, b1, b2, u1, m), (t1, t2, b1, b2, m, u2)]
        for k in kids:
            rr = evalcell(k); nev += 1
            if rr: best_L = max(best_L, rr[1])
            pend.append((k, rr))
    comp = final
    # tails (whole beta-cell and shift range; cells subdivided in beta if Lemma B too crude)
    tl = 0.0
    for side in (1, -1):
        T = TR if side == 1 else TL
        while T < 1e6:
            T2 = min(T*1.08, 1e6)
            a, b = (T, T2) if side == 1 else (-T2, -T)
            nb = 1; rmax = None
            while nb <= 64:
                vals = [evalcell((a, b, be_a + (be_b - be_a)*k/nb, be_a + (be_b - be_a)*(k + 1)/nb, 0.0, 1.0), tail=True) for k in range(nb)]
                if all(v is not None for v in vals):
                    rmax = max(v[0] for v in vals)
                    if rmax <= max(comp, best_L + tol) or nb >= 16: break
                nb *= 2
            if rmax is None: return dict(cell=(be_a, be_b), ok=False, why='tail', at=(a, b))
            tl = max(tl, rmax); T = T2
    (c1, h1), (cg, hg), r1, rg = ranges(be_a, be_b, 0.0, 1.0)
    far = far_bound(1e6, abs(c1) + h1, abs(cg) + hg, kr, r1[1], rg[1])
    return dict(cell=(be_a, be_b), ok=True, Theta=max(comp, tl, far), comp=comp, tail=tl, far=far, at=worst, nev=nev, L=best_L)

ONE3 = arb(1)/3
if __name__ == '__main__':
    N, J, G, smin = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    bmax = float(sys.argv[5]) if len(sys.argv) > 5 else 3.3
    Na, Ja, Ga = arb(N), arb(J), arb(G)
    g = Ga/Ja; nu1 = 1 + Na
    kapc = 4*Na/Ja; cc = 2*Na/(1 + Na - g)
    ds = kapc*(cc - 1)                                     # s1 - s_g
    # store exact-ish decimal strings (enclosing balls are re-created from strings; widen by tiny relative amount)
    nu1s, gs, dss = nu1.mid().str(20, radius=False), g.mid().str(20, radius=False), ds.mid().str(20, radius=False)
    assert abs(float(g.mid()) - float(gs)) < 1e-15
    ED = [0, 0.002, 0.005, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.8, 1, 1.3, 1.7, 2.2, 3, 4, 5, 6.5]
    cells = []
    for i in range(len(ED) - 1):
        a, b = ED[i], ED[i + 1]
        if a >= bmax: break
        # sub-cells: eta_1 width <= 0.25
        cells.append((a, b, i))
    tol = lambda b: 0.008 if b <= 0.1 else (0.015 if b <= 0.3 else (0.04 if b <= 1 else 0.1))
    jobs = [(a, b, nu1s, gs, dss, smin, 100.0, 60.0, tol(b)) for (a, b, i) in cells]
    jobs = jobs[::-1]; cells = cells[::-1]          # expensive cells first
    import time; T0 = time.time()
    with Pool(2) as P:
        res = P.map(cell_sup, jobs, chunksize=1)
    out = dict(NJG=[N, J, G], sigma_min=float(smin), edges=ED, cells=[])
    for (a, b, i), r in zip(cells, res):
        out['cells'].append(dict(r, edge=i))
    # per driver edge: max over sub-cells
    tab = {}
    for c in out['cells']:
        if not c['ok']: tab[c['edge']] = math.inf; continue
        tab[c['edge']] = max(tab.get(c['edge'], 0.0), c['Theta'])
    out['Theta_edge'] = [[ED[i], ED[i + 1], tab[i]] for i in sorted(tab)]
    out['seconds'] = time.time() - T0
    json.dump(out, open(_wout('g0_constants', 'T8_cert_%s_%s_%s.json' % (N, J, G)), 'w'), indent=1)
    for r in out['Theta_edge']: print('beta in [%g,%g]: Theta <= %.4f' % tuple(r))
    for c in out['cells']:
        if not c['ok']: print('FAIL', c)
    print('seconds', out['seconds'])
