"""
G0: rigorous (ball-arithmetic) verification of the sector lemma W_beta.
Core identity (contour through t = w0, any w0 with Re w0 > 0):
   Ai(z)   = (1/2pi) e^{w0^3/3 - z w0} N(z),      N(z)   = int_R Phi(s;z) ds
   L(z,b)  = (1/2pi) e^{w0^3/3 - z w0} D(z,b),    D(z,b) = int_R Phi(s;z)/(w0+b+is) ds
   Phi(s;z)= exp(-w0 s^2 - i s^3/3 + i (w0^2 - z) s)
so  G(z,b) = Ai(z)/L(z,b) = N(z)/D(z,b).
All quantities below are arb/acb balls (python-flint); every bound is rigorous.
"""
import math, cmath
from flint import arb, acb, ctx
ctx.prec = 64
I = acb(0, 1)
PI = arb.pi()
EM = (-I*PI/6).exp()      # e^{-i pi/6}
EP = (I*PI/6).exp()       # e^{+i pi/6}
E3 = (I*PI/3).exp()       # e^{+i pi/3}
SQ3H = arb(3).sqrt()/2    # cos(pi/6)

def ub(x):
    """rigorous float upper bound of real arb x (rounded outward)"""
    return math.nextafter(float(x.upper()), math.inf)
def lb(x):
    """rigorous float lower bound of real arb x (rounded outward)"""
    return math.nextafter(float(x.lower()), -math.inf)

def A(v):
    """exact/enclosing arb from a decimal string or float"""
    return arb(v) if not isinstance(v, str) else arb(v)

def disk(c, r):
    """acb rectangle containing the closed disk B(c, r) (c acb, r arb/float >= 0)."""
    r = arb(r)
    rr = arb(0, ub(r))
    return c + acb(rr, rr)

def moments(a, g):
    """upper bounds for m_k = int_R |s|^k exp(-a s^2 + g|s|) ds, k=0,1,2 (a>0, g>=0)."""
    a = arb(a); g = arb(g)
    I0 = (PI/a).sqrt()/2 * (g*g/(4*a)).exp() * (1 + (g/(2*a.sqrt())).erf())
    I1 = (1 + g*I0)/(2*a)
    I2 = (I0 + g*I1)/(2*a)
    return 2*I0, 2*I1, 2*I2

def ND_center(zc, bc, w0):
    """N(zc), D(zc,bc) as acb balls, including the truncation error of |s|>S."""
    a = float(w0.real.mid())
    q = I*(w0*w0 - zc)
    g0 = ub(abs((w0*w0 - zc).imag))
    # choose S with a S^2 - g0 S >= 46
    S = (g0 + math.sqrt(g0*g0 + 4*a*46))/(2*a)
    S = math.ceil(S*8)/8
    fN = lambda x, an: (-(w0*x*x) - I*x*x*x/3 + q*x).exp()
    fD = lambda x, an: (-(w0*x*x) - I*x*x*x/3 + q*x).exp()/(w0 + bc + I*x)
    N = acb.integral(fN, -S, S)
    D = acb.integral(fD, -S, S)
    aa = arb(a); Sa = arb(S)
    tail = 2*(-(aa*Sa*Sa) + arb(g0)*Sa).exp()/(2*aa*Sa - arb(g0))
    N = disk(N, tail)
    D = disk(D, tail/aa)
    return N, D

def choose_w0(zc_c, amin):
    w = cmath.sqrt(zc_c)
    if w.real < amin:
        w = complex(amin, w.imag)
    return acb(w.real, w.imag)

def line_point(t, eta):
    """exact acb enclosure of z(t) = t e^{-i pi/6} - eta e^{i pi/3}  (t, eta floats/arbs)."""
    return acb(t)*EM - arb(eta)*E3

def box_G(t1, t2, r1, r2, eta, amin=0.6):
    """Enclosure of G(z,b) for z on the (shifted) line, t in [t1,t2], rho in [r1,r2].
    Returns (Gball, contrib_t, contrib_rho)."""
    tc = 0.5*(t1+t2); ht = 0.5*(t2-t1)
    rc = 0.5*(r1+r2); hr = 0.5*(r2-r1)
    zt = line_point(tc, eta); zc_c = complex(float(zt.real.mid()), float(zt.imag.mid()))
    zc = acb(zc_c.real, zc_c.imag)
    rZ = arb(ht) + arb(ub(abs(zt - zc)))
    bt = acb(rc)*EP; bc_c = complex(float(bt.real.mid()), float(bt.imag.mid()))
    bc = acb(bc_c.real, bc_c.imag)
    rB = arb(hr) + arb(ub(abs(bt - bc)))
    w0 = choose_w0(zc_c, amin)
    a = w0.real
    gam = arb(ub(abs((w0*w0 - zc).imag))) + rZ
    # inf Re b along the segment [bc, b]: >= min(Re bc, r1 cos(pi/6)) and >= 0
    reb = min(lb(bc.real), lb(arb(r1)*SQ3H)); reb = max(reb, 0.0)
    ab = a + arb(reb)
    m0, m1, m2 = moments(a, gam)
    N, D = ND_center(zc, bc, w0)
    eN = rZ*m1
    eD = rZ*m1/ab + rB*m0/(ab*ab)
    Nb = disk(N, eN); Db = disk(D, eD)
    G = Nb/Db
    absN = abs(N); absD = abs(D)
    ct = ub(rZ*m1/absN + rZ*m1/(ab*absD)) if lb(absN) > 0 and lb(absD) > 0 else float('inf')
    cr = ub(rB*m0/(ab*ab*absD)) if lb(absD) > 0 else float('inf')
    return G, ct, cr, (N, D)

def cone_ok(G, th1, th2, mmin=None):
    """rigorous check: arg G in (th1*pi, th2*pi) via two half-plane tests; |G|>mmin."""
    u1 = (G*(-I*PI*arb(th1)).exp()).imag     # must be > 0
    u2 = (G*(-I*PI*arb(th2)).exp()).imag     # must be < 0
    ok = lb(u1) > 0 and ub(u2) < 0
    if mmin is not None:
        ok = ok and lb(abs(G)) > mmin
    return ok

# ---------------- Tail lemma B (|z| large, all rho >= 0) ----------------
def tailB_eps(r, kappa):
    """|eps| bound for |w|=r, Re w >= kappa*r:  eps <= (3/4) kappa^{-5/2} / (r^3 - (5/48) kappa^{-7/2})."""
    r = arb(r); k = arb(kappa)
    den = r**3 - arb(5)/48*k**(-arb(7)/2)
    assert lb(den) > 0
    return arb(3)/4*k**(-arb(5)/2)/den

# ---------------- Tail lemma A (rho >= R, z in a t-box) ----------------
def tailA_box(t1, t2, R, eta, amin=0.6):
    """For t in [t1,t2], rho in [R, inf): returns (arg interval/pi lower, upper, |G| lower)."""
    tc = 0.5*(t1+t2); ht = max(ub(arb(t2)-arb(tc)), ub(arb(tc)-arb(t1)))
    zt = line_point(tc, eta); zc_c = complex(float(zt.real.mid()), float(zt.imag.mid()))
    zc = acb(zc_c.real, zc_c.imag)
    rZ = arb(ht) + arb(ub(abs(zt - zc)))
    w0 = choose_w0(zc_c, amin); a = w0.real
    gam = arb(ub(abs((w0*w0 - zc).imag))) + rZ
    m0, m1, m2 = moments(a, gam)
    N, _ = ND_center(zc, acb(0), w0)
    Nlow = abs(N) - rZ*m1
    if lb(Nlow) <= 0: return None
    delta = m1/((a + arb(R)*SQ3H)*Nlow)
    dl = ub(delta)
    if dl >= 1: return None
    s = arb(0).union(1/arb(R))               # s = 1/rho in [0, 1/R]  (exact hull; audit fix)
    u = 1 + w0*EM*s
    argu = u.arg()
    asd = arb(dl).asin()
    lo = (PI/6 + argu - asd)/PI; hi = (PI/6 + argu + asd)/PI
    absc0 = arb(R) - abs(w0)                 # |c0| >= rho - |w0|
    Glow = absc0/(1 + arb(dl))
    return lb(lo), ub(hi), lb(Glow)

def tailB_check(eta, TL, TR, th1, th2, mmin):
    """Tail lemma B on the shifted line, |t|>=T, all rho>=0.  Returns list of (side, lo, hi, Glow, ok)."""
    eta = arb(eta); TL = arb(TL); TR = arb(TR); th1 = arb(th1); th2 = arb(th2)
    out = []
    for side, T in (('right', TR), ('left', TL)):
        al = (eta/T).atan()/2
        r = T.sqrt()
        if side == 'right':
            k = (PI/12 + al).cos(); lo_c = -PI/12 - al; hi_c = PI/6
        else:
            k = (5*PI/12 + al).cos(); lo_c = PI/6; hi_c = 5*PI/12 + al
        eps = tailB_eps(r, k)
        assert ub(eps) < 1
        asn = eps.asin()
        lo = (lo_c - asn)/PI; hi = (hi_c + asn)/PI
        Glow = r/(1 + eps)
        ok = lb(lo - th1) > 0 and ub(hi - th2) < 0 and (mmin is None or lb(Glow) > mmin)
        out.append((side, lb(lo), ub(hi), lb(Glow), ub(eps), ok))
    return out

# ---------------- second-order (affine) box enclosure ----------------
def centre_data(zc, bc, w0):
    """N0, N1=N'(zc), D0, Db=dD/db at (zc,bc); balls incl. truncation error."""
    a = float(w0.real.mid())
    q = I*(w0*w0 - zc)
    g0 = ub(abs((w0*w0 - zc).imag))
    S = (g0 + math.sqrt(g0*g0 + 4*a*46))/(2*a)
    S = math.ceil(S*8)/8
    c = w0 + bc
    ph = lambda x: (-(w0*x*x) - I*x*x*x/3 + q*x).exp()
    N0 = acb.integral(lambda x, an: ph(x), -S, S)
    N1 = acb.integral(lambda x, an: -I*x*ph(x), -S, S)
    D0 = acb.integral(lambda x, an: ph(x)/(c + I*x), -S, S)
    Db = acb.integral(lambda x, an: -ph(x)/((c + I*x)**2), -S, S)
    aa = arb(a); Sa = arb(S); gg = arb(g0)
    e0 = 2*(-(aa*Sa*Sa) + gg*Sa).exp()/(2*aa*Sa - gg)
    e1 = 2*Sa*(-(aa*Sa*Sa) + gg*Sa).exp()/(2*aa*Sa - gg - 1/Sa)
    return disk(N0, e0), disk(N1, e1), disk(D0, e0/aa), disk(Db, e0/(aa*aa)), c

def box_affine(t1, t2, r1, r2, eta, amin=0.6):
    """Second-order enclosure on the box t in [t1,t2], rho in [r1,r2] of the (shifted) line.
    G(z(t), rho e^{i pi/6}) in G0*(1 + k1*dt + k2*drho + B(0,rem)), |dt|<=ht, |drho|<=hr.
    Returns dict."""
    tc = 0.5*(t1+t2); ht = max(ub(arb(t2)-arb(tc)), ub(arb(tc)-arb(t1)))
    rc = 0.5*(r1+r2); hr = max(ub(arb(r2)-arb(rc)), ub(arb(rc)-arb(r1)))
    zt = line_point(tc, eta); zc_c = complex(float(zt.real.mid()), float(zt.imag.mid()))
    zc = acb(zc_c.real, zc_c.imag); ez = arb(ub(abs(zt - zc)))
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
    if lb(aN0) <= 0 or lb(aD0) <= 0:
        return None
    RN = rZ*rZ*m2/2
    RD = (rZ*rZ*m2/ab + 2*rZ*rB*m1/(ab*ab) + 2*rB*rB*m0/(ab**3))/2
    X = (abs(N1)*rZ + RN)/aN0
    Y = (abs(Dz)*rZ + abs(Db)*rB + RD)/aD0
    if ub(Y) >= 1:
        cross = rZ*rB*m1/(ab*ab*aD0)
        st = ub(arb(ht)*(abs(N1)/aN0 + abs(Dz)/aD0) + RN/aN0 + rZ*rZ*m2/(2*ab*aD0) + cross/2)
        sr = ub(arb(hr)*abs(Db)/aD0 + rB*rB*m0/(ab**3*aD0) + cross/2)
        return dict(ok=False, G0=G0, st=st, sr=sr, k1=None)
    lz = N1/N0 - Dz/D0
    k1 = lz*EM
    k2 = -(Db/D0)*EP
    rem = abs(lz)*ez + abs(Db/D0)*eb + RN/aN0 + RD/aD0 + (Y*Y + X*Y)/(1 - Y)
    # splitting heuristic (not part of the proof): first-order + second-order contributions per direction
    cross = rZ*rB*m1/(ab*ab*aD0)
    st = ub(arb(ht)*(abs(N1)/aN0 + abs(Dz)/aD0) + RN/aN0 + rZ*rZ*m2/(2*ab*aD0) + cross/2)
    sr = ub(arb(hr)*abs(Db)/aD0 + rB*rB*m0/(ab**3*aD0) + cross/2)
    return dict(G0=G0, k1=k1, k2=k2, rem=arb(ub(rem)), ht=arb(ht), hr=arb(hr), st=st, sr=sr)

def halfplane_low(d, th):
    """rigorous lower bound of Im(e^{-i th pi} G) over the box (th in units of pi)."""
    e = (-I*PI*arb(th)).exp()
    A = e*d['G0']
    return lb(A.imag - abs((A*d['k1']).imag)*d['ht'] - abs((A*d['k2']).imag)*d['hr'] - abs(A)*d['rem'])

def halfplane_up(d, th):
    e = (-I*PI*arb(th)).exp()
    A = e*d['G0']
    return ub(A.imag + abs((A*d['k1']).imag)*d['ht'] + abs((A*d['k2']).imag)*d['hr'] + abs(A)*d['rem'])

def absG_low(d):
    # |G0 (1+u)| >= |G0| Re(1+u),  u = k1 dt + k2 drho + e,  |e| <= rem
    return lb(abs(d['G0'])*(1 - abs(d['k1'].real)*d['ht'] - abs(d['k2'].real)*d['hr'] - d['rem']))

def check_box(d, th1, th2, mmin):
    """returns (ok, margins) ; margins in units of |G0| (rough)"""
    if d is None or d.get('k1') is None: return False, None
    l1 = halfplane_low(d, th1); u2 = halfplane_up(d, th2)
    ok = l1 > 0 and u2 < 0
    m = None
    if mmin is not None:
        m = absG_low(d); ok = ok and m > mmin
    return ok, (l1, u2, m)
