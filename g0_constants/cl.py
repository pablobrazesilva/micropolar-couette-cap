"""G0/consts: rigorous (Arb ball) enclosures on the line z_eta(t) = t e^{-i pi/6} - eta e^{i pi/3}, b = 0, of
   F(z)     = Ai(z)/int_z^inf Ai            (= G(z,0) of g0lib)
   u(z)     = Re(e^{-i pi/6} F(z))          (-> a(eta))
   ut(z)    = Re(-e^{-i pi/6} Ai'(z)/Ai(z)) (-> at(eta))
Representation (Lemma 1 of the G0 report):  Ai'/Ai = -w0 + N'(z)/N(z),  F = N/D.
"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, math
sys.path.insert(0, _os.path.join(_R, 'g0_sector_lemma'))
from g0lib import *

def moments4(a, g):
    """upper bounds m_k = int |s|^k exp(-a s^2 + g|s|) ds, k=0..3 (recursion 2a I_{k+1} - g I_k = k I_{k-1})."""
    a = arb(a); g = arb(g)
    I0 = (PI/a).sqrt()/2*(g*g/(4*a)).exp()*(1 + (g/(2*a.sqrt())).erf())
    I1 = (1 + g*I0)/(2*a)
    I2 = (I0 + g*I1)/(2*a)
    I3 = (2*I1 + g*I2)/(2*a)
    return 2*I0, 2*I1, 2*I2, 2*I3

def centre_N012D(zc, w0):
    """N0, N1=N', N2=N'' at zc and D0 = D(zc,0), each a ball including the |s|>S truncation error."""
    a = float(w0.real.mid())
    q = I*(w0*w0 - zc)
    g0 = ub(abs((w0*w0 - zc).imag))
    S = (g0 + math.sqrt(g0*g0 + 4*a*46))/(2*a); S = math.ceil(S*8)/8
    ph = lambda x: (-(w0*x*x) - I*x*x*x/3 + q*x).exp()
    N0 = acb.integral(lambda x, an: ph(x), -S, S)
    N1 = acb.integral(lambda x, an: -I*x*ph(x), -S, S)
    N2 = acb.integral(lambda x, an: -x*x*ph(x), -S, S)
    D0 = acb.integral(lambda x, an: ph(x)/(w0 + I*x), -S, S)
    aa = arb(a); Sa = arb(S); gg = arb(g0)
    E = (-(aa*Sa*Sa) + gg*Sa).exp()
    assert lb(2*aa*Sa - gg - 2/Sa) > 0
    e0 = 2*E/(2*aa*Sa - gg)
    e1 = 2*Sa*E/(2*aa*Sa - gg - 1/Sa)
    e2 = 2*Sa*Sa*E/(2*aa*Sa - gg - 2/Sa)
    return disk(N0, e0), disk(N1, e1), disk(N2, e2), disk(D0, e0/aa)

def box_line(t1, t2, eta, amin=0.6):
    """Rigorous bounds over t in [t1,t2] on the line, b=0.
    returns dict(Fup, Flo, ulo, utlo) (floats, rounded outward) or None."""
    tc = 0.5*(t1 + t2); ht = max(ub(arb(t2) - arb(tc)), ub(arb(tc) - arb(t1)))
    zt = line_point(tc, eta); zc_c = complex(float(zt.real.mid()), float(zt.imag.mid()))
    zc = acb(zc_c.real, zc_c.imag); ez = arb(ub(abs(zt - zc)))
    rZ = arb(ht) + ez
    w0 = choose_w0(zc_c, amin); a = w0.real
    gam = arb(ub(abs((w0*w0 - zc).imag))) + rZ
    m0, m1, m2, m3 = moments4(a, gam)
    N0, N1, N2, D0 = centre_N012D(zc, w0)
    aN0 = abs(N0); aD0 = abs(D0)
    if lb(aN0) <= 0 or lb(aD0) <= 0: return None
    out = {}
    # ---- F = N/D : first-order model with second-order remainder (Lemma C with r_B = 0)
    Dz = w0*D0 - N0
    RN = rZ*rZ*m2/2; RD = rZ*rZ*m2/(2*a)
    X = (abs(N1)*rZ + RN)/aN0
    Y = (abs(Dz)*rZ + RD)/aD0
    if ub(Y) >= 1 or ub(X) >= 1: return None
    G0 = N0/D0
    lz = N1/N0 - Dz/D0
    k1 = lz*EM
    rem = abs(lz)*ez + RN/aN0 + RD/aD0 + (Y*Y + X*Y)/(1 - Y)
    A = EM*G0
    out['Fup'] = ub(abs(G0)*(1 + abs(k1)*arb(ht) + rem))
    out['Flo'] = lb(abs(G0)*(1 - abs(k1)*arb(ht) - rem))
    out['ulo'] = lb(A.real - abs((A*k1).real)*arb(ht) - abs(A)*rem)
    # ---- f = N'/N ;  -Ai'/Ai = w0 - f
    f0 = N1/N0
    k = N2/N0 - f0*f0
    Yn = (abs(N1)*rZ + RN)/aN0
    if ub(Yn) >= 1: return None
    rho = (rZ*rZ*m3/2)/aN0 + abs(f0)*RN/aN0
    Ef = (rho + abs(k)*rZ*Yn)/(1 - Yn)
    B = EM*(w0 - f0)
    out['utlo'] = lb(B.real - abs((EM*EM*k).real)*arb(ht) - abs(k)*ez - Ef)
    out['F0'] = G0; out['ut0'] = B.real; out['u0'] = A.real
    return out

def point(t, eta, prec=128):
    """rigorous point values (F, u, ut) at a single t (decimal string ok)."""
    old = ctx.prec; ctx.prec = prec
    z = acb(arb(t))*EM - arb(eta)*E3
    zm = complex(float(z.real.mid()), float(z.imag.mid())); w0 = choose_w0(zm, 0.6)
    N0, N1, N2, D0 = centre_N012D(z, w0)
    F = N0/D0; u = (EM*F).real; ut = (EM*(w0 - N1/N0)).real
    ctx.prec = old
    return F, u, ut

# ---------------- tails (|t| >= T), b = 0 ----------------
def eps_D(r, kappa):
    """|N'/N| <= kappa^{-5/2} / (4 r^2 (1 - (5/48) kappa^{-7/2} r^{-3}))  (w0 = sqrt z)."""
    r = arb(r); k = arb(kappa)
    den = 1 - arb(5)/48*k**(-arb(7)/2)/r**3
    assert lb(den) > 0
    return k**(-arb(5)/2)/(4*r*r*den)

def tail_side(T, eta, side):
    """for |t| >= T on the given side: rigorous (ulo, utlo, Fup_coef) where |F| <= Fup_coef*(t^2+eta^2)^{1/4}."""
    T = arb(T); eta = arb(eta)
    al = (eta/T).atan()/2
    r = T.sqrt()
    k = (PI/12 + al).cos() if side == 'right' else (5*PI/12 + al).cos()
    eB = tailB_eps(r, k); eD = eps_D(r, k)
    assert ub(eB) < 1
    ang = PI/4 + al + eB.asin()
    assert ub(ang) < ub(PI/2)
    ulo = r/(1 + eB)*ang.cos()
    utlo = r*(PI/4 + al).cos() - eD
    return lb(ulo), lb(utlo), ub(1/(1 - eB)), ub(eB), ub(eD)
