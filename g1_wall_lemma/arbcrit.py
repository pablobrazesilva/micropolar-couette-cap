"""Arb (ball-arithmetic) port of the coupling criterion crit2.py (G1', Props. C5-C9), for BOXES
   sigma in [s1,s2], beta in [b1,b2] (b1 may be 0).  Returns a rigorous LOWER bound of
   margin = min(thK, 1 - leak)  (criterion holds on the whole box iff margin > 0).

   Weights (homogeneity: Mo(w/e2) = Mo(w)/e2):
     'm'   : w = e^{-beta s}                      (W_m)
     'one' : w = 1                                (L1 route)
     'd0'  : w = e^{-beta(S-s)}                   (lambda = 0, e2 = 1)
     'd1n' : w = 2 sinh(beta s)/(e^{beta S}-1)    (lambda = 1, NORMALISED by e2 = 1-e^{-beta S}; analytic at beta=0)
   Kernel-block integrals of the weights are closed forms with monotone primitives (arbprim);
   the t-integrals (double integral of route 1, h' integral of route 2) are rigorous upper sums over cells.
   Modular: crit_box(model, box, choice) -- a new criterion (crit3) can be plugged by writing another
   function with the same signature."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import math, json, os
from arbprim import *

HERE = os.path.dirname(os.path.abspath(__file__))
CP = _os.path.join(_CERT, 'g1_wall_lemma') + _os.sep
_cp = json.load(open(CP + 'cert_pi_0.05.json')); _cq = json.load(open(CP + 'cert_pipsi.json'))
# certified constants (upper bounds for Pi's, lower bounds for a, at, cos theta, m_F)
K = dict(Pi=arb(_cp['Pi']), Pi1=arb(_cp['Pi1']), Pipsi=arb(_cq['Pi_psi']), Pi1psi=arb(_cq['Pi1_psi']),
         a=arb('0.44452699'), at=arb('0.35346270'), ct=arb('0.545117'), mF=arb('0.6083'))
# P(S') upper bounds (G0 constantes, eta=0.05), P nonincreasing; values rounded UP
PT = [(2, '0.2494'), (6, '1.583e-2'), (10, '7.19e-5'), (11.5, '6.512e-6'), (14, '8.48e-8'),
      (18.2, '2.192e-11'), (20, '4.54e-13'), (23.6, '1.14e-16')]
_PTF = _os.path.join(_CERT, 'g0_constants', 'ptable_0.05.json')
if os.path.exists(_PTF):                       # fine certified table (G0-bis, ptable.py), step 0.05-0.25
    PT = PT + [(float(a), b) for a, b in json.load(open(_PTF))]
PT.sort()
def Pbound(Sp_lo):
    """upper bound of P(S') for all S' >= Sp_lo (float lower end); P nonincreasing."""
    if Sp_lo < 2: return arb('inf')
    v = arb('0.5051')*(-arb('0.3534')*arb(Sp_lo)).exp()
    for s, p in PT:
        if Sp_lo >= s: v = amin(v, arb(p))
    return hi(v)

def params(N, J, G):
    N, J, G = arb(N), arb(J), arb(G)
    g = G/J; kap = 4*N/J; nu1 = 1 + N; c = 2*N/(1 + N - g)
    return dict(N=N, J=J, G=G, g=g, kap=kap, nu1=nu1, c=c, k12=2*N*kap*(1 - g)/(1 + N - g)**2, k21=kap/2)

# ------------------------------------------------------------------ weights
class Weight:
    def __init__(s, kind, beta, S):
        s.k, s.b, s.S = kind, beta, S
        if kind == 'd1n': s.Kn = 2/(S*psi(beta*S))           # = 2 beta/(e^{beta S}-1)
    def val(s, T):
        b, S = s.b, s.S
        if s.k == 'm': return (-b*T).exp()
        if s.k == 'one': return ONE
        if s.k == 'd0': return (-b*pos(S - T)).exp()
        return s.Kn*T*shc(b*T)
    def wS(s):
        b, S = s.b, s.S
        return {'m': (-b*S).exp(), 'one': ONE, 'd0': ONE, 'd1n': 1 + (-b*S).exp()}[s.k]
    def w0(s): return {'m': ONE, 'one': ONE, 'd0': (-s.b*s.S).exp(), 'd1n': ZERO}[s.k]
    def wmax(s): return {'m': ONE, 'one': ONE, 'd0': ONE, 'd1n': 1 + (-s.b*s.S).exp()}[s.k]
    def G(s, T, r):     # int_T^S |w'(u)| e^{-r(u-T)} du
        b, S = s.b, s.S; L = pos(S - T)
        if s.k == 'm': return b*(-b*T).exp()*E(-(b + r), L)
        if s.k == 'one': return ZERO
        if s.k == 'd0': return b*(-r*L).exp()*E(r - b, L)
        return s.Kn*((b*T).exp()*E(b - r, L) + (-b*T).exp()*E(-b - r, L))/2
    def H(s, T, r):     # int_0^T |w'(u)| e^{-r(T-u)} du
        b, S = s.b, s.S
        if s.k == 'm': return b*(-r*T).exp()*E(r - b, T)
        if s.k == 'one': return ZERO
        if s.k == 'd0': return b*(-b*S - r*T).exp()*E(b + r, T)
        return s.Kn*(-r*T).exp()*(E(r + b, T) + E(r - b, T))/2
    def ThR(s, T, at):  # int_T^S w(u) e^{-at(u-T)} du
        b, S = s.b, s.S; L = pos(S - T)
        if s.k == 'm': return (-b*T).exp()*E(-(b + at), L)
        if s.k == 'one': return E(-at, L)
        if s.k == 'd0': return (-at*L).exp()*E(at - b, L)
        return s.Kn*SHE(b, -at, T, L)
    def ThX(s, T, at):  # int_0^T w(u) e^{-at(T-u)} du
        b, S = s.b, s.S
        if s.k == 'm': return (-at*T).exp()*E(at - b, T)
        if s.k == 'one': return E(-at, T)
        if s.k == 'd0': return (-b*pos(S - T)).exp()*E(-(b + at), T)
        return s.Kn*SHX(b, -at, T)

class WeightRefl:
    """lambda=0 far-wall weight w(s) = e^{-beta(S-s)}, evaluated in REFLECTED coordinate t' = S - t
       (avoids the S-correlation e^{beta dS}); it is the reflection of the 'm' weight."""
    refl = True
    def __init__(s, beta, S): s.b, s.S = beta, S; s.m = Weight('m', beta, S)
    def wmax(s): return ONE
    # values at t = S - Tp:
    def val(s, Tp): return (-s.b*Tp).exp()
    def GR(s, Tp, r):  return (-s.b*Tp).exp() + (-r*Tp).exp() + s.m.H(Tp, r)
    def GX(s, Tp, r):  return (-s.b*Tp).exp() + (-s.b*s.S - r*pos(s.S - Tp)).exp() + s.m.G(Tp, r)
    def ThR(s, Tp, at): return s.m.ThX(Tp, at)
    def ThX(s, Tp, at): return s.m.ThR(Tp, at)

# ------------------------------------------------------------------ t-cells
def cells(S, h0=0.02, grow=0.25, tmax_fine=None):
    """partition of [0, S_lo] (+ final cell [S_lo, S_hi]); width h0*(1+grow*t)."""
    Slo = float(lo(S).mid()) if False else float(S.lower())
    out = []; t = 0.0
    while t < Slo:
        t2 = min(Slo, t + h0*(1 + grow*t))
        out.append((t, t2)); t = t2
    if float(S.upper()) > Slo: out.append((Slo, float(S.upper())))
    return out

def int_exp_cell(r, t1, t2):   # int_{t1}^{t2} e^{-r t} dt, rigorous (t1,t2 floats)
    a = arb(t1); return (-r*a).exp()*E(-r, arb(t2) - a)

def route_integrals_refl(Wt, kb, rf, need2, cl):
    S = Wt.S; at = kb['at']; r = kb['r']; nu = kb['nu']
    k0 = kb['C']*K['Pi']; k1 = kb['C']*K['Pi1']*nu**(-ONE/3)
    P1s, Pps = K['Pi1psi'], K['Pipsi']
    dbl = ZERO; I2 = ZERO
    for (t1, t2) in cl:                       # cells in t' = S - t
        Tp = iv(t1, t2)
        ie = (-rf*pos(S - arb(t2))).exp()*E(-rf, arb(t2) - arb(t1))     # int_cell e^{-rf (S - t')} dt'
        ThR, ThX = Wt.ThR(Tp, at), Wt.ThX(Tp, at)
        dbl += ie*hi(ThR + ThX)
        if need2:
            hp = kb['C']*(amin(P1s*Wt.GX(Tp, r), nu**(-ONE/3)*K['Pi1']*ThX) + amin(P1s*Wt.GR(Tp, r), nu**(-ONE/3)*K['Pi1']*ThR))
            I2 += ie*hi(hp)
    R1 = k0*(Wt.ThR(S, at) + (-rf*S).exp()*Wt.ThX(ZERO, at)) + k1*dbl
    if not need2: return hi(R1), arb('inf')
    h0 = kb['C']*amin(nu**(ONE/3)*Pps*Wt.GR(S, r), K['Pi']*Wt.ThR(S, at))
    hS = kb['C']*amin(nu**(ONE/3)*Pps*Wt.GX(ZERO, r), K['Pi']*Wt.ThX(ZERO, at))
    return hi(R1), hi(h0 + (-rf*S).exp()*hS + I2)

def route_integrals(Wt, kb, rf, need2, cl):
    if getattr(Wt, 'refl', False): return route_integrals_refl(Wt, kb, rf, need2, cl)
    """returns (MU1, MU2) upper bounds for |int w y|/|phi_f|, y = Kernel_kb[R_f-hat], w = Wt.
       kb: kernel block dict (nu, r, at, C); rf: decay rate r of the forcing block."""
    S = Wt.S; at = kb['at']; r = kb['r']; nu = kb['nu']
    k0 = kb['C']*K['Pi']; k1 = kb['C']*K['Pi1']*nu**(-ONE/3)
    # ---- route 1:  k0 [int w e^{-at s} + e^{-rf S} int_0^S w(u) e^{-at(S-u)} du] + k1 dbl
    dbl = ZERO; I2 = ZERO
    P1s, Pps = K['Pi1psi'], K['Pipsi']
    for (t1, t2) in cl:
        T = iv(t1, t2)
        ie = int_exp_cell(rf, t1, t2)
        ThR, ThX = Wt.ThR(T, at), Wt.ThX(T, at)
        dbl += ie*hi(ThR + ThX)
        if need2:
            v = Wt.val(T)
            GR = v + Wt.wS()*(-r*pos(S - T)).exp() + Wt.G(T, r)
            GX = v + Wt.w0()*(-r*T).exp() + Wt.H(T, r)
            hp = kb['C']*(amin(P1s*GX, nu**(-ONE/3)*K['Pi1']*ThX) + amin(P1s*GR, nu**(-ONE/3)*K['Pi1']*ThR))
            I2 += ie*hi(hp)
    R1 = k0*(Wt.ThR(ZERO, at) + (-rf*S).exp()*Wt.ThX(S, at)) + k1*dbl
    if not need2: return hi(R1), arb('inf')
    h0 = kb['C']*amin(nu**(ONE/3)*Pps*(Wt.val(ZERO) + Wt.wS()*(-r*S).exp() + Wt.G(ZERO, r)), K['Pi']*Wt.ThR(ZERO, at))
    hS = kb['C']*amin(nu**(ONE/3)*Pps*(Wt.val(S) + Wt.w0()*(-r*S).exp() + Wt.H(S, r)), K['Pi']*Wt.ThX(S, at))
    R2 = h0 + (-rf*S).exp()*hS + I2
    return hi(R1), hi(R2)

# ------------------------------------------------------------------ the criterion (crit2, Props. C5-C9)
def crit_box(P, sig, beta, lam, th_rho, cone, route2=True, h0=0.02, grow=0.25, coupling=1):
    """P = params(...); sig, beta: arb intervals.  Returns dict(margin=float lower bound or -inf, why=...)."""
    a, at, ct, mF, Pi = K['a'], K['at'], K['ct'], K['mF'], K['Pi']
    g, nu1, c = P['g'], P['nu1'], P['c']
    k12, k21 = coupling*P['k12'], coupling*P['k21']
    S = 2*sig**(ONE/3); s23 = sig**(-TWO/3)
    th = ONE/3
    B1 = dict(nu=nu1, r=a*nu1**(-th), at=at*nu1**(-th), C=2*PI*nu1**(-2*th)*s23)
    Bg = dict(nu=g, r=a*g**(-th), at=at*g**(-th), C=2*PI*g**(-2*th)*s23)
    r1, rg, a1, ag = B1['r'], Bg['r'], B1['at'], Bg['at']
    k0 = {1: B1['C']*Pi, 'g': Bg['C']*Pi}; k1 = {1: B1['C']*K['Pi1']*nu1**(-th), 'g': Bg['C']*K['Pi1']*g**(-th)}
    rho = arb(th_rho)*a1
    I = lambda av: 1/(av - rho) + 1/(av + rho)
    N1, Ng = k0[1]*I(a1), k0['g']*I(ag)
    tau2 = k12*k21*N1*Ng
    out = dict(lam=lam, th_rho=th_rho, cone=cone)
    def fail(w): return dict(out, margin=-math.inf, why=w)
    if not tau2.upper() < 1: return fail('tau2')
    # ---- column U_1 (zeta)
    u1, u0 = k0['g']*I(ag), k0['g']/(ag + a1)
    B_1 = k21*u1/(1 - tau2)
    B20 = k12*k21*k0['g']/(ag + rho)*N1*B_1; A0 = k12*k0[1]/(a1 + rho)*B_1
    P1 = A0 + (2 - c)*(k21*u0 + B20)
    Ub = k0['g']*(1 + ((rho - r1)*S).exp()) + k1['g']*(1/(ag - rho) + 1/(ag + r1))
    Bp = k21*Ub/(1 - tau2); Ap = k12*N1*Bp
    Q1 = (-rho*S).exp()*(Ap + (2 - c)*Bp)
    K1n, Kgn = k0[1]*2/a1, k0['g']*2/ag
    tau1 = k12*k21*K1n*Kgn
    if not tau1.upper() < 1: return fail('tau1')
    cl = cells(S, h0, grow)
    W1 = Weight('one', beta, S)
    L1u = route_integrals(W1, Bg, r1, False, cl)[0]
    L1B = k21*L1u/(1 - tau1); L1A = k12*K1n*L1B
    Wm = Weight('m', beta, S)
    Wd = Weight('d1n', beta, S) if lam == 1 else WeightRefl(beta, S)
    def wrho(Wt):   # int_0^S w e^{-rho s}
        return Wt.ThR(S, rho) if getattr(Wt, 'refl', False) else Wt.ThR(ZERO, rho)
    Mo1 = {}
    for nm, Wt in (('m', Wm), ('d', Wd)):
        R1, R2 = route_integrals(Wt, Bg, r1, route2, cl)
        MU = amin(R1, R2)
        hi_sup = (Ap + c*tau2*Bp)*wrho(Wt)
        hi_L1 = Wt.wmax()*(L1A + c*tau1*L1B)
        Mo1[nm] = hi(c*k21*MU + amin(hi_sup, hi_L1))
    # ---- column U_g (w)
    v1, v0 = k0[1]*I(a1), k0[1]/(a1 + ag)
    A_g = k12*v1/(1 - tau2)
    A20 = k12*k21*k0[1]/(a1 + rho)*Ng*A_g; B0 = k21*k0['g']/(ag + rho)*A_g
    Pg = k12*v0 + A20 + (2 - c)*B0
    Vb = k0[1]*(1 + ((rho - rg)*S).exp()) + k1[1]*(1/(a1 - rho) + 1/(a1 + rg))
    Agp = k12*Vb/(1 - tau2); Bgp = k21*Ng*Agp
    Qg = (-rho*S).exp()*(Agp + (2 - c)*Bgp)
    L1v = route_integrals(W1, B1, rg, False, cl)[0]
    L1Ag = k12*L1v/(1 - tau1); L1Bg = k21*Kgn*L1Ag
    Mog = {}
    for nm, Wt in (('m', Wm), ('d', Wd)):
        R1, R2 = route_integrals(Wt, B1, rg, route2, cl)
        MV = amin(R1, R2)
        hi_sup = (tau2*Agp + c*Bgp)*wrho(Wt)
        hi_L1 = Wt.wmax()*(tau1*L1Ag + c*L1Bg)
        Mog[nm] = hi(k12*MV + amin(hi_sup, hi_L1))
    if not P1.upper() < 1: return fail('P1')
    # ---- uncoupled pieces
    b = beta
    def Bm(r): return 1 + (-(b + r)*S).exp() + b/(b + r)*(1 - (-(b + r)*S).exp())
    def Dl_e2(r):   # e^{-2al} Delta_nu(lam) / e2
        if lam == 1: return (1 + (-b*S).exp())*(-r*S).exp() + Wd.Kn*(E(b - r, S) + E(-b - r, S))/2
        t1 = b*(-b*S).exp()*E(b - r, S); t2 = b*(-r*S).exp()*E(r - b, S)   # same quantity, two forms
        return (-b*S).exp() + (-r*S).exp() + (t1 if t1.upper() < t2.upper() else t2)
    epx = ((2 - c)*P1 + Pg)/((1 - P1)*(2 - c))
    if not epx.upper() < 1: return fail('epx')
    ph1 = nu1**th/mF
    xb = (2 - c)*(1 + epx)
    kmg = ct*rg/(b + rg) - 2*(-rg*S).exp()
    D1, Dg = Dl_e2(r1), Dl_e2(rg)
    if cone:
        thp = ct.acos() + epx.asin()
        if not thp.upper() < (PI/2).lower(): return fail('cone')
        km1 = (1 - epx)*thp.cos()*r1/(b + r1) - 2*(1 + epx)*(-r1*S).exp()
        if not (km1.lower() > 0 and kmg.lower() > 0): return fail('km')
        th_max = 1 - amax((1 + epx)*D1/km1, Dg/kmg)
        eW = amax(xb*Mo1['m']/((2 - c)*km1), Mog['m']/(c*kmg))
        eD = amax(xb*Mo1['d']/((2 - c)*km1), Mog['d']/(c*kmg))
    else:
        km1 = ct*r1/(b + r1) - 2*(-r1*S).exp()
        if not (km1.lower() > 0 and kmg.lower() > 0): return fail('km')
        th_max = 1 - amax(D1/km1, Dg/kmg)
        dx = (2 - c)*epx
        eW = amax((dx*Bm(r1) + xb*Mo1['m'])/((2 - c)*km1), Mog['m']/(c*kmg))
        eD = amax((dx*D1 + xb*Mo1['d'])/((2 - c)*km1), Mog['d']/(c*kmg))
    thK = th_max - eW - eD
    out.update(P1=P1, Pg=Pg, eW=eW, eD=eD, th_max=th_max, thK=thK)
    if not (thK.lower() > 0 and eW.upper() < 1): return fail('thK')
    # ---- leakage (Prop. C9)
    Sp1 = float((S*nu1**(-th)).lower()); Spg = float((S*g**(-th)).lower())
    P_1, P_g = Pbound(Sp1), Pbound(Spg)
    a4 = xb*(nu1**(-th)*P_1 + Q1)/(2 - c); b4 = ((2 - c)*g**(-th)*P_g + Qg)/c
    E4 = amax(a4/km1, b4/kmg)
    e3 = ph1*(nu1**(-th)*P_1 + Q1)
    lamA = arb(lam)
    y_m = Bm(r1) + Mo1['m']; y_l = D1 + Mo1['d']
    dR = (1 - th_max + eD)/(1 - eW)
    Y = y_m*dR + y_l; Z2 = 2*y_m + lamA*y_l + lamA*y_m*dR + dR*y_l; Z3 = 2*y_m**2 + 2*lamA*y_m*y_l + y_l**2
    q = 1 - P1
    leak = 2*E4*ph1*Y/(q*thK) + e3**2/q**2 + 2*E4*e3*ph1*Z2/(q**2*thK) + E4**2*ph1**2*Z3/(q**2*thK*(1 - eW))
    m = amin(thK, 1 - leak)
    out.update(E4=E4, leak=leak)
    return dict(out, margin=float(m.lower()), why='ok' if m.lower() > 0 else 'leak')
