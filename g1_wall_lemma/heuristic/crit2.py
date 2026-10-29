"""Criterio de acoplamento (Teorema C), versao uniforme melhorada (A'): ponto flutuante.
   Pesos 'd' escalados por e^{-2alpha} (sem overflow).  Ver relatorio, secoes 3-5."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import numpy as np, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
_cp = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pi_0.05.json'))); _cq = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pipsi.json')))
CONST = dict(Pi=_cp['Pi'], Pi1=_cp['Pi1'], Pipsi=_cq['Pi_psi'], Pi1psi=_cq['Pi1_psi'],
             a=0.44452699, at=0.35346270, ct=0.545117, mF=0.6083)
PT = [(2,0.2494),(6,1.583e-2),(10,7.19e-5),(11.5,6.512e-6),(14,8.48e-8),(18.2,2.192e-11),(20,4.54e-13),(23.6,1.14e-16)]
def Pbound(Sp):
    if Sp < 2: return np.inf
    v = 0.5051*np.exp(-0.3534*Sp)
    for s, p in PT:
        if Sp >= s: v = min(v, p)
    return v
def E(p, L):                       # int_0^L e^{p u} du  (vectorized in L)
    L = np.asarray(L, float)
    if abs(p) < 1e-14: return L.copy() if L.ndim else float(L)
    return np.expm1(p*L)/p
def X(p, q, L):                    # (e^{pL}-e^{qL})/(p-q) = e^{qL} E(p-q, L)
    return np.exp(q*np.asarray(L, float))*E(p - q, L)
def params(N, J, G):
    g = G/J; kap = 4*N/J; nu1 = 1 + N; c = 2*N/(1 + N - g)
    return dict(g=g, kap=kap, nu1=nu1, c=c, k12=2*N*kap*(1 - g)/(1 + N - g)**2, k21=kap/2)
def EL(lc, q, L):
    """e^{lc} * int_0^L e^{q u} du, estavel (lc, L arrays ou escalares; q escalar)"""
    lc = np.asarray(lc, float); L = np.asarray(L, float)
    if abs(q)*max(np.max(L), 1e-300) < 1e-9: return np.exp(lc)*L
    if q > 0: return (np.exp(lc + q*L) - np.exp(lc))/q
    return np.exp(lc)*np.expm1(q*L)/q
class Weight:
    """w(s) on [0,S]: mode 'm': e^{-beta s}; mode 'd': e^{-beta(S-s)} - lam e^{-beta(S+s)} (= e^{-2al}(e^{beta s}-lam e^{-beta s}))
       componentes (lc, sgn, k): w = sum sgn e^{lc + k s}"""
    def __init__(s, mode, beta, lam, S): s.m, s.b, s.l, s.S = mode, beta, lam, S
    def comps(s):
        b, S = s.b, s.S
        if s.m == 'm': return [(0.0, 1.0, -b)]
        c = [(-b*S, 1.0, b)]
        if s.l != 0: c.append((-b*S + np.log(s.l), -1.0, -b))
        return c
    def w(s, t): return sum(sg*np.exp(lc + k*np.asarray(t, float)) for lc, sg, k in s.comps())
    def wmax(s): return 1.0 if s.m == 'm' else float(s.w(s.S))       # w monotona: max em 0 (m) ou S (d)
    def wint(s, p, shift=0.0): return sum(sg*EL(lc + shift, k + p, s.S) for lc, sg, k in s.comps())   # e^{shift} int w e^{p s}
    def Gam_R(s, t, r):
        S = s.S; out = np.abs(s.w(t)) + abs(s.w(S))*np.exp(-r*(S - t))
        for lc, sg, k in s.comps(): out = out + abs(k)*EL(lc + k*t, k - r, S - t)
        return out
    def Gam_X(s, t, r):
        out = np.abs(s.w(t)) + abs(s.w(0.0))*np.exp(-r*t)
        for lc, sg, k in s.comps(): out = out + abs(k)*EL(lc - r*t, k + r, t)
        return out
    def The_X(s, t, at): return sum(sg*EL(lc - at*t, k + at, t) for lc, sg, k in s.comps())
    def The_R(s, t, at): return sum(sg*EL(lc + k*t, k - at, s.S - t) for lc, sg, k in s.comps())
    def dbl(s, r, at):
        S = s.S; tt = np.linspace(0, S, 6001); tot = 0.0
        for lc, sg, k in s.comps():
            f = EL(lc + (k - r)*tt, k - at, S - tt) + EL(lc - (r + at)*tt, k + at, tt)
            tot += sg*np.trapezoid(f, tt)
        return tot
def col_bounds(kb, fb, Wt, K):
    """first-order response y = Kernel_kb [ r_f ] of a forcing column with block fb.
       kb, fb: dicts with nu, r (psi rate), at (Ai rate), C (=2 pi nu^{-2/3} sigma^{-2/3}).
       returns (route1, route2) bounds of |int w y| / |phi_f|  (w = weight Wt)"""
    S = Wt.S; Pi, Pi1, Pps, P1ps = K['Pi'], K['Pi1'], K['Pipsi'], K['Pi1psi']
    k0 = kb['C']*Pi; k1 = kb['C']*Pi1*kb['nu']**(-1/3)
    # route 1: pointwise kernel IBP then integrate against w
    r1 = k0*(Wt.wint(-kb['at']) + np.exp(-fb['r']*S)*Wt.wint(kb['at'], -kb['at']*S)) + k1*Wt.dbl(fb['r'], kb['at'])
    # route 2: dual function h = Kernel^T w, IBP against Psi_f
    t = np.linspace(0, S, 6001); nu, r, at = kb['nu'], kb['r'], kb['at']
    GR, GX, TR, TX = Wt.Gam_R(t, r), Wt.Gam_X(t, r), Wt.The_R(t, at), Wt.The_X(t, at)
    h0 = kb['C']*min(nu**(1/3)*Pps*GR[0], Pi*TR[0]); hS = kb['C']*min(nu**(1/3)*Pps*GX[-1], Pi*TX[-1])
    hp = kb['C']*(np.minimum(P1ps*GX, nu**(-1/3)*Pi1*TX) + np.minimum(P1ps*GR, nu**(-1/3)*Pi1*TR))
    r2 = h0 + np.exp(-fb['r']*S)*hS + np.trapezoid(np.exp(-fb['r']*t)*hp, t)*1.001
    return r1, r2
def crit(N, J, G, sig, beta, lam, th_rho=0.5, K=CONST, coupling=1.0, route='best', cone=True):
    P = params(N, J, G); g, nu1, c = P['g'], P['nu1'], P['c']; k12, k21 = coupling*P['k12'], coupling*P['k21']
    a, at, ct, mF, Pi = K['a'], K['at'], K['ct'], K['mF'], K['Pi']
    S = 2*sig**(1/3); al = beta*S/2; s23 = sig**(-2/3)
    B1 = dict(nu=nu1, r=a*nu1**(-1/3), at=at*nu1**(-1/3), C=2*np.pi*nu1**(-2/3)*s23)
    Bg = dict(nu=g, r=a*g**(-1/3), at=at*g**(-1/3), C=2*np.pi*g**(-2/3)*s23)
    r1, rg, a1, ag = B1['r'], Bg['r'], B1['at'], Bg['at']
    k0 = {1: B1['C']*Pi, 'g': Bg['C']*Pi}; k1 = {1: B1['C']*K['Pi1']*nu1**(-1/3), 'g': Bg['C']*K['Pi1']*g**(-1/3)}
    rho = th_rho*a1
    I = lambda av: 1/(av - rho) + 1/(av + rho)
    N1, Ng = k0[1]*I(a1), k0['g']*I(ag)
    tau2 = k12*k21*N1*Ng
    out = dict(S=S, al=al, tau2=tau2, lam=lam, th_rho=th_rho)
    if tau2 >= 1: return dict(out, ok=False, why='tau2')
    Wm, Wd = Weight('m', beta, lam, S), Weight('d', beta, lam, S)
    e2 = -np.expm1(-beta*S) if lam == 1 else 1 - lam*np.exp(-beta*S)     # e^{-2al}(e^{2al}-lam)
    # ---------- coluna U_1 (zeta): units 1
    u1, u0 = k0['g']*I(ag), k0['g']/(ag + a1)
    B_1 = k21*u1/(1 - tau2)
    B20 = k12*k21*k0['g']/(ag + rho)*N1*B_1; A0 = k12*k0[1]/(a1 + rho)*B_1
    P1 = A0 + (2 - c)*(k21*u0 + B20)
    Ub = k0['g']*(1 + np.exp((rho - r1)*S)) + k1['g']*(1/(ag - rho) + 1/(ag + r1))
    Bp = k21*Ub/(1 - tau2); Ap = k12*N1*Bp
    Q1 = np.exp(-rho*S)*(Ap + (2 - c)*Bp)
    # rota L1 para termos de ordem superior: ||K_nu||_{L1->L1} <= k0_nu 2/at_nu
    K1n, Kgn = k0[1]*2/a1, k0['g']*2/ag
    tau1 = k12*k21*K1n*Kgn
    if tau1 >= 1: return dict(out, ok=False, why='tau1')
    W1 = Weight('m', 0.0, lam, S)
    L1u = col_bounds(Bg, B1, W1, K)[0]           # int |u| / |phi_1|  (rota 1, w = 1)
    L1B = k21*L1u/(1 - tau1); L1A = k12*K1n*L1B
    Mo1 = {}
    for nm, Wt in (('m', Wm), ('d', Wd)):
        R1, R2 = col_bounds(Bg, B1, Wt, K)
        MU = {'1': R1, '2': R2, 'best': min(R1, R2)}[route]
        hi_sup = (Ap + c*k12*k21*Ng*N1*Bp)*Wt.wint(-rho)
        hi_L1 = Wt.wmax()*(L1A + c*tau1*L1B)
        Mo1[nm] = c*k21*MU + min(hi_sup, hi_L1)
    # ---------- coluna U_g (w)
    v1, v0 = k0[1]*I(a1), k0[1]/(a1 + ag)
    A_g = k12*v1/(1 - tau2)
    A20 = k12*k21*k0[1]/(a1 + rho)*Ng*A_g; B0 = k21*k0['g']/(ag + rho)*A_g
    Pg = k12*v0 + A20 + (2 - c)*B0
    Vb = k0[1]*(1 + np.exp((rho - rg)*S)) + k1[1]*(1/(a1 - rho) + 1/(a1 + rg))
    Agp = k12*Vb/(1 - tau2); Bgp = k21*Ng*Agp
    Qg = np.exp(-rho*S)*(Agp + (2 - c)*Bgp)
    L1v = col_bounds(B1, Bg, W1, K)[0]           # int |v| / |phi_g|
    L1Ag = k12*L1v/(1 - tau1); L1Bg = k21*Kgn*L1Ag
    Mog = {}
    for nm, Wt in (('m', Wm), ('d', Wd)):
        R1, R2 = col_bounds(B1, Bg, Wt, K)
        MV = {'1': R1, '2': R2, 'best': min(R1, R2)}[route]
        hi_sup = (tau2*Agp + c*Bgp)*Wt.wint(-rho)
        hi_L1 = Wt.wmax()*(tau1*L1Ag + c*L1Bg)
        Mog[nm] = k12*MV + min(hi_sup, hi_L1)
    if P1 >= 1: return dict(out, ok=False, why='P1')
    # ---------- desacoplado (pesos d escalados por e^{-2al})
    def Bm(r): return 1 + np.exp(-(beta + r)*S) + beta*(-np.expm1(-(beta + r)*S))/(beta + r)
    def Dl(r): # e^{-2al} * [ (1-lam) + (e^{bS}-lam e^{-bS}) e^{-rS} + beta int (e^{bs}+lam e^{-bs}) e^{-rs} ]
        return (1 - lam)*np.exp(-beta*S) + (1 - lam*np.exp(-2*beta*S))*np.exp(-r*S) + beta*(EL(-beta*S, beta - r, S) + lam*EL(-beta*S, -beta - r, S))
    # x = (2-c)(1+eps_x)
    epx = ((2 - c)*P1 + Pg)/((1 - P1)*(2 - c))
    if epx >= 1: return dict(out, ok=False, why='epx')
    ph1 = nu1**(1/3)/mF
    xb = (2 - c)*(1 + epx)
    kmg = ct*rg/(beta + rg) - 2*np.exp(-rg*S)
    if cone:   # variante (i): perturbacao pontual absorvida no cone (rotacao <= arcsin eps_x)
        thp = np.arccos(ct) + np.arcsin(epx)
        if thp >= np.pi/2: return dict(out, ok=False, why='cone')
        km1 = (1 - epx)*np.cos(thp)*r1/(beta + r1) - 2*(1 + epx)*np.exp(-r1*S)
        if km1 <= 0 or kmg <= 0: return dict(out, ok=False, why='km')
        th_max = 1 - max((1 + epx)*Dl(r1)/(e2*km1), Dl(rg)/(e2*kmg))
        eW = max(xb*Mo1['m']/((2 - c)*km1), Mog['m']/(c*kmg))
        eD = max(xb*Mo1['d']/((2 - c)*e2*km1), Mog['d']/(c*e2*kmg))
    else:      # variante (ii): termo (x-(2-c)) M_1 tratado como perturbacao
        km1 = ct*r1/(beta + r1) - 2*np.exp(-r1*S)
        if km1 <= 0 or kmg <= 0: return dict(out, ok=False, why='km')
        th_max = 1 - max(Dl(r1)/(e2*km1), Dl(rg)/(e2*kmg))
        dx = (2 - c)*epx
        eW = max((dx*Bm(r1) + xb*Mo1['m'])/((2 - c)*km1), Mog['m']/(c*kmg))
        eD = max((dx*Dl(r1) + xb*Mo1['d'])/((2 - c)*e2*km1), Mog['d']/(c*e2*kmg))
    thK = th_max - eW - eD
    out.update(P1=P1, Pg=Pg, eW=eW, eD=eD, th_max=th_max, thK=thK, Mo1=Mo1['m'], Mog=Mog['m'])
    if thK <= 0 or eW >= 1: return dict(out, ok=False, why='thK')
    Sp1, Spg = 2*(sig/nu1)**(1/3), 2*(sig/g)**(1/3)
    P_1, P_g = Pbound(Sp1), Pbound(Spg)
    a4 = xb*(nu1**(-1/3)*P_1 + Q1)/(2 - c); b4 = ((2 - c)*g**(-1/3)*P_g + Qg)/c
    E4 = max(a4/km1, b4/kmg)
    e3 = ph1*(nu1**(-1/3)*P_1 + Q1)
    y_m = Bm(r1) + Mo1['m']; y_l = (Dl(r1) + Mo1['d'])/e2
    dR = (1 - th_max + eD)/(1 - eW)
    Y = y_m*dR + y_l; Z2 = 2*y_m + lam*y_l + lam*y_m*dR + dR*y_l; Z3 = 2*y_m**2 + 2*lam*y_m*y_l + y_l**2
    q = 1 - P1
    leak = 2*E4*ph1*Y/(q*thK) + e3**2/q**2 + 2*E4*e3*ph1*Z2/(q**2*thK) + E4**2*ph1**2*Z3/(q**2*thK*(1 - eW))
    out.update(E4=E4, e3=e3, leak=leak, ok=bool(leak < 1))
    return out
def best(N, J, G, sig, beta, lams=(0.0, 1.0), ths=(0.35, 0.5, 0.65, 0.8), cones=(True, False), **kw):
    bv = None
    for lam in lams:
        for th in ths:
          for cn in cones:
            o = crit(N, J, G, sig, beta, lam, th, cone=cn, **kw); o['cone'] = cn
            sc = (1 - o['leak'])*o['thK'] if o.get('ok') else -1 - (o.get('thK', -9) if np.isfinite(o.get('thK', -9)) else -9)*0
            if bv is None or sc > bv[0]: bv = (sc, o)
    return bv
if __name__ == '__main__':
    import sys
    N, J, G = map(float, sys.argv[1:4])
    for sig in [1e3, 2e3, 3e3, 5e3, 1e4, 3e4]:
        for beta in [1e-6, 0.02, 0.1, 0.3, 1.0, 3.0]:
            sc, o = best(N, J, G, sig, beta)
            print('sig=%7.0f beta=%7.1e ok=%-5s lam=%g th=%.2f P1=%.3f eW=%.3f eD=%.3f thmax=%.3f thK=%.3f leak=%.1e %s' % (
                sig, beta, o.get('ok'), o['lam'], o['th_rho'], o.get('P1', np.nan), o.get('eW', np.nan), o.get('eD', np.nan),
                o.get('th_max', np.nan), o.get('thK', np.nan), o.get('leak', np.nan), o.get('why', '')))
