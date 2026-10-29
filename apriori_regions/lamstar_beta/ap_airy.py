"""Criterio 'Airy exato' (AE) da regiao a priori, terna (1;1;0,5).

   Para cada componente nu, f'' = q_nu f com q_nu(t) = (Lambda_nu - i t)/nu LINEAR em t: as solucoes sao EXATAMENTE
   funcoes de Airy:  z(t) = c (t + i Lambda_nu),  c = nu^{-1/3} e^{-i pi/6}  (c^3 = -i/nu),
       A(t) = Ai(z(t))                (recessiva para t crescente),
       B(t) = Ai(e^{-2 pi i/3} z(t))  (recessiva para t decrescente),     W = A B' - A' B  (constante).
   Nao ha erro de Olver: as colunas descompostas R = A/A(0), X = B/B(S) e seus momentos
       m_R(f) = int_0^S e^{-beta t} f,   m_L(f) = int_0^S e^{-beta (S-t)} f
   sao calculados com Arb (acb.airy_ai, acb.integral: quadratura rigorosa para integrandos analiticos).
   O acoplamento (k_12, k_21 ~ sigma^{-2/3}) entra como no criterio de Schur (Prop. 2.6), com o nucleo de Green EXATO
   |G(s,t)| = |A(max) B(min)| / |W|  cotado por celula (Airy em bolas), em normas sup, na parede e em momentos ponderados.
   Linhas (pt_R, pt_L, m_R, m_L) (beta > 0).  D_m^{(0)} = matriz de momentos EXATA do problema desacoplado."""
import math, numpy as np
from flint import arb, acb, ctx
import ap_arb as AA
from ap_arb import up, lo, iv, A as Af, MARG, INF, COMPS, nu1, g, c, cm, k12, k21, I
ctx.prec = 96
W3 = acb(-0.5, -arb(3).sqrt()/2)                 # e^{-2 pi i/3}
def cfac(nu): return nu**(-arb(1)/3)*acb(arb(3).sqrt()/2, arb(-1)/2)     # nu^{-1/3} e^{-i pi/6}

def _moment(cc, Lnu, bet, S, side, which, rel=1e-20):
    """int_0^S w(t) f(t) dt, w = e^{-beta t} (side 'R') ou e^{-beta (S-t)} (side 'L'), f = A ('A') ou B ('B')"""
    def fA(x, an):
        z = cc*(x + I*Lnu)
        if which == 'B': z = W3*z
        if side == 'R': w = (-bet*x).exp()
        elif side == 'L': w = (-bet*(S - x)).exp()
        else: w = (S - x) if bet == 0 else (bet*(S - x)).sinh()/bet       # r_4' = e^{beta S} r_4 (linha reescalada)
        return w*z.airy_ai()
    return acb.integral(fA, acb(0), acb(S), rel_tol=rel, eval_limit=200000)

def comp_airy(Lnu, nu, S, bet, n):
    cc = cfac(nu)
    def Aat(t): return (cc*(t + I*Lnu)).airy_ai()
    def Bat(t): return (W3*cc*(t + I*Lnu)).airy_ai()
    z0 = cc*(I*Lnu)
    a0, ap0 = z0.airy_ai(), z0.airy_ai(derivative=1)
    b0, bp0 = (W3*z0).airy_ai(), (W3*z0).airy_ai(derivative=1)
    W = cc*(a0*W3*bp0 - ap0*b0)                     # Wronskiano em t
    A0 = a0; AS = Aat(S); B0 = b0; BS = Bat(S)
    Amax = np.empty(n); Bmax = np.empty(n)
    for k in range(n):
        t = S*iv(arb(k)/n, arb(k + 1)/n)
        Amax[k] = up(abs(Aat(t))); Bmax[k] = up(abs(Bat(t)))
    return dict(cc=cc, W=W, A0=A0, AS=AS, B0=B0, BS=BS, Amax=Amax, Bmax=Bmax)

def airy_eval(sig, bet, X, Y, n=512, detail=False, rel=1e-20):
    if not lo(bet) > 0: return dict(ok=False, why='beta=0')
    S = 2*sig**(arb(1)/3); s23 = sig**(-arb(2)/3)
    Lam = acb(X, Y)
    Dh = math.nextafter(up(S)/n, INF)*(1 + 1e-15)
    k = np.arange(n, dtype=float)
    Sl, Sh_ = lo(S), up(S); bl = max(lo(bet), 0.0)
    q = up((-(Af(bl)*Af(Sl))/n).exp()); qk = AA.cumprod_up(np.full(n, q))
    ebt = qk[:n]; eL = qk[n - 1::-1][:n]
    C = {}; main = {}; err = {}
    for key, nu, sh in COMPS:
        Lnu = Lam + nu*bet**2 + sh*s23
        d = comp_airy(Lnu, nu, S, bet, n)
        aW = lo(abs(d['W'])); aA0 = lo(abs(d['A0'])); aBS = lo(abs(d['BS']))
        if aW <= 0 or aA0 <= 0 or aBS <= 0: return dict(ok=False, why='W/A0/BS')
        Am, Bm = d['Amax'], d['Bmax']
        if not (np.all(np.isfinite(Am)) and np.all(np.isfinite(Bm))): return dict(ok=False, why='airy')
        # Green exato por celula: |G(s,t)| <= Amax[max] Bmax[min] / |W|
        cA = np.cumsum(Am[::-1])[::-1]                 # sum_{j>=i} Amax_j
        cB = np.concatenate([[0.0], np.cumsum(Bm)])[:-1]   # sum_{j<i} Bmax_j
        rows = (Bm*cA + Am*cB)*(1 + 1e-12)
        Ninf = float(np.max(rows))*Dh/aW*MARG
        envR = Am/aA0*MARG; envX = Bm/aBS*MARG
        RS = up(abs(d['AS']/d['A0'])); X0 = up(abs(d['B0']/d['BS']))
        # momentos exatos do problema desacoplado
        mRA = _moment(d['cc'], Lnu, bet, S, 'R', 'A', rel); mLA = _moment(d['cc'], Lnu, bet, S, 'L', 'A', rel)
        mRB = _moment(d['cc'], Lnu, bet, S, 'R', 'B', rel); mLB = _moment(d['cc'], Lnu, bet, S, 'L', 'B', rel)
        main[('mR', 'R', key)] = mRA/d['A0']; main[('mL', 'R', key)] = mLA/d['A0']
        main[('mR', 'X', key)] = mRB/d['BS']; main[('mL', 'X', key)] = mLB/d['BS']
        for kk in [('mR', 'R', key), ('mL', 'R', key), ('mR', 'X', key), ('mL', 'X', key)]: err[kk] = 0.0
        # Green ponderado e na parede (para cross bounds)
        g0 = Am*up(abs(d['B0']))/aW*MARG; gS = Bm*up(abs(d['AS']))/aW*MARG
        def weighted(w):
            a = np.cumsum((w*Bm))                           # sum_{s<=t} w_s Bmax_s  (A(t) e' o max)
            b = np.concatenate([np.cumsum((w*Am)[::-1])[::-1][1:], [0.0]])   # sum_{s>t} w_s Amax_s
            return (Am*a + Bm*b)*Dh/aW*MARG*(1 + 1e-12)
        C[key] = dict(nu=nu, Ninf=Ninf, envR=envR, envX=envX, Rinf=float(np.max(envR)), Xinf=float(np.max(envX)),
                      R1=float(np.sum(envR))*Dh*MARG, X1n=float(np.sum(envX))*Dh*MARG, RS=RS, X0=X0, Dh=Dh,
                      G0=float(np.sum(g0))*Dh*MARG, GS=float(np.sum(gS))*Dh*MARG, g0=g0, gS=gS,
                      WR=weighted(ebt), WL=weighted(eL))
    xover = cross_airy(C, s23)
    return AA.assemble(C, main, err, s23, 1.0, detail, rows=('mR', 'mL'), xover=xover)

def cross_airy(C, s23):
    """como ap_arb.cross_bounds, com o nucleo de Green exato (Airy) por celula"""
    k1 = up(k12*s23/nu1); kg = up(k21*s23/g); kap = {'1': k1, 'g': kg}
    tau = k1*kg*C['1']['Ninf']*C['g']['Ninf']*MARG
    if tau >= 1: return None
    out = {}
    for own in ('1', 'g'):
        o = 'g' if own == '1' else '1'; Co = C[o]; Cw = C[own]; Dh = Co['Dh']; ko = kap[o]; kw = kap[own]
        for fam in 'RX':
            f = Cw['envR'] if fam == 'R' else Cw['envX']
            finf = Cw['Rinf'] if fam == 'R' else Cw['Xinf']; f1 = Cw['R1'] if fam == 'R' else Cw['X1n']
            cinf = ko*Co['Ninf']*finf/(1 - tau)*MARG; c1n = ko*Co['Ninf']*f1/(1 - tau)*MARG
            out[(fam, own, 0)] = ko*(float(np.sum(Co['g0']*f))*Dh + Co['G0']*kw*Cw['Ninf']*cinf)*MARG
            out[(fam, own, 'S')] = ko*(float(np.sum(Co['gS']*f))*Dh + Co['GS']*kw*Cw['Ninf']*cinf)*MARG
            out[(fam, own, 'mR')] = ko*(float(np.sum(Co['WR']*f))*Dh + float(np.max(Co['WR']))*kw*Cw['Ninf']*c1n)*MARG
            out[(fam, own, 'mL')] = ko*(float(np.sum(Co['WL']*f))*Dh + float(np.max(Co['WL']))*kw*Cw['Ninf']*c1n)*MARG
            if 'W4' in Co:
                out[(fam, own, 'r4')] = ko*(float(np.sum(Co['W4']*f))*Dh + float(np.max(Co['W4']))*kw*Cw['Ninf']*c1n)*MARG
    return out
