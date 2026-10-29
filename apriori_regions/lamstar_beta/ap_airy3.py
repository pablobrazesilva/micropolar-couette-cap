"""Criterio 'Airy exato' (AE) em CAIXAS, versao com DERIVADA LOGARITMICA (Riccati) -- terna (1;1;0,5).

   Para cada componente nu:  f'' = q f,  q = (Lambda_nu - i t)/nu (linear em t)  =>  solucoes EXATAS
       A(t) = Ai(z),  B(t) = Ai(w),  z = c (t + i Lambda_nu),  w = e^{-2 pi i/3} z,  c = nu^{-1/3} e^{-i pi/6},
       W = A B' - A' B = nu^{-1/3}/(2 pi)   (DLMF 9.2.8).
   Colunas R = A/A(0), X = B/B(S).  Nada e' aproximado (sem Olver).  Cotas na caixa via a derivada logaritmica
       Y(u) = Ai'(u)/Ai(u),   Y' = u - Y^2  (Riccati),
   cuja serie de Taylor em u0 tem majorante b_k (b_1 = |u0 - y0^2|, b_{k+1} = (d_{k1} + 2|y0| b_k + sum b_j b_{k-j})/(k+1))
   com cauda geometrica rigorosa (Lema R no relatorio).  Assim, na celula t x caixa de parametros (disco de raio rho):
       |Y(u) - Y(u0)| <= dY(u0, rho)   (dY pequeno: ~ |u0 - y0^2| rho),
   e  |R(t)| = exp(-int_0^t Re(A'/A)),  |X(t)| = exp(-int_t^S Re(B'/B)),  |A B| = |P| com (log P)' = Y(z) + w/z Y(w)
   sao cotados por somas por celula -- correlacionados (sem 'wrapping' da normalizacao).
   Termos principais (momentos exatos do problema desacoplado): centro da caixa (acb.integral) + valor medio.
   Linhas (pt_R, pt_L, m_R, m_L), beta > 0."""
import math, numpy as np
from flint import arb, acb, ctx
import ap_arb as AA
from ap_arb import up, lo, iv, A as Af, MARG, INF, COMPS, nu1, g, c, cm, k12, k21, I
import ap_airy as AY
ctx.prec = 96
W3 = AY.W3
KR = 40

def ydisc(u0s, rho):
    """u0s: lista de pontos acb; rho: array de raios.  Devolve (y0 lista acb, |y0| up, dY up) com
       |Y(u) - y0| <= dY para |u - u0| <= rho (e Y analitica no disco)."""
    n = len(u0s); y0 = []; ay = np.empty(n); b1 = np.empty(n); ai = []
    for i, u in enumerate(u0s):
        a0 = u.airy_ai(); a1 = u.airy_ai(derivative=1)
        if not (lo(abs(a0)) > 0): raise ValueError('Ai(u0) ~ 0')
        y = a1/a0; y0.append(y); ai.append(a0)
        ay[i] = up(abs(y)); b1[i] = up(abs(u - y*y))
    UP = 1 + 1e-15
    rho = np.asarray(rho, float)
    b = [np.zeros(n), b1]                               # b_0 (nao usado), b_1
    for k in range(1, KR):                              # b_{k+1}
        s = 2*ay*b[k]
        for j in range(1, k):
            s = s + b[j]*b[k - j]
        if k == 1: s = s + 1.0
        b.append(s/(k + 1)*UP*UP)
    sv = np.zeros(n)
    for j in range(1, KR + 1):
        sv = (sv + b[j]*rho**j*UP)*UP
    # cauda: com lam = theta/rho, se C = max_{j<=KR} b_j/lam^j <= lam/2 e 4|y0|/(KR+1) <= lam/2, entao
    # b_k <= C lam^k para todo k (inducao) e sum_{k>KR} b_k rho^k <= C theta^{KR+1}/(1-theta)
    best = np.full(n, np.inf)
    for theta in (0.5, 0.7, 0.85, 0.95):
        lam = theta/rho
        Cc = np.zeros(n)
        for j in range(1, KR + 1):
            Cc = np.maximum(Cc, b[j]*(rho/theta)**j*UP)
        ok = (Cc <= lam/2) & (4*ay/(KR + 1) <= lam/2)
        tail = np.where(ok, Cc*theta**(KR + 1)/(1 - theta)*UP, np.inf)
        best = np.minimum(best, tail)
    if not np.all(np.isfinite(best)): raise ValueError('disco de Riccati grande demais')
    dY = (sv + best)*MARG
    return y0, ay, dY, ai

def _exp_up(x):            # >= e^{x}, x float qualquer (Arb)
    return up(Af(x).exp())

def airy3_box(sig, bet, X, Y, n=256, detail=False, rel=1e-20, rows=('mR', 'mL')):
    if rows == ('mR', 'mL') and not lo(bet) > 0: return dict(ok=False, why='beta=0')
    S = 2*sig**(arb(1)/3); s23 = sig**(-arb(2)/3)
    Lam = acb(X, Y)
    sc, bc, xc, yc = arb(sig.mid()), arb(bet.mid()), arb(X.mid()), arb(Y.mid())
    Sc = 2*sc**(arb(1)/3); s23c = sc**(-arb(2)/3); Lamc = acb(xc, yc)
    Shi, Slo = up(S), lo(S)
    Dh = math.nextafter(Shi/n, INF)*(1 + 1e-15)
    bl, bh = max(lo(bet), 0.0), up(bet); rb = up(bet - bc) if True else 0.0
    rb = max(up(bet) - float(bc), float(bc) - lo(bet), 0.0)*MARG
    rs = max(up(sig) - float(sc), float(sc) - lo(sig), 0.0)*MARG
    rS = up(2*arb(lo(sig))**(-arb(2)/3)/3)*rs*MARG
    rLam = up(abs(acb(X, Y) - Lamc))*MARG
    # malha absoluta
    tl_ = np.array([k*Dh for k in range(n)]); th_ = np.array([(k + 1)*Dh for k in range(n)])
    tm_ = 0.5*(tl_ + th_)
    full = th_ <= Slo                                    # celulas inteiras dentro de [0, S] para todo S
    inS = tl_ < Shi                                      # celulas que tocam [0, S_hi]
    ebt = np.array([_exp_up(-bl*a) for a in tl_]); eL = np.array([_exp_up(-bl*max(Slo - b_, 0.0)) for b_ in th_])
    t_up = th_.copy(); St_up = np.maximum(Shi - tl_, 0.0)
    # linha r_4' = e^{beta S} r_4 = int_0^S sinh(beta(S-t))/beta F dt  (reescalar uma linha nao altera o criterio)
    #   peso w = sinh(beta u)/beta, u = S - t:  w <= sinh(bh u)/bh;  |dw/dbeta| <= bh u^3 cosh(bh u)/3;  dw/dS = cosh(beta u)
    def _sh(u): return up(Af(u).sinh()/Af(bh)) if bh > 0 else u
    w4 = np.array([_sh(bh*u)/1.0 if bh > 0 else u for u in St_up])*MARG if False else          np.array([(up((Af(bh)*Af(u)).sinh()/Af(bh)) if bh > 0 else u) for u in St_up])*MARG
    ch = np.array([up((Af(bh)*Af(u)).cosh()) for u in St_up])
    dw4 = bh*St_up**3*ch/3*MARG
    e2S = ch                                                           # dw/dS = cosh(beta (S-t))
    W4s = float(np.max(w4))
    C = {}; T = {}; dT = {}; err = {}
    for key, nu, sh in COMPS:
        Lnuc = Lamc + nu*bc**2 + sh*s23c
        rL = (rLam + up(nu)*(2*float(bc)*rb + rb*rb) + up(sh)*up(2*arb(lo(sig))**(-arb(5)/3)/3)*rs)*MARG
        cc = AY.cfac(nu); acc = up(abs(cc)); W = nu**(-arb(1)/3)/(2*arb.pi()); aW = lo(W)
        rho = acc*(Dh/2 + rL)*MARG
        zs = [cc*(acb(float(t)) + I*Lnuc) for t in tm_]; ws = [W3*z for z in zs]
        yA, ayA, dA, aiA = ydisc(zs, np.full(n, rho)); yB, ayB, dB, aiB = ydisc(ws, np.full(n, rho))
        # taxas: Re(-A'/A) >= aA,  Re(B'/B) >= aB   (A'/A = c Y(z), B'/B = c w Y(w), w = e^{-2pi i/3})
        aA = np.array([lo((-cc*y).real) for y in yA]) - acc*dA
        aB = np.array([lo((cc*W3*y).real) for y in yB]) - acc*dB
        lA = acc*(ayA + dA); lB = acc*(ayB + dB)        # |A'/A|, |B'/B| <=
        uA = np.array([_exp_up(-a*Dh) for a in aA]); fA = np.maximum(uA, 1.0)
        uB = np.array([_exp_up(-a*Dh) for a in aB]); fB = np.maximum(uB, 1.0)
        # |R(t)| (t na celula k) <= prod_{j<k} uA_j * fA_k
        cpA = np.concatenate([[1.0], np.cumprod(uA)])*(1 + 1e-12)
        envR = cpA[:-1]*fA*MARG
        # |X(t)| <= prod_{j>k, celula j inteira em [0,S_lo]} uB_j * prod_{celulas em [S_lo,S_hi]} fB * fB_k
        uBe = np.where(full, uB, np.where(inS, fB, 1.0))
        cpB = np.concatenate([np.cumprod(uBe[::-1])[::-1][1:], [1.0]])*(1 + 1e-12)
        envX = cpB*fB*MARG
        envR = np.where(inS, envR, 0.0); envX = np.where(inS, envX, 0.0)
        RS = float(np.prod(np.where(full, uA, np.where(inS, fA, 1.0))))*MARG
        X0 = float(cpB[0]*uBe[0])*MARG
        # |A B| por celula (disco): |P(u)| <= |Ai(z_k) Ai(w_k)| exp(rho (|yA + e^{-2pi i/3} yB| + dA + dB))
        Pk = np.array([up(abs(a*b_)*(Af(rho)*(Af(up(abs(ya + W3*yb))) + Af(da) + Af(db))).exp())
                       for a, b_, ya, yb, da, db in zip(aiA, aiB, yA, yB, dA, dB)])
        Pk = np.where(inS, Pk, 0.0)
        # paredes: P(0), P(S), Y_A(z(0)), Y_B(w(S))
        z0 = cc*(I*Lnuc); zS = cc*(Sc + I*Lnuc)
        rho0 = acc*rL*MARG; rhoS = acc*(rL + rS)*MARG
        yw, ayw, dw, aiw = ydisc([z0, W3*z0, zS, W3*zS], np.array([rho0, rho0, rhoS, rhoS]))
        P0 = up(abs(aiw[0]*aiw[1])*(Af(rho0)*(Af(up(abs(yw[0] + W3*yw[1]))) + Af(dw[0]) + Af(dw[1]))).exp())
        PS = up(abs(aiw[2]*aiw[3])*(Af(rhoS)*(Af(up(abs(yw[2] + W3*yw[3]))) + Af(dw[2]) + Af(dw[3]))).exp())
        # Green: |G(s,t)| <= P_min f_min prod_between u f_max / |W|
        def gweighted(wt):
            a = wt*Pk*fA; n_ = len(a); ul = uA.tolist(); al = a.tolist(); fl = fA.tolist(); wl = (wt*fA).tolist()
            Lf = [0.0]*n_; Q = [0.0]*n_
            for j in range(1, n_): Lf[j] = Lf[j - 1]*ul[j - 1] + al[j - 1]
            for j in range(n_ - 2, -1, -1): Q[j] = wl[j + 1] + ul[j + 1]*Q[j + 1]
            Lf = np.array(Lf); Q = np.array(Q)
            return (fA*Lf + wt*Pk*fA*fA + Pk*fA*Q)*Dh/aW*MARG*(1 + 1e-12)
        ones = np.where(inS, 1.0, 0.0)
        Grow = gweighted(ones)
        Ninf = float(np.max(Grow))*MARG
        g0 = P0*envR/aW*MARG; gS = PS*envX/aW*MARG
        C[key] = dict(nu=nu, Ninf=Ninf, envR=envR, envX=envX, Rinf=float(np.max(envR)), Xinf=float(np.max(envX)),
                      R1=float(np.sum(envR))*Dh*MARG, X1n=float(np.sum(envX))*Dh*MARG, RS=RS, X0=X0, Dh=Dh,
                      G0=float(np.sum(g0))*Dh*MARG, GS=float(np.sum(gS))*Dh*MARG, g0=g0, gS=gS,
                      WR=gweighted(np.where(inS, ebt, 0.0)), WL=gweighted(np.where(inS, eL, 0.0)))
        if rows == ('mL', 'r4'): C[key]['W4'] = gweighted(np.where(inS, w4, 0.0))
        # termos principais no centro
        N_A = z0.airy_ai(); N_B = (W3*zS).airy_ai()
        side1 = 'R' if rows[0] == 'mR' else 'r4'
        r0, r1 = rows if rows == ('mR', 'mL') else ('r4', 'mL')
        m1A = AY._moment(cc, Lnuc, bc, Sc, side1, 'A', rel); mLA = AY._moment(cc, Lnuc, bc, Sc, 'L', 'A', rel)
        m1B = AY._moment(cc, Lnuc, bc, Sc, side1, 'B', rel); mLB = AY._moment(cc, Lnuc, bc, Sc, 'L', 'B', rel)
        T[(r0, 'R', key)] = m1A/N_A; T[('mL', 'R', key)] = mLA/N_A
        T[(r0, 'X', key)] = m1B/N_B; T[('mL', 'X', key)] = mLB/N_B
        # variacao (valor medio): d/dLnu log R(t) = i c (Y(z_t) - Y(z_0)),  d/dLnu log X(t) = i c w (Y(w_t) - Y(w_S))
        DA = np.array([up(abs(ya - yw[0])) for ya in yA]) + dA + dw[0]
        DB = np.array([up(abs(yb - yw[3])) for yb in yB]) + dB + dw[3]
        lBS = acc*(ayw[3] + dw[3])
        rowlist = (('mR', ebt, t_up), ('mL', eL, St_up)) if rows == ('mR', 'mL') else (('r4', w4, None), ('mL', eL, St_up))
        for row, wv, wt in rowlist:
            for fam, env, Dd in (('R', envR, DA), ('X', envX, DB)):
                m_up = float(np.sum(wv*env))*Dh*MARG
                dL = float(np.sum(wv*env*acc*Dd))*Dh*MARG
                db_ = float(np.sum((dw4 if row == 'r4' else wv*wt)*env))*Dh*MARG
                if row == 'r4':
                    if fam == 'R': dS = float(np.sum(e2S*env))*Dh*MARG
                    else: dS = W4s*X0 + float(np.sum(wv*env*acc*Dd))*Dh*MARG     # tau = S - t (peso so' de tau)
                elif fam == 'R':
                    dS = (RS*_exp_up(-bl*Slo) if row == 'mR' else RS + bh*m_up)
                else:
                    # tau = S - t:  m(X)(S) = int_0^S w X_S(S - tau) d tau,  d/dS X_S(u) = X(u) (B'/B(u) - B'/B(S))
                    corr = float(np.sum(wv*env*acc*Dd))*Dh*MARG
                    dS = (X0 + bh*m_up + corr if row == 'mR' else X0*_exp_up(-bl*Slo) + corr)
                dT[(row, fam, key)] = (dL*rL + db_*rb + dS*rS)*MARG
                err[(row, fam, key)] = 0.0
    R2 = rows if rows == ('mR', 'mL') else ('mL', 'r4')
    Dc = [[c*T[(row, fam, 'g')] + cm*T[(row, fam, '1')] for fam in 'RX'] for row in R2]
    dD = np.array([[up(c)*dT[(row, fam, 'g')] + up(cm)*dT[(row, fam, '1')] for fam in 'RX'] for row in R2])*MARG
    detc = Dc[0][0]*Dc[1][1] - Dc[0][1]*Dc[1][0]
    absD = [[up(abs(Dc[i][j]))*MARG for j in range(2)] for i in range(2)]
    dlo = lo(abs(detc))*(1 - 1e-9)
    absDb = [[absD[i][j] + dD[i][j] for j in range(2)] for i in range(2)]
    ddet = (absD[0][0]*dD[1][1] + absD[1][1]*dD[0][0] + dD[0][0]*dD[1][1] +
            absD[0][1]*dD[1][0] + absD[1][0]*dD[0][1] + dD[0][1]*dD[1][0])*MARG
    main = {}
    for k_ in T:
        main[k_] = (up(abs(T[k_])) + dT[k_])*MARG if k_[2] == '1' else T[k_]
    xover = AY.cross_airy(C, s23)
    if xover is not None and R2 == ('mL', 'r4'): xover = {k: v for k, v in xover.items() if k[2] != 'mR'}
    r = AA.assemble(C, main, err, s23, (W4s if R2 == ('mL', 'r4') else 1.0), detail, Dbounds=(absD, dlo, absDb, dlo - ddet),
                    Eadd=dD, rows=R2, xover=xover)
    r['rows'] = 'airy3' + ('' if R2 == ('mR', 'mL') else '-r4')
    return r

def airy3_safe(*a, **k):
    try:
        return airy3_box(*a, **k)
    except ValueError as e:
        return dict(ok=False, why=str(e))
