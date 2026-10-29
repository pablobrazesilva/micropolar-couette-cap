"""Porte RIGOROSO (Arb / python-flint) do criterio da Prop. 2.6 do G2 v2 (regiao a priori em Lambda), terna (1;1;0,5).

   box_eval(sig, bet, X, Y): caixa em (sigma, beta, Re Lambda, Im Lambda) -> cota superior de rho  (ok se rho<1 e rA<1)
   ann_eval(sig, bet, w_lo, w_hi): "anel" {Re Lambda >= -delta', |Lambda_1| in [w_hi^-2, w_lo^-2]} (w_lo pode ser 0:
       cauda infinita), via cotas assintoticas uniformes na variavel w = |Lambda_1|^{-1/2} (secao 3 do relatorio).

   Formulas: as de apriori.check (G2 v2 §2.3), com:
     * sup/inf por celula de [0,S] obtidos avaliando as funcoes em bolas (t-celula x caixa de parametros);
     * Xi(t)=Re int_0^t sqrt(q) >= soma das cotas inferiores de Re sqrt(q) nas celulas anteriores (sem dependencia);
     * integrais = somas superiores (sup na celula x comprimento);
     * a dupla soma de Green em float (numpy) com margem relativa 1e-9 (ver relatorio, §2.4).
   Todos os valores float vem de extremos de bolas arredondados para fora (nextafter)."""
import math, numpy as np
from flint import arb, acb, ctx
ctx.prec = 96
INF = float('inf')
def up(x):
    v = float(x.upper()); return math.nextafter(v, INF)
def lo(x):
    v = float(x.lower()); return math.nextafter(v, -INF)
def iv(a, b):
    a = arb(a) if not isinstance(a, arb) else a; b = arb(b) if not isinstance(b, arb) else b
    return a.union(b)
def A(x): return arb(x)          # float -> arb exato
I = acb(0, 1)
# ---- parametros (1;1;0,5) exatos
g = arb(1)/2; nu1 = arb(2); kap = arb(4); c = arb(4)/3; cm = 2 - c          # c = 2N/(1+N-g) = 4/3
s1 = kap*1/(1 + 1 - g); sgs = kap*(1 - g)/(1 + 1 - g)                        # 8/3, 4/3
k12 = 2*kap*(1 - g)/(1 + 1 - g)**2; k21 = kap/2                                # |k12| = 16/9, |k21| = 2
DLT = arb('0.05')*g**(arb(1)/3)                                                # delta' = 0,05 g^{1/3}
COMPS = (('1', nu1, s1), ('g', g, sgs))
MARG = 1 + 1e-9

def phi_iv(x):
    """phi(x)=(1-e^{-x})/x (phi(0)=1), decrescente em x>=0; x bola com x>=0"""
    xl, xh = max(lo(x), 0.0), up(x)
    def ph(v):
        if v == 0: return arb(1)
        if v < 1e-6: return iv(1 - A(v)/2, arb(1))        # 1 - v/2 <= phi(v) <= 1
        V = A(v); return (1 - (-V).exp())/V
    return iv(ph(xh).lower(), ph(xl).upper())

def dphi_iv(x):
    """phi'(x) = (e^{-x}(1+x) - 1)/x^2, crescente de -1/2 (x=0) a 0;  x bola com x >= 0"""
    xl, xh = max(lo(x), 0.0), up(x)
    def dp(v):
        if v < 1e-3: return iv(arb(-1)/2, arb(-1)/2 + A(v)/3)      # -1/2 <= phi'(v) <= -1/2 + v/3
        V = A(v); return ((-V).exp()*(1 + V) - 1)/V**2
    return iv(dp(xl).lower(), min(dp(xh).upper(), arb(0)))

def comp_cells(L, nu, S, bet, n):
    """quantidades por celula e de fronteira para uma componente"""
    sq_nu = nu.sqrt()
    r = np.empty(n); alo = np.empty(n); sqh = np.empty(n); aqlo = np.empty(n); aplo = np.empty(n)
    rat0 = np.empty(n); ratS = np.empty(n); eL = np.empty(n); ebt = np.empty(n); St_hi = np.empty(n)
    LS = L - I*S; absL = abs(L); absLS = abs(LS); b2 = bet**2
    for k in range(n):
        xk = iv(arb(k)/n, arb(k + 1)/n)
        t = S*xk; St = S*iv(arb(n - k - 1)/n, arb(n - k)/n)
        z = L - I*t; a = abs(z); q = z/nu; sq = q.sqrt()
        r[k] = max(lo(sq.real), 0.0); alo[k] = lo(a); sqh[k] = up(abs(sq)); aqlo[k] = lo(a/nu)
        aplo[k] = max(lo(abs(q - b2)), 0.0)
        rat0[k] = up(absL/a); ratS[k] = up(absLS/a)
        eL[k] = up((-bet*St).exp()); ebt[k] = up((-bet*t).exp()); St_hi[k] = up(St)
    return dict(r=r, alo=alo, sqh=sqh, aqlo=aqlo, aplo=aplo, rat0=rat0, ratS=ratS, eL=eL, ebt=ebt, St_hi=St_hi,
                absL=absL, absLS=absLS)

# ---- versao vetorizada (numpy) de comp_cells: extremos EXATOS das funcoes monotonas na celula x caixa,
#      com arredondamento para fora apos cada operacao IEEE (+,-,*,/,sqrt corretamente arredondadas);
#      exp com margem relativa 1e-13 (libm com erro < 1 ulp).
_DN = lambda v: np.nextafter(v, -np.inf)
_UP = lambda v: np.nextafter(v, np.inf)
def comp_cells_f(L, nu, S, bet, n):
    nuf = float(nu.mid()); assert arb(nuf) == nu                    # nu = 2 ou 1/2 (exato em float)
    xl, xh = lo(L.real), up(L.real); yl, yh = lo(L.imag), up(L.imag)
    Sl, Sh = lo(S), up(S); bl, bh = max(lo(bet), 0.0), up(bet)
    k = np.arange(n, dtype=float)
    tl = np.maximum(_DN(_DN(Sl*k)/n), 0.0); th = _UP(_UP(Sh*(k + 1))/n)
    # y - t em [yl - th, yh - tl]
    dlo_ = _DN(yl - th); dhi_ = _UP(yh - tl)
    dmin = np.where((dlo_ <= 0) & (dhi_ >= 0), 0.0, np.minimum(np.abs(dlo_), np.abs(dhi_)))
    dmax = np.maximum(np.abs(dlo_), np.abs(dhi_))
    x2min = 0.0 if xl <= 0 <= xh else min(xl*xl, xh*xh)
    x2min = _DN(x2min) if x2min > 0 else 0.0
    x2max = _UP(max(xl*xl, xh*xh))
    az_lo = _DN(np.sqrt(np.maximum(_DN(x2min + _DN(dmin*dmin)), 0.0)))
    az_hi = _UP(np.sqrt(_UP(x2max + _UP(dmax*dmax))))
    # Re sqrt(q) = sqrt((|z| + Re z)/(2 nu)), crescente em Re z e em |Im z|
    azl_x = _DN(np.sqrt(np.maximum(_DN(_DN(xl*xl) + _DN(dmin*dmin)), 0.0)))         # <= sqrt(xl^2 + dmin^2)
    num = _DN(azl_x + xl)
    r = np.where(num > 0, _DN(np.sqrt(np.maximum(_DN(num/(2*nuf)), 0.0))), 0.0)
    alo = np.maximum(az_lo, 0.0)
    sqh = _UP(np.sqrt(_UP(az_hi/nuf)))
    aqlo = _DN(alo/nuf)
    b2l, b2h = _DN(bl*bl), _UP(bh*bh)
    pl = _DN(_DN(xl/nuf) - b2h); ph = _UP(_UP(xh/nuf) + 0.0 - b2l)
    p2min = 0.0 if pl <= 0 <= ph else _DN(min(pl*pl, ph*ph))
    imq = _DN(dmin/nuf)
    aplo = _DN(np.sqrt(np.maximum(_DN(p2min + _DN(imq*imq)), 0.0)))
    absL = abs(L); absLS = abs(L - I*S)
    aLh, aLSh = up(absL), up(absLS)
    with np.errstate(divide='ignore'):
        rat0 = np.where(alo > 0, _UP(aLh/alo), np.inf); ratS = np.where(alo > 0, _UP(aLSh/alo), np.inf)
    Stl = np.maximum(_DN(_DN(Sl*(n - k - 1))/n), 0.0); Sth = _UP(_UP(Sh*(n - k))/n)
    # e^{-beta t} <= q^k em t >= k S_lo/n, q = e^{-beta_lo S_lo/n} (Arb);  potencias por produto acumulado (para cima)
    q = up((-(A(bl)*A(Sl))/n).exp()); qk = cumprod_up(np.full(n, q))
    ebt = qk[:n]; eL = qk[n - 1::-1][:n]                                    # eL[k] = q^{n-k-1}
    return dict(r=r, alo=alo, sqh=sqh, aqlo=aqlo, aplo=aplo, rat0=rat0, ratS=ratS, eL=eL, ebt=ebt, St_hi=Sth,
                absL=absL, absLS=absLS)

# ---- operacoes so com +,-,*,/,sqrt IEEE (corretamente arredondadas) + nextafter para fora; nenhuma libm
def p14_up(a): return _UP(np.sqrt(_UP(np.sqrt(a))))                       # >= a^{1/4}
def p14_dn(a): return _DN(np.sqrt(np.maximum(_DN(np.sqrt(a)), 0.0)))      # <= a^{1/4}
def pm14_up(a): return _UP(1.0/p14_dn(a))                                  # >= a^{-1/4}
def p15_dn(a): return _DN(a*_DN(np.sqrt(a)))                               # <= a^{3/2}
def pm25_up(a): return _UP(1.0/_DN(_DN(a*a)*_DN(np.sqrt(a))))              # >= a^{-5/2}
def expneg_up(x):
    """>= e^{-x} para x >= 0:  e^x >= 1 + x + x^2/2 + x^3/6  (denominador arredondado para baixo)"""
    x = np.maximum(x, 0.0)
    x2 = _DN(x*x); x3 = _DN(x2*x)
    den = _DN(_DN(_DN(1.0 + x) + _DN(x2/2)) + _DN(x3/6))
    return _UP(1.0/den)
def cumprod_up(u):
    """[1, u0, u0 u1, ...] (n+1 termos), cota superior: erro relativo do produto acumulado <= n*2^-53, margem 1e-12"""
    c = np.concatenate([[1.0], np.cumprod(u)])*(1 + 1e-12) + 1e-300       # + 1e-300: cobre arredondamento em subnormais
    return c
USE_F = True
EMAX = 0.6          # Olver: E < 2/3 basta (Green: 2 - 3E > 0; 1 - E > 0)
def box_eval(sig, bet, X, Y, n=None, detail=False):
    S = 2*sig**(arb(1)/3); s23 = sig**(-arb(2)/3)
    Lam = acb(X, Y)
    Dhi, Dlo = up(S)/1.0, lo(S)
    out = dict()
    C = {}
    if n is None:                                      # mesma malha para as duas componentes (cotas cruzadas de Green)
        nn = 1024
        for key, nu, sh in COMPS:
            L = Lam + nu*bet**2 + sh*s23
            aL = up(abs(L)) + up(S); rmax = math.sqrt(aL/lo(nu)); nn = max(nn, int(min(4096, math.ceil(32*rmax*up(S)))))
    else: nn = n
    for key, nu, sh in COMPS:
        L = Lam + nu*bet**2 + sh*s23
        d = (comp_cells_f if USE_F else comp_cells)(L, nu, S, bet, nn)
        Dh = math.nextafter(up(S)/nn, INF)*(1 + 1e-15); Dl = math.nextafter(lo(S)/nn, -INF)*(1 - 1e-15)
        nuf_lo, nuf_hi = lo(nu), up(nu)
        if (d['r'] <= 0).any(): return dict(ok=False, why='q')
        p_ok = not (d['aplo'] <= 0).any()
        # Olver: V <= 5/16 sqrt(nu) int a^{-5/2}
        V = up(arb(5)/16*nu.sqrt()*A(float(np.sum(pm25_up(d['alo'])))*(1 + 1e-12))*A(Dh))
        E = up(A(V).exp() - 1)
        if not (E < EMAX): return dict(ok=False, why="E")
        # u_j >= exp(-(r_j Dl)) por celula;  exp(-Xi(t)) <= prod_{celulas anteriores} u_j ; exp(-(Xi(S)-Xi(t))) <= prod_{seguintes}
        rd = _DN(d['r']*Dl)
        u = expneg_up(rd)
        eXi = cumprod_up(u)[:-1]                                   # eXi[k] >= exp(-Xi(t)), t na celula k
        eXX = cumprod_up(u[::-1])[:-1][::-1]                       # eXX[k] >= exp(-(Xi(S)-Xi(t))), t na celula k
        XiS_lo = float(np.sum(rd))*(1 - 1e-12)                     # <= Xi(S)
        eta_h = _UP(1.0/_DN(4*nuf_lo*p15_dn(d['aqlo'])))*MARG
        fE = (1 + E)/(1 - E)
        r0q = p14_up(d['rat0']); rSq = p14_up(d['ratS'])
        # ---- erro de Olver LOCAL (Olver cap. 6, prova do Teor. 11.1): |d eps/d xi (t)| <= (1+E) int e^{-2|Xi(v)-Xi(t)|} dF(v),
        #      integral para o lado da referencia (A: v in [t,S]; B: v in [0,t]);  dF = 5/16 sqrt(nu) |z|^{-5/2} dt.
        phD = pm25_up(d['alo'])*(5/16*math.sqrt(nuf_hi)*(1 + 1e-15))*Dh            # >= int_celula dF
        u2 = u*u*(1 + 1e-15)
        pl = phD.tolist(); ul = u2.tolist(); nn0 = len(pl)
        ZA = [0.0]*(nn0 + 1); ZB = [0.0]*(nn0 + 1)                                 # ZA[k] = sum_{j>=k} phD_j prod_{i=k}^{j-1} u_i^2
        for k in range(nn0 - 1, -1, -1): ZA[k] = pl[k] + ul[k]*ZA[k + 1]
        for k in range(nn0): ZB[k + 1] = pl[k] + ul[k]*ZB[k]                      # ZB[k+1] = sum_{j<=k} phD_j prod_{i=j+1}^{k} u_i^2
        SA = phD + np.array(ZA[1:]); SB = phD + np.array(ZB[:-1])               # celula propria sem decaimento
        EpA = np.minimum((1 + E)*SA*MARG, E); EpB = np.minimum((1 + E)*SB*MARG, E)  # >= |eps'|/|sqrt q| na celula
        DR = np.minimum(np.cumsum(d['sqh']*EpA*Dh)*MARG, 2*E)                     # >= |eps_A(t) - eps_A(0)|
        DX = np.minimum(np.cumsum((d['sqh']*EpB*Dh)[::-1])[::-1]*MARG, 2*E)       # >= |eps_B(t) - eps_B(S)|
        gR = 1 + DR/(1 - E); gX = 1 + DX/(1 - E)                                  # >= |1+eps(t)|/|1+eps(ref)|
        envR = r0q*eXi*gR*MARG
        envRp = r0q*eXi*d['sqh']*(gR*(1 + eta_h) + EpA/(1 - E))*MARG
        envX = rSq*eXX*gX*MARG
        envXp = rSq*eXX*d['sqh']*(gX*(1 + eta_h) + EpB/(1 - E))*MARG
        # |f' + sqrt(q) f| (R) e |f' - sqrt(q) f| (X):  <= env0 * ( |1+eps|/(4 nu |q|) + |eps'| ) / |1+eps(ref)|
        i4q = _UP(1.0/_DN(4*nuf_lo*d['aqlo']))
        envRr = r0q*eXi*(gR*i4q + d['sqh']*EpA/(1 - E))*MARG
        envXr = rSq*eXX*(gX*i4q + d['sqh']*EpB/(1 - E))*MARG
        bq = _DN(max(lo(bet), 0.0) + d['r'])                         # <= |beta + sqrt q|  (Re sqrt q >= r)
        sqlo = _DN(np.sqrt(d['aqlo']))                               # <= |sqrt q|
        # fronteira
        q0 = L/nu; qS = (L - I*S)/nu; sq0 = q0.sqrt(); sqS = qS.sqrt()
        aq0 = abs(q0); aqS = abs(qS)
        eta0 = up(1/(4*nu*A(lo(aq0))**(arb(3)/2))); etaS = up(1/(4*nu*A(lo(aqS))**(arb(3)/2)))
        etaR = eta0 + E/(1 - E)*MARG; etaX = etaS + E/(1 - E)*MARG
        rq = up(aq0/aqS)
        eXS = up((-A(XiS_lo)).exp())
        RS = up(A(rq)**(arb(1)/4))*eXS/(1 - E)*MARG
        RpS = RS*up(abs(sqS))*(1 + etaS)*MARG
        X0 = up((1/A(lo(aq0/aqS)))**(arb(1)/4))*eXS/(1 - E)*MARG
        Xp0 = X0*up(abs(sq0))*(1 + eta0)*MARG
        # Green: G_ij <= (1+E)^2/(2-3E) hh_i hh_j exp(-(soma estrita entre i e j))
        hh = pm14_up(d['aqlo'])
        nn_ = len(hh)
        # somas de linha  sum_j hh_j exp(-(soma estrita das celulas entre i e j))  em O(n) por recorrencia
        #   R_i = hh_{i+1} + e^{-rd_{i+1}} R_{i+1},   L_i = hh_{i-1} + e^{-rd_{i-1}} L_{i-1}
        er = u
        hl = hh.tolist(); el = er.tolist()
        Rr = [0.0]*nn_; Ll = [0.0]*nn_
        for i in range(nn_ - 2, -1, -1): Rr[i] = hl[i + 1] + el[i + 1]*Rr[i + 1]
        for i in range(1, nn_): Ll[i] = hl[i - 1] + el[i - 1]*Ll[i - 1]
        rows = hh*(hh + np.array(Rr) + np.array(Ll))
        CE = (1 + E)*(1 + E)/(2 - 3*E)*MARG
        Ninf = CE*np.max(rows)*Dh*MARG*(1 + 1e-9)
        C[key] = dict(nu=nu, L=L, d=d, E=E, envR=envR, envRp=envRp, envX=envX, envXp=envXp, sq0=sq0, sqS=sqS, p_ok=p_ok,
                      envRr=envRr, envXr=envXr, bq=bq, sqlo=sqlo, hh=hh, u=u, CE=CE,
                      q0=q0, qS=qS, etaR=etaR, etaX=etaX, RS=RS, RpS=RpS, X0=X0, Xp0=Xp0, Ninf=Ninf, Dh=Dh, n=nn_)
    return finish(C, sig, bet, S, s23, detail, box=(sig, bet, X, Y))

ROWS_MODE = 'both'          # 'old' (m_L, r_4), 'LR' (m_R, m_L) ou 'both' (a melhor das duas; ambas validas)
def finish(C, sig, bet, S, s23, detail=False, box=None):
    for Cc in C.values():
        Cc['Rinf'] = float(np.max(Cc['envR'])); Cc['Xinf'] = float(np.max(Cc['envX']))
        Cc['R1'] = float(np.sum(Cc['envR'])*Cc['Dh']*MARG); Cc['X1n'] = float(np.sum(Cc['envX'])*Cc['Dh']*MARG)
    res = []
    if ROWS_MODE in ('old', 'both') and all(Cc['p_ok'] for Cc in C.values()):
        r = finish_old(C, sig, bet, S, s23, detail, box); r['rows'] = 'mL,r4'; res.append(r)
    if ROWS_MODE in ('LR', 'both') and lo(bet) > 0:
        r = finish_LR(C, sig, bet, S, s23, detail, box); r['rows'] = 'mR,mL'; res.append(r)
    if not res: return dict(ok=False, why='p')
    ok = [r for r in res if r.get('ok')]
    if ok: return min(ok, key=lambda r: r['rho'])
    return min(res, key=lambda r: r.get('rho', INF))

def finish_LR(C, sig, bet, S, s23, detail=False, box=None):
    """linhas m_R = int e^{-beta t} F, m_L = int e^{-beta(S-t)} F (equivalentes a Phi=Phi'=0 nas paredes se beta>0).
       Identidade: para f'=-sqrt(q) f + r,  int e^{-beta t} f = f(0)/(beta+sqrt q0) - e^{-beta S} f(S)/(beta+sqrt qS)
                   - int e^{-beta t} f q'/(2 sqrt q (beta+sqrt q)^2) + int e^{-beta t} r/(beta+sqrt q);  |q'| = 1/nu.
       Denominadores |beta + sqrt q| >= beta + Re sqrt q > 0: nenhuma singularidade (ao contrario de p = q - beta^2)."""
    b_lo = max(lo(bet), 0.0)
    eb = (-bet*S).exp(); ebh = up(eb)
    main = {}; err = {}
    for key, Cc in C.items():
        nu = Cc['nu']; d = Cc['d']; Dh = Cc['Dh']; bq = Cc['bq']
        w2 = 1.0/(2*lo(nu)*Cc['sqlo']*bq*bq)
        sq0, sqS = Cc['sq0'], Cc['sqS']
        bq0 = lo(abs(bet + sq0)); bqS = lo(abs(bet + sqS))
        if bq0 <= 0 or bqS <= 0 or (bq <= 0).any(): return dict(ok=False, why='bq')
        main[('mR', 'R', key)] = 1/(bet + sq0)
        main[('mL', 'X', key)] = 1/(bet + sqS)
        main[('mL', 'R', key)] = acb(0)
        main[('mR', 'X', key)] = acb(0)
        err[('mR', 'R', key)] = (ebh*Cc['RS']/bqS + np.sum(d['ebt']*(Cc['envR']*w2 + Cc['envRr']/bq))*Dh)*MARG
        err[('mL', 'X', key)] = (ebh*Cc['X0']/bq0 + np.sum(d['eL']*(Cc['envX']*w2 + Cc['envXr']/bq))*Dh)*MARG
        err[('mL', 'R', key)] = np.sum(d['eL']*Cc['envR'])*Dh*MARG
        err[('mR', 'X', key)] = np.sum(d['ebt']*Cc['envX'])*Dh*MARG
    xover = cross_bounds(C, s23)
    Dbounds = None; Eadd = None
    if box is not None:
        Dbounds, Cabs, Eadd = centered_bounds(*box, rows=('mR', 'mL'))
        for k, v in Cabs.items(): main[k] = v
    return assemble(C, main, err, s23, 1.0, detail, Dbounds=Dbounds, Eadd=Eadd, rows=('mR', 'mL'), xover=xover)

def _green_weighted(hh, u, w):
    """V_k >= sum_s w_s hh_s exp(-(Xi entre as celulas s e k)) (sem o fator hh_k, CE, Delta), recorrencia O(n)"""
    n = len(hh); a = (w*hh).tolist(); ul = u.tolist()
    Lf = [0.0]*n; Rt = [0.0]*n
    for k in range(1, n): Lf[k] = a[k - 1] + ul[k - 1]*Lf[k - 1]
    for k in range(n - 2, -1, -1): Rt[k] = a[k + 1] + ul[k + 1]*Rt[k + 1]
    return np.array(a) + np.array(Lf) + np.array(Rt)

def cross_bounds(C, s23):
    """cotas pontuais/ponderadas da correcao de acoplamento na OUTRA componente (Lema 2.5 refinado):
       c_o = kap_o K_o[f] + kap_o K_o[kap_own K_own c_o],  |G_o(s,t)| <= CE hh_s hh_t e^{-|Xi(s)-Xi(t)|}.
       parede 0: |c_o(0)| <= k_o P0_o[f] + k_o G0_o k_own Ninf_own ||c_o||_inf ;  idem em S.
       momento m_R: int e^{-beta s}|c_o| <= k_o sum_t W_o(t) f(t) + k_o sup W_o k_own Ninf_own ||c_o||_1 ."""
    k1 = up(k12*s23/nu1); kg = up(k21*s23/g)
    kap = {'1': k1, 'g': kg}
    c1, cg = C['1'], C['g']
    if len(c1['hh']) != len(cg['hh']): return None
    tau = k1*kg*c1['Ninf']*cg['Ninf']*MARG
    if tau >= 1: return None
    out = {}
    for own in ('1', 'g'):
        o = 'g' if own == '1' else '1'; Co = C[o]; Cw = C[own]; Dh = Co['Dh']
        ko = kap[o]; kw = kap[own]
        hh = Co['hh']; CE = Co['CE']; d = Co['d']
        # parede 0 (s na celula 0) e parede S (s na ultima celula)
        u = Co['u']
        dec0 = cumprod_up(np.concatenate([[1.0], u[1:]]))[:-1]           # prod_{i=1}^{k-1} u_i  (k>=1); k=0 -> 1
        dec0[0] = 1.0
        decS = cumprod_up(np.concatenate([[1.0], u[::-1][1:]]))[:-1][::-1]
        decS[-1] = 1.0
        g0 = CE*hh[0]*hh*dec0; gS = CE*hh[-1]*hh*decS                     # >= |G_o(0,t)|, |G_o(S,t)|, t na celula
        G0 = float(np.sum(g0))*Dh*MARG; GS = float(np.sum(gS))*Dh*MARG
        WR = CE*hh*_green_weighted(hh, u, d['ebt'])*Dh*MARG               # >= int e^{-beta s}|G_o(s,t)| ds
        WL = CE*hh*_green_weighted(hh, u, d['eL'])*Dh*MARG
        for fam in 'RX':
            f = Cw['envR'] if fam == 'R' else Cw['envX']
            finf = Cw['Rinf'] if fam == 'R' else Cw['Xinf']; f1 = Cw['R1'] if fam == 'R' else Cw['X1n']
            cinf = ko*Co['Ninf']*finf/(1 - tau)*MARG; c1n = ko*Co['Ninf']*f1/(1 - tau)*MARG
            out[(fam, own, 0)] = ko*(float(np.sum(g0*f))*Dh + G0*kw*Cw['Ninf']*cinf)*MARG
            out[(fam, own, 'S')] = ko*(float(np.sum(gS*f))*Dh + GS*kw*Cw['Ninf']*cinf)*MARG
            out[(fam, own, 'mR')] = ko*(float(np.sum(WR*f))*Dh + float(np.max(WR))*kw*Cw['Ninf']*c1n)*MARG
            out[(fam, own, 'mL')] = ko*(float(np.sum(WL*f))*Dh + float(np.max(WL))*kw*Cw['Ninf']*c1n)*MARG
    return out

def finish_old(C, sig, bet, S, s23, detail=False, box=None):
    """monta o criterio (comum a box_eval); C[key] contem as cotas por componente"""
    b_lo, b_hi = max(lo(bet), 0.0), up(bet)
    eb = (-bet*S).exp(); ebh = up(eb)
    sb = S*phi_iv(2*bet*S); cb = (1 + (-2*bet*S).exp())/2
    sbh = up(sb); cbh = up(cb)
    W4 = min(up(S), 1/(2*b_lo)) if b_lo > 0 else up(S)
    main = {}; err = {}
    for key, Cc in C.items():
        nu = Cc['nu']; d = Cc['d']; Dh = Cc['Dh']
        p0 = Cc['q0'] - bet**2; pS = Cc['qS'] - bet**2
        ap0 = lo(abs(p0)); apS = lo(abs(pS))
        if ap0 <= 0 or apS <= 0: return dict(ok=False, why='p0')
        ip2 = 1/(lo(nu)*d['aplo']*d['aplo'])
        ws = d['ebt']*np.minimum(d['St_hi'], (1/(2*b_lo) if b_lo > 0 else INF))
        wc = d['ebt']
        rhoL_R = np.sum(d['eL']*(Cc['envRp'] + b_hi*Cc['envR'])*ip2)*Dh*MARG
        rho4_R = np.sum((ws*Cc['envRp'] + wc*Cc['envR'])*ip2)*Dh*MARG
        rhoL_X = np.sum(d['eL']*(Cc['envXp'] + b_hi*Cc['envX'])*ip2)*Dh*MARG
        rho4_X = np.sum((ws*Cc['envXp'] + wc*Cc['envX'])*ip2)*Dh*MARG
        sq0, sqS = Cc['sq0'], Cc['sqS']
        main[('mL', 'R', key)] = eb*(sq0 + bet)/p0
        main[('r4', 'R', key)] = (sb*sq0 - cb)/p0
        main[('mL', 'X', key)] = (sqS - bet)/pS
        main[('r4', 'X', key)] = eb/pS
        asq0 = up(abs(sq0)); asqS = up(abs(sqS))
        err[('mL', 'R', key)] = ((Cc['RpS'] + b_hi*Cc['RS'])/apS + ebh*asq0*Cc['etaR']/ap0 + rhoL_R)*MARG
        err[('r4', 'R', key)] = (ebh*Cc['RS']/apS + sbh*asq0*Cc['etaR']/ap0 + rho4_R)*MARG
        err[('mL', 'X', key)] = (asqS*Cc['etaX']/apS + ebh*(Cc['Xp0'] + b_hi*Cc['X0'])/ap0 + rhoL_X)*MARG
        err[('r4', 'X', key)] = ((sbh*Cc['Xp0'] + cbh*Cc['X0'])/ap0 + rho4_X)*MARG
    Dbounds = None; Eadd = None
    if box is not None:
        Dbounds, Cabs, Eadd = centered_bounds(*box)
        for k, v in Cabs.items(): main[k] = v
    xover = cross_bounds(C, s23)
    if xover is not None: xover = {k: v for k, v in xover.items() if k[2] != 'mR'}     # paredes e m_L (mesmo peso)
    return assemble(C, main, err, s23, W4, detail, Dbounds=Dbounds, Eadd=Eadd, xover=xover)

def assemble(C, main, err, s23, W4, detail=False, rowscale=None, Dbounds=None, Eadd=None, rows=('mL', 'r4'), xover=None):
    """Schur 2x2 (Prop. 2.6).  main: bolas acb (ou dict de cotas para o anel); err: floats >= 0"""
    cf = {'1': 1.0, 'g': up(c)}; cfo = {'1': up(c), 'g': 1.0}
    k1 = up(k12*s23/nu1); kg = up(k21*s23/g)
    c1, cg = C['1'], C['g']
    tauI = k1*kg*c1['Ninf']*cg['Ninf']*MARG; tau1 = tauI
    if tauI >= 1: return dict(ok=False, why='tau')
    cmh = up(cm)
    def col_errs(fam, own):
        o = C[own]; ot = C['g' if own == '1' else '1']; kt = kg if own == '1' else k1
        finf = o['Rinf'] if fam == 'R' else o['Xinf']; f1 = o['R1'] if fam == 'R' else o['X1n']
        return (tauI*finf/(1 - tauI)*MARG, kt*ot['Ninf']*finf/(1 - tauI)*MARG,
                tau1*f1/(1 - tau1)*MARG, kt*ot['Ninf']*f1/(1 - tau1)*MARG)
    E = {(f, o): col_errs(f, o) for f in 'RX' for o in ('1', 'g')}
    def pt_err(fam, own, wall):
        o = C[own]; oi, ti, _, _ = E[(fam, own)]
        coef = 1.0 if own == '1' else cmh; coefo = cmh if own == '1' else 1.0
        if xover is not None and (fam, own, wall) in xover: ti = min(ti, xover[(fam, own, wall)])
        e = coef*oi + coefo*ti
        own_wall = (fam == 'R' and wall == 0) or (fam == 'X' and wall == 'S')
        if own_wall: return e*MARG
        vabs = o['RS'] if fam == 'R' else o['X0']
        return (coef*vabs + e)*MARG
    def mom(fam, own, row):
        _, _, o1, t1 = E[(fam, own)]
        wsup = W4 if row == 'r4' else 1.0                 # sup dos pesos: e^{-beta t}, e^{-beta(S-t)} <= 1
        tt = wsup*t1
        if xover is not None and (fam, own, row) in xover: tt = min(tt, xover[(fam, own, row)])
        e = (cf[own]*err[(row, fam, own)] + wsup*cf[own]*o1 + cfo[own]*tt)*MARG
        return e
    # A - I (so' erros), B (so' erros), C, D
    dA = np.array([[pt_err('R', '1', 0), pt_err('X', '1', 0)], [pt_err('R', '1', 'S'), pt_err('X', '1', 'S')]])
    Bt = np.zeros((2, 2))
    for i, wall in enumerate((0, 'S')):
        for j, fam in enumerate('RX'):
            Bt[i, j] = (pt_err(fam, 'g', wall) + cmh*pt_err(fam, '1', wall))*MARG
    if rowscale is None: rowscale = (1.0, 1.0)
    Ct = np.zeros((2, 2)); De = np.zeros((2, 2)); Dm = [[None, None], [None, None]]
    for i, row in enumerate(rows):
        for j, fam in enumerate('RX'):
            mabs = main[(row, fam, '1')] if isinstance(main[(row, fam, '1')], float) else up(abs(main[(row, fam, '1')]))
            Ct[i, j] = (mabs + mom(fam, '1', row))*MARG
            De[i, j] = (mom(fam, 'g', row) + cmh*mom(fam, '1', row))*MARG
            if Dbounds is None: Dm[i][j] = c*main[(row, fam, 'g')] + cm*main[(row, fam, '1')]
    # rho(dA) e (I-dA)^{-1} (2x2 nao negativa)
    a, b_, c_, d_ = dA[0, 0], dA[0, 1], dA[1, 0], dA[1, 1]
    detIA = (1 - a)*(1 - d_) - b_*c_
    if not (a < 1 and d_ < 1 and detIA > 0): return dict(ok=False, why='A')
    Ainv = np.array([[1 - d_, b_], [c_, 1 - a]])/detIA*(1 + 1e-9)
    rA = (0.5*(a + d_) + math.sqrt(0.25*(a - d_)*(a - d_) + b_*c_))*(1 + 1e-9)
    Delta = (De + Ct @ Ainv @ Bt)*(1 + 1e-9)
    def rho_of(absD, dlo, Dl):
        if dlo <= 0: return INF
        Dinv = np.array([[absD[1][1], absD[0][1]], [absD[1][0], absD[0][0]]])/dlo*(1 + 1e-9)
        Mx = (Dinv @ Dl)*(1 + 1e-9)
        a, b_, c_, d_ = Mx[0, 0], Mx[0, 1], Mx[1, 0], Mx[1, 1]
        return (0.5*(a + d_) + math.sqrt(0.25*(a - d_)*(a - d_) + b_*c_))*(1 + 1e-9), Dinv
    if Dbounds is None:
        det0 = Dm[0][0]*Dm[1][1] - Dm[0][1]*Dm[1][0]
        cands = [([[up(abs(Dm[i][j])) for j in range(2)] for i in range(2)], lo(abs(det0)), Delta)]
    elif len(Dbounds) == 2:
        cands = [(Dbounds[0], Dbounds[1], Delta)]
    else:
        # (i) D_m(centro) com a variacao de D_m como perturbacao extra;  (ii) sup|D_m|, inf|det| na caixa.  Ambas validas.
        absD, dlo, absDb, dlob = Dbounds
        cands = [(absD, dlo, (Delta + Eadd)*(1 + 1e-9)), (absDb, dlob, Delta)]
    best = (INF, None)
    for absD_, dlo_, Dl in cands:
        if dlo_ <= 0: continue
        rr = rho_of(absD_, dlo_, Dl)
        if rr[0] < best[0]: best = rr
    if best[1] is None: return dict(ok=False, why='det0')
    rho, Dinv = best
    res = dict(ok=bool(rho < 1 and rA < 1), rho=rho, rA=rA, tau=tauI)
    if detail: res.update(Dinv=Dinv, Delta=Delta, De=De, Ct=Ct, Bt=Bt, dA=dA)
    return res


# ---------------------------------------------------------------------------------------------------------
def ann_eval(sig, bet, w_lo, w_hi):
    """cotas uniformes para todo Lambda com Re Lambda >= -delta' e |Lambda_1| in [w_hi^-2, w_lo^-2] (w_lo=0: cauda).
       Linhas de momento escaladas por 1/w (rho e' invariante; ver relatorio §3)."""
    S = 2*sig**(arb(1)/3); s23 = sig**(-arb(2)/3)
    w = iv(A(w_lo), A(w_hi)); wh = up(w); th = iv(-1, 1)
    d = (nu1 - g)*bet**2 + (s1 - sgs)*s23                       # Lambda_g = Lambda_1 - d, d >= 0
    zmf = 1 - (d + S)*w**2                                       # |z_nu(s)| >= zmf / w^2
    if lo(zmf) <= 0: return dict(ok=False, why='zmf')
    epsA = DLT*w**2/zmf
    if up(epsA) >= 1: return dict(ok=False, why='eps')
    eps = epsA.asin()
    rS = S*w**2/(1 - d*w**2)
    if up(rS) >= 1: return dict(ok=False, why='rS')
    klo = 1 - rS; khi = 1 + rS
    cang = (arb.pi()/4 + eps/2).cos()
    gd = (d*w**2/zmf)
    if up(gd) >= 1: return dict(ok=False, why='gd')
    gam_d = gd.asin()
    b_lo, b_hi = max(lo(bet), 0.0), up(bet)
    eb = (-bet*S).exp(); sb = S*phi_iv(2*bet*S); cb = (1 + (-2*bet*S).exp())/2
    W4 = min(up(S), 1/(2*b_lo)) if b_lo > 0 else up(S)
    C = {}; main = {}; err = {}; T = {}
    for key, nu, sh in COMPS:
        sh_ = (1/nu).sqrt() if key == '1' else ((1 + th*d*w**2)/nu).sqrt()     # w sqrt(Q_nu) (exato em |q(0)|)
        shi = sh_*khi.sqrt()                                                    # >= w |sqrt q(t)|
        rh = sh_*klo.sqrt()*cang                                                # <= w Re sqrt q(t)
        if lo(rh) <= 0: return dict(ok=False, why='rh')
        V = arb(5)/16*nu.sqrt()*S*w**5/zmf**(arb(5)/2)
        E = V.exp() - 1
        if up(E) >= 0.3: return dict(ok=False, why='E')
        etab = w**3*nu.sqrt()/(4*zmf**(arb(3)/2))
        Fexp = (-(A(lo(rh))*A(lo(S)))/A(wh)).exp() if wh > 0 else arb(0)
        Fexp = Fexp*MARG
        hr = klo**(-arb(1)/4)
        hr = arb(max(up(hr), up(khi**(arb(1)/4))))
        KE = hr*(1 + E)/(1 - E); KEp = hr*((1 + 2*E) + (1 + E)*etab)/(1 - E)
        RSp = hr*Fexp/(1 - E)
        Ninf = (1 + E)**2/(2 - 3*E)*(nu/zmf).sqrt()*w*2*w/rh
        P2 = (zmf - nu*bet**2*w**2)
        if lo(P2) <= 0: return dict(ok=False, why='P2')
        iP2w = nu*w**3/P2**2                                    # (1/(nu|p|^2))/w  <= nu w^3/P2^2
        D0 = sh_**2 - bet**2*w**2; DS = sh_**2*klo - bet**2*w**2
        if lo(D0) <= 0 or lo(DS) <= 0: return dict(ok=False, why='D')
        etaR = etab + E/(1 - E)
        rhoL = iP2w*(KEp*shi/rh + bet*KE*w/rh)
        rho4 = iP2w*(W4*KEp*shi/rh + KE*w/rh)
        err[('mL', 'R', key)] = up(RSp*(shi*(1 + etab) + bet*w)/DS + eb*sh_*etaR/D0 + rhoL)
        err[('r4', 'R', key)] = up(eb*RSp*w/DS + sb*sh_*etaR/D0 + rho4)
        err[('mL', 'X', key)] = up(shi*etaR/DS + eb*RSp*(sh_*(1 + etab) + bet*w)/D0 + rhoL)
        err[('r4', 'X', key)] = up(RSp*(sb*sh_*(1 + etab) + cb*w)/D0 + rho4)
        # cotas dos termos principais (escalados por 1/w)
        ef = (cb*w/(sb*sh_) + bet**2*w**2/sh_**2)/(1 - bet**2*w**2/sh_**2)
        if up(ef) >= 1: return dict(ok=False, why='ef')
        if up(bet*w/(sh_*klo.sqrt())) >= 1 or lo(sh_ - bet*w) <= 0: return dict(ok=False, why='bw')
        T[key] = dict(mLR_hi=up(eb/(sh_ - bet*w)), r4R_hi=up(sb*(1 + ef)/sh_), r4R_lo=lo(sb*(1 - ef)/sh_),
                      mLX_hi=up(1/(sh_*klo.sqrt())), mLX_lo=lo(1/(shi + bet*w)), r4X_hi=up(eb*w/DS),
                      ef=ef, bang=(bet*w/(sh_*klo.sqrt())).asin())
        for k_ in ('mLR', 'r4R', 'mLX', 'r4X'): pass
        main[('mL', 'R', key)] = T[key]['mLR_hi']; main[('r4', 'R', key)] = T[key]['r4R_hi']
        main[('mL', 'X', key)] = T[key]['mLX_hi']; main[('r4', 'X', key)] = T[key]['r4X_hi']
        C[key] = dict(nu=nu, Ninf=up(Ninf), Rinf=up(KE), Xinf=up(KE), R1=up(KE/rh), X1n=up(KE/rh),   # R1, X1n: ja divididos por w
                      RS=up(RSp), X0=up(RSp))
    ch = up(c); cmh_ = up(cm); cl = lo(c); cml = lo(cm)
    # angulos entre os termos nu=1 e nu=g
    A_S = gam_d/2 + T['1']['bang'] + T['g']['bang']
    A_0 = gam_d/2 + T['1']['ef'].asin() + T['g']['ef'].asin()
    if up(A_S) >= 3 or up(A_0) >= 3: return dict(ok=False, why='ang')
    mLS_lo = lo((A_S/2).cos())*(cl*T['g']['mLX_lo'] + cml*T['1']['mLX_lo'])
    r40_lo = lo((A_0/2).cos())*(cl*T['g']['r4R_lo'] + cml*T['1']['r4R_lo'])
    absD = [[ch*T['g']['mLR_hi'] + cmh_*T['1']['mLR_hi'], ch*T['g']['mLX_hi'] + cmh_*T['1']['mLX_hi']],
            [ch*T['g']['r4R_hi'] + cmh_*T['1']['r4R_hi'], ch*T['g']['r4X_hi'] + cmh_*T['1']['r4X_hi']]]
    absD = [[x*MARG for x in r] for r in absD]
    dlo = (mLS_lo*r40_lo - absD[0][0]*absD[1][1]*MARG)*(1 - 1e-9)
    # o 'main' passado a assemble e' |.| (float) para as colunas R_1, X_1 (linha escalada)
    return assemble(C, main, err, s23, W4, Dbounds=(absD, dlo))

# ---------------------------------------------------------------------------------------------------------
# Forma centrada (valor medio) para os termos principais D_m e as colunas C: f(caixa) c f(centro) + sum_i d_i f(caixa) [-r_i, r_i]
class Du:
    """numero dual com 4 direcoes (Re Lam, Im Lam, beta, sigma); valores acb (bolas)"""
    __slots__ = ('v', 'g')
    def __init__(s, v, g=None):
        s.v = v if isinstance(v, acb) else acb(v); s.g = g if g is not None else [acb(0)]*4
    @staticmethod
    def lift(x): return x if isinstance(x, Du) else Du(x)
    def __add__(s, o): o = Du.lift(o); return Du(s.v + o.v, [a + b for a, b in zip(s.g, o.g)])
    __radd__ = __add__
    def __neg__(s): return Du(-s.v, [-a for a in s.g])
    def __sub__(s, o): return s + (-Du.lift(o))
    def __rsub__(s, o): return Du.lift(o) + (-s)
    def __mul__(s, o):
        o = Du.lift(o); return Du(s.v*o.v, [a*o.v + s.v*b for a, b in zip(s.g, o.g)])
    __rmul__ = __mul__
    def inv(s):
        iv_ = 1/s.v; return Du(iv_, [-a*iv_*iv_ for a in s.g])
    def __truediv__(s, o): return s*Du.lift(o).inv()
    def __rtruediv__(s, o): return Du.lift(o)*s.inv()
    def sqrt(s):
        r = s.v.sqrt(); return Du(r, [a/(2*r) for a in s.g])
    def exp(s):
        e = s.v.exp(); return Du(e, [a*e for a in s.g])
    def pw(s, p):                       # real power of a positive real dual
        r = s.v**p; return Du(r, [a*p*s.v**(p - 1) for a in s.g])

def main_entries(sig, bet, Lam, use_dual, rows=('mL', 'r4')):
    """D_m (2x2) e |termos de C| como Du (use_dual) ou acb (centro)"""
    if use_dual:
        Lam_ = Du(Lam, [acb(1), acb(0, 1), acb(0), acb(0)])
        bet_ = Du(acb(bet), [acb(0), acb(0), acb(1), acb(0)])
        sig_ = Du(acb(sig), [acb(0), acb(0), acb(0), acb(1)])
    else:
        Lam_, bet_, sig_ = Du(Lam), Du(acb(bet)), Du(acb(sig))
    S = 2*sig_.pw(arb(1)/3); s23 = sig_.pw(-arb(2)/3)
    x2 = 2*bet_*S
    eb = (-(bet_*S)).exp(); cb = (1 + (-x2).exp())*(arb(1)/2)
    # sb = S phi(2 beta S): phi via phi_iv (valor) e phi' em [-1/2, 0]
    xr = x2.v.real
    phv = acb(phi_iv(xr))
    dph = acb(dphi_iv(xr)) if use_dual else acb(0)
    sb = Du(S.v*phv, [a*phv + S.v*dph*b for a, b in zip(S.g, x2.g)])
    T = {}
    for key, nu, sh in COMPS:
        L = Lam_ + nu*bet_*bet_ + sh*s23
        q0 = L/nu; qS = (L - I*S)/nu; sq0 = q0.sqrt(); sqS = qS.sqrt()
        if rows == ('mR', 'mL'):
            T[('mR', 'R', key)] = (bet_ + sq0).inv(); T[('mL', 'X', key)] = (bet_ + sqS).inv()
            T[('mL', 'R', key)] = Du(acb(0)); T[('mR', 'X', key)] = Du(acb(0))
            continue
        p0 = q0 - bet_*bet_; pS = qS - bet_*bet_
        T[('mL', 'R', key)] = eb*(sq0 + bet_)/p0
        T[('r4', 'R', key)] = (sb*sq0 - cb)/p0
        T[('mL', 'X', key)] = (sqS - bet_)/pS
        T[('r4', 'X', key)] = eb/pS
    Dm = [[c*T[(row, fam, 'g')] + cm*T[(row, fam, '1')] for fam in 'RX'] for row in rows]
    det = Dm[0][0]*Dm[1][1] - Dm[0][1]*Dm[1][0]
    return T, Dm, det

def centered_bounds(sig, bet, X, Y, rows=('mL', 'r4')):
    """|D_m| sup, |det| inf, |C principal| sup pela forma do valor medio"""
    sc, bc, xc, yc = arb(sig.mid()), arb(bet.mid()), arb(X.mid()), arb(Y.mid())
    rads = [float(X.rad()), float(Y.rad()), float(bet.rad()), float(sig.rad())]
    rads = [math.nextafter(r*(1 + 1e-12), INF) for r in rads]
    Tc, Dc, detc = main_entries(sc, bc, acb(xc, yc), False, rows)
    Tb, Db, detb = main_entries(sig, bet, acb(X, Y), True, rows)
    def hi(cv, bv):
        return up(abs(cv.v)) + sum(up(abs(gi))*r for gi, r in zip(bv.g, rads))
    def lw(cv, bv):
        return lo(abs(cv.v)) - sum(up(abs(gi))*r for gi, r in zip(bv.g, rads))
    # D_m(x) = D_m(c) + E(x),  |E_jk| <= sum_i |d_i D_jk(caixa)| r_i  -> entra como perturbacao adicional (Eadd)
    absD = [[up(abs(Dc[i][j].v))*MARG for j in range(2)] for i in range(2)]
    dlo = lo(abs(detc.v))*(1 - 1e-9)
    Eadd = np.array([[sum(up(abs(gi))*r for gi, r in zip(Db[i][j].g, rads))*MARG for j in range(2)] for i in range(2)])
    # alternativa (ii): sup|D_m| e inf|det| na caixa (forma centrada e bola direta; a menor/maior)
    absDb = [[min(hi(Dc[i][j], Db[i][j])*MARG, up(abs(Db[i][j].v))) for j in range(2)] for i in range(2)]
    dlob = max(lw(detc, detb)*(1 - 1e-9), lo(abs(detb.v)))
    Cabs = {k: min(hi(Tc[k], Tb[k])*MARG, up(abs(Tb[k].v))) for k in Tc if k[2] == '1'}
    return (absD, dlo, absDb, dlob), Cabs, Eadd
