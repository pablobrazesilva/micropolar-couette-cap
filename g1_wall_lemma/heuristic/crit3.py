"""Criterio do Teorema C, versao fina (crit3): funcoes de dois pontos (tp.py) + rota 3 (funcao G deslizante, Lema F4)
   + Lema 4.1 com OmPsi.  Para cada quantidade usa o minimo entre a cota nova e a cota de crit2 (ambas validas).
   [PROVADO* condicionado a certificacao das tabelas de tp.py]"""
import numpy as np, functools
import crit2 as C2, tp
K = C2.CONST
MARG = 1.01                         # margem nas quadraturas novas (termos somados)
MSUB = 1.0005                       # margem para integrais que entram subtraidas (erro de trapezio ~ (beta ds)^2/12)
def trap(y, x): return np.trapezoid(y, x)*MARG
import os, json
BGRID = [0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3, 5, 8, 12, 20, 40]
_JF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Jcache.json')
_JC = json.load(open(_JF)) if os.path.exists(_JF) else {}
def Jc(nu_k, nu_f, kind, beta=0.0):
    key = '%.6f|%.6f|%s|%g' % (nu_k, nu_f, kind, beta)
    if key not in _JC:
        _JC[key] = tp.J(nu_k, nu_f, kind, beta=beta); json.dump(_JC, open(_JF, 'w'))
    return _JC[key]*1.005          # J calculado com OmPsi numerico (ate 0,45% abaixo do certificado)
def JGb(nu_k, nu_f, beta):
    """J_G com peso e^{-beta s}: nao-crescente em beta -> valor no no' da grade abaixo de beta"""
    bk = max(b for b in BGRID if b <= beta)
    return Jc(nu_k, nu_f, '1G', bk)
class Prof:
    """perfis independentes de beta, para (N,J,G,sigma)"""
    def __init__(s, N, J, G, sig, ds=0.01):
        P = C2.params(N, J, G); s.P = P; s.sig = sig
        s.g, s.nu1, s.c = P['g'], P['nu1'], P['c']
        s.S = 2*sig**(1/3); s23 = sig**(-2/3)
        s.Sc = min(s.S, 90.0)                      # integracao numerica ate Sc (alem: caudas desprezaveis, ver relatorio)
        s.x = np.arange(0, s.Sc + ds/2, ds); s.ds = ds
        s.B = {}
        for nm, nu in (('1', s.nu1), ('g', s.g)):
            s.B[nm] = dict(nu=nu, r=K['a']*nu**(-1/3), at=K['at']*nu**(-1/3), C=2*np.pi*nu**(-2/3)*s23)
        x = s.x
        s.OmA = {k: tp.F('OmA', b['nu']**(-1/3)*x, rate_tail=K['at']) for k, b in s.B.items()}
        s.OmP = {k: tp.F('OmPsi', b['nu']**(-1/3)*x, rate_tail=K['a']) for k, b in s.B.items()}
        s.KP = {k: tp.F('KPi', b['nu']**(-1/3)*x, rate_tail=K['at']) for k, b in s.B.items()}
        s.KP1 = {k: tp.F('KPi1', b['nu']**(-1/3)*x, rate_tail=K['at']) for k, b in s.B.items()}
        # minimo com as cotas exponenciais uniformes (ambas validas): OmA<=e^{-at l}, OmPsi<=e^{-a l}, KPi<=Pi e^{-at l}, KPi1<=Pi1 e^{-at l}
        for k, b in s.B.items():
            l = b['nu']**(-1/3)*x
            s.OmA[k] = np.minimum(s.OmA[k], np.exp(-K['at']*l)); s.OmP[k] = np.minimum(s.OmP[k], np.exp(-K['a']*l))
            s.KP[k] = np.minimum(s.KP[k], K['Pi']*np.exp(-K['at']*l)); s.KP1[k] = np.minimum(s.KP1[k], K['Pi1']*np.exp(-K['at']*l))
        s.Q = {k: tp.F('Q', b['nu']**(-1/3)*x, rate_tail=K['at']) for k, b in s.B.items()}
        s.Qs = {k: tp.F('Qs', b['nu']**(-1/3)*x, rate_tail=K['at']) for k, b in s.B.items()}
        # pares (kernel k, forcing f): ('g','1') -> u ;  ('1','g') -> v
        s.pair = {}
        for k, f in (('g', '1'), ('1', 'g')):
            bk, bf = s.B[k], s.B[f]; Ck = bk['C']
            # pontual (unidade 1): |y(0)| <= Ck int KPi_k(t) OmA_f(t) dt
            y0 = Ck*trap(s.KP[k]*s.OmA[f], x)
            # perfil phi-unidades (rota 1 pontual): U(s) = Ck[KPi_k(s) + OmP_f(S) KPi_k(S-s) + nu_k^{-1/3} int KPi1_k(|s-t|) OmP_f(t) dt]
            conv = np.convolve(s.OmP[f], np.r_[s.KP1[k][::-1], s.KP1[k][1:]], mode='full')[len(x)-1:2*len(x)-1]*ds*MARG
            OmPS = min(tp.F('OmPsi', bf['nu']**(-1/3)*s.S, rate_tail=K['a']), np.exp(-K['a']*bf['nu']**(-1/3)*s.S))
            KPS = tp.F('KPi', bk['nu']**(-1/3)*(s.S - x), rate_tail=K['at'])
            U = Ck*(s.KP[k] + OmPS*KPS + bk['nu']**(-1/3)*conv)
            # perfis de correlacao para a rota 3
            PhiQ = np.correlate(s.OmP[f], s.Q[k], mode='full')[len(x)-1:]*ds*MARG          # int_u^S OmP_f(t) Q_k(t-u) dt
            PhiQs = np.convolve(s.OmP[f], s.Qs[k], mode='full')[:len(x)]*ds*MARG             # int_0^u OmP_f(t) Qs_k(u-t) dt
            s.pair[k+f] = dict(k=k, f=f, y0=y0, U=U, PhiQ=PhiQ, PhiQs=PhiQs, OmPS=OmPS,
                               JG=Jc(bk['nu'], bf['nu'], '1G'), JH=Jc(bk['nu'], bf['nu'], '1'), JsH=Jc(bk['nu'], bf['nu'], 's'))
    # ---- normas do nucleo via KPi
    def IK(s, k, rho):            # sup_s e^{rho s} int KPi_k(|s-t|) e^{-rho t} dt  <=  int_0^inf KPi_k(l)(e^{rho l}+e^{-rho l})
        x = s.x; return trap(s.KP[k]*(np.exp(rho*x) + np.exp(-rho*x)), x)
    def K0(s, k, rho):            # int_0^inf KPi_k(t) e^{-rho t} dt   ( (K h)(0) <= C K0 |h|_rho )
        return trap(s.KP[k]*np.exp(-rho*s.x), s.x)
def wvec(x, beta, lam, S, mode):
    if mode == 'm': return np.exp(-beta*x), beta*np.exp(-beta*x)
    a = np.exp(-beta*(S - x)); b = lam*np.exp(-beta*(S + x))
    return a - b, beta*(a + b)
def route3(pr, key, beta, lam, mode, h0, hS):
    """bound de |int w y|/|phi_f| pela rota 3 (Lema F4).  Valida para mode 'm' (qualquer beta) e mode 'd' com lam=1."""
    p = pr.pair[key]; k, f = p['k'], p['f']; bk, bf = pr.B[k], pr.B[f]; x = pr.x; S = pr.S
    if mode == 'd' and lam != 1: return np.inf
    w, wp = wvec(x, beta, lam, S, mode)
    wS = float(wvec(np.array([S]), beta, lam, S, mode)[0][0])
    if mode == 'm':
        main = JGb(bk['nu'], bf['nu'], beta) + trap(pr.OmP[f]*(1 - np.exp(-beta*x))*pr.Q[k], x)
    else:
        main = min(2*beta*p['JsH'], wS*p['JH'])
    QsS = tp.F('Qs', bk['nu']**(-1/3)*(S - x), rate_tail=K['at'])
    tS = wS*trap(pr.OmP[f]*QsS, x)
    tw = trap(wp*p['PhiQ'], x) + trap(wp*p['PhiQs'], x)
    return h0 + p['OmPS']*hS + bk['C']*(main + tS + tw)
def crit(N, J, G, sig, beta, lam, th_rho=0.5, variant='i', pr=None, coupling=1.0, override=None):
    if 2*sig**(1/3) > 90: return C2.crit(N, J, G, sig, beta, lam, th_rho, cone=(variant in ('i', 'ir')), coupling=coupling)
    if pr is None: pr = Prof(N, J, G, sig)
    P = pr.P; g, nu1, c = pr.g, pr.nu1, pr.c; k12, k21 = coupling*P['k12'], coupling*P['k21']; S = pr.S; x = pr.x
    a, at, ct, mF = K['a'], K['at'], K['ct'], K['mF']
    B1, Bg = pr.B['1'], pr.B['g']; r1, rg, a1, ag = B1['r'], Bg['r'], B1['at'], Bg['at']
    rho = th_rho*a1
    N1, Ng = B1['C']*pr.IK('1', rho), Bg['C']*pr.IK('g', rho)
    tau2 = k12*k21*N1*Ng
    out = dict(S=S, al=beta*S/2, tau2=tau2, lam=lam, th_rho=th_rho, cone=variant)
    if tau2 >= 1: return dict(out, ok=False, why='tau2')
    e2 = -np.expm1(-beta*S) if lam == 1 else 1 - lam*np.exp(-beta*S)
    K1n = B1['C']*2*trap(pr.KP['1'], x); Kgn = Bg['C']*2*trap(pr.KP['g'], x)   # v4 (Arb): sem fator nu^{1/3} extra
    tau1 = k12*k21*K1n*Kgn
    if tau1 >= 1: return dict(out, ok=False, why='tau1')
    pu, pv = pr.pair['g1'], pr.pair['1g']
    # ---- coluna U_1: pontual
    u1 = Bg['C']*pr.IK('g', rho); B_1 = k21*u1/(1 - tau2)
    B20 = k12*k21*Bg['C']*pr.K0('g', rho)*N1*B_1; A0 = k12*B1['C']*pr.K0('1', rho)*B_1
    P1 = A0 + (2 - c)*(k21*pu['y0'] + B20)
    Ub = (np.exp(rho*x)*pu['U']).max()
    if S > pr.Sc: Ub = max(Ub, C2_Ub(pr, 'g1', rho))
    Bp = k21*Ub/(1 - tau2); Ap = k12*N1*Bp
    Q1 = min(np.exp(-rho*S)*(Ap + (2 - c)*Bp), (2 - c)*k21*pu['U'][-1] + np.exp(-rho*S)*(Ap + (2 - c)*tau2*Bp)) if S <= pr.Sc else np.exp(-rho*S)*(Ap + (2 - c)*Bp)
    L1u = trap(pu['U'], x); L1B = k21*L1u/(1 - tau1); L1A = k12*K1n*L1B
    # ---- coluna U_g
    v1 = B1['C']*pr.IK('1', rho); A_g = k12*v1/(1 - tau2)
    A20 = k12*k21*B1['C']*pr.K0('1', rho)*Ng*A_g; B0 = k21*Bg['C']*pr.K0('g', rho)*A_g
    Pg = k12*pv['y0'] + A20 + (2 - c)*B0
    Vb = (np.exp(rho*x)*pv['U']).max()
    Agp = k12*Vb/(1 - tau2); Bgp = k21*Ng*Agp
    Qg = min(np.exp(-rho*S)*(Agp + (2 - c)*Bgp), k12*pv['U'][-1] + np.exp(-rho*S)*(tau2*Agp + (2 - c)*Bgp)) if S <= pr.Sc else np.exp(-rho*S)*(Agp + (2 - c)*Bgp)
    L1v = trap(pv['U'], x); L1Ag = k12*L1v/(1 - tau1); L1Bg = k21*Kgn*L1Ag
    # ---- momentos: min(rota1 nova, rotas 1/2 de crit2, rota3)
    Mo1, Mog = {}, {}
    for mode in ('m', 'd'):
        w, _ = wvec(x, beta, lam, S, mode)
        Wt = C2.Weight(mode, beta, lam, S)
        o1 = C2.col_bounds(Bg, B1, Wt, K); o2 = C2.col_bounds(B1, Bg, Wt, K)
        # h(0), h(S) da rota 2 antiga (constantes uniformes) -- reaproveitadas
        h0u, hSu = hbounds(Bg, Wt); h0v, hSv = hbounds(B1, Wt)
        MU = min(trap(w*pu['U'], x), o1[0], o1[1], route3(pr, 'g1', beta, lam, mode, h0u, hSu), route4(pr, 'g1', beta, lam, mode, N, J, G))
        MV = min(trap(w*pv['U'], x), o2[0], o2[1], route3(pr, '1g', beta, lam, mode, h0v, hSv), route4(pr, '1g', beta, lam, mode, N, J, G))
        if override is not None:           # projecao crit3+ (valores exatos maximizados na reta; [NUMERICO])
            MU = min(MU, override[0 if mode == 'm' else 1]); MV = min(MV, override[2 if mode == 'm' else 3])
        wint_rho = trap(w*np.exp(-rho*x), x); wmax = w.max()
        Mo1[mode] = c*k21*MU + min((Ap + c*k12*k21*Ng*N1*Bp)*wint_rho, wmax*(L1A + c*tau1*L1B))
        Mog[mode] = k12*MV + min((tau2*Agp + c*Bgp)*wint_rho, wmax*(tau1*L1Ag + c*L1Bg))
    if P1 >= 1: return dict(out, ok=False, why='P1')
    # ---- desacoplado com OmPsi (Lema F2)
    def blow(nu):    return 1 - beta*trap(np.exp(-beta*x)*pr.OmP[nm(nu)], x)
    def tail(nu):    return 2*pr.pair['g1' if nu == '1' else '1g']['OmPS'] if False else 2*min(tp.F('OmPsi', pr.B[nu]['nu']**(-1/3)*S, rate_tail=a), np.exp(-a*pr.B[nu]['nu']**(-1/3)*S))
    def Bm(nu):      return 1 + np.exp(-beta*S)*min(tp.F('OmPsi', pr.B[nu]['nu']**(-1/3)*S, rate_tail=a), np.exp(-a*pr.B[nu]['nu']**(-1/3)*S)) + beta*trap(np.exp(-beta*x)*pr.OmP[nu], x)
    def Dl(nu):      # e^{-2al}[(1-lam) + (e^{bS}-lam e^{-bS}) OmPsi(S') + beta int (e^{bs}+lam e^{-bs}) OmPsi]
        OS = min(tp.F('OmPsi', pr.B[nu]['nu']**(-1/3)*S, rate_tail=a), np.exp(-a*pr.B[nu]['nu']**(-1/3)*S))
        return (1 - lam)*np.exp(-beta*S) + (1 - lam*np.exp(-2*beta*S))*OS + beta*trap((np.exp(-beta*(S - x)) + lam*np.exp(-beta*(S + x)))*pr.OmP[nu], x)
    def nm(nu): return nu
    def blow2(nu):   # max(cota exponencial exata, cota de dois pontos)
        r = pr.B[nu]['r']; return max(r/(beta + r), 1 - beta*np.trapezoid(np.exp(-beta*x)*pr.OmP[nu], x)*MSUB - pr.OmP[nu][-1]*np.exp(-beta*x[-1]))
    bl1 = blow2('1'); blg = blow2('g')
    epx = ((2 - c)*P1 + Pg)/((1 - P1)*(2 - c))
    if epx >= 1: return dict(out, ok=False, why='epx')
    ph1 = nu1**(1/3)/mF; xb = (2 - c)*(1 + epx)
    D1, Dg = Dl('1'), Dl('g')
    rel = variant.endswith('r'); variant = variant[:-1] if rel else variant
    if rel:      # Lema F8: cone relativo entre as duas familias (vale para todo beta; util para beta pequeno)
        d = (nu1 - g)*beta**2 + (P['kap']*N/(1 + N - g) - P['kap']*(1 - g)/(1 + N - g))*sig**(-2/3)
        e1, eg = 1 - bl1, 1 - blg
        if max(e1, eg) >= 1: return dict(out, ok=False, why='rel')
        # v5: so' o ramo de tabela (enunciado do Lema F8: Theta*(beta) calculado com b_nu, SEM termos arcsen(1-b_low));
        # o antigo ramo de 'fallback' (Theta*(0)+arcsen(1-b_low)) foi removido -- sem tabela, a variante falha ('noT8').
        threl = np.inf
        tb = rc_table(N, J, G)
        if tb is not None and beta <= tb[-1][0]:
            bsx = [r[0] for r in tb]; k = int(np.searchsorted(bsx, beta + 1e-12, side='right') - 1)
            r0 = tb[k]; r1 = tb[min(k + 1, len(tb) - 1)]; hstep = max(r1[0] - r0[0], 1e-12)
            mx = lambda j: 1.05*max(r0[j], r1[j])
            Lth = mx(2)*2*nu1**(2/3)*(beta + hstep) + mx(3)*nu1**(1/3) + mx(4)*2*g**(2/3)*(beta + hstep) + mx(5)*g**(1/3)
            ds = abs(P['kap']*N/(1 + N - g) - P['kap']*(1 - g)/(1 + N - g))*sig**(-2/3)
            th2 = 1.01*r0[1] + Lth*(beta - r0[0]) + g**(-1/3)*ds*mx(4) + (np.arcsin(epx) if variant == 'i' else 0.0)
            threl = th2
        if not np.isfinite(threl): return dict(out, ok=False, why='noT8')
        if threl >= np.pi: return dict(out, ok=False, why='rel')
        crel = np.cos(threl/2)
        kmg = crel*blg - tail('g')
    else:
        kmg = ct*blg - tail('g')
    if variant == 'i':
        if rel:
            km1 = (1 - epx)*crel*bl1 - (1 + epx)*tail('1')
        else:
            thp = np.arccos(ct) + np.arcsin(epx)
            if thp >= np.pi/2: return dict(out, ok=False, why='cone')
            km1 = (1 - epx)*np.cos(thp)*bl1 - (1 + epx)*tail('1')
        if km1 <= 0 or kmg <= 0: return dict(out, ok=False, why='km')
        th_max = 1 - max((1 + epx)*D1/(e2*km1), Dg/(e2*kmg))
        eW = max(xb*Mo1['m']/((2 - c)*km1), Mog['m']/(c*kmg))
        eD = max(xb*Mo1['d']/((2 - c)*e2*km1), Mog['d']/(c*e2*kmg))
    else:
        km1 = (crel if rel else ct)*bl1 - tail('1')
        if km1 <= 0 or kmg <= 0: return dict(out, ok=False, why='km')
        th_max = 1 - max(D1/(e2*km1), Dg/(e2*kmg))
        dx = (2 - c)*epx
        eW = max((dx*Bm('1') + xb*Mo1['m'])/((2 - c)*km1), Mog['m']/(c*kmg))
        eD = max((dx*D1 + xb*Mo1['d'])/((2 - c)*e2*km1), Mog['d']/(c*e2*kmg))
    thK = th_max - eW - eD
    out.update(P1=P1, Pg=Pg, eW=eW, eD=eD, th_max=th_max, thK=thK, Mo1=Mo1['m'], Mog=Mog['m'])
    if thK <= 0 or eW >= 1: return dict(out, ok=False, why='thK')
    Sp1, Spg = 2*(sig/nu1)**(1/3), 2*(sig/g)**(1/3)
    P_1, P_g = Pfine(Sp1), Pfine(Spg)
    a4 = xb*(nu1**(-1/3)*P_1 + Q1)/(2 - c); b4 = ((2 - c)*g**(-1/3)*P_g + Qg)/c
    E4 = max(a4/km1, b4/kmg); e3 = ph1*(nu1**(-1/3)*P_1 + Q1)
    y_m = Bm('1') + Mo1['m']; y_l = (D1 + Mo1['d'])/e2
    dR = (1 - th_max + eD)/(1 - eW)
    Y = y_m*dR + y_l; Z2 = 2*y_m + lam*y_l + lam*y_m*dR + dR*y_l; Z3 = 2*y_m**2 + 2*lam*y_m*y_l + y_l**2
    q = 1 - P1
    leak = 2*E4*ph1*Y/(q*thK) + e3**2/q**2 + 2*E4*e3*ph1*Z2/(q**2*thK) + E4**2*ph1**2*Z3/(q**2*thK*(1 - eW))
    out.update(E4=E4, e3=e3, leak=leak, ok=bool(leak < 1))
    return out
def hbounds(kb, Wt):
    """|h(0)|, |h(S)| (Lema C3(c), constantes uniformes)"""
    S = Wt.S; t = np.array([0.0, S]); nu, r, at = kb['nu'], kb['r'], kb['at']
    GR, GX, TR, TX = Wt.Gam_R(t, r), Wt.Gam_X(t, r), Wt.The_R(t, at), Wt.The_X(t, at)
    h0 = kb['C']*min(nu**(1/3)*K['Pipsi']*GR[0], K['Pi']*TR[0]); hS = kb['C']*min(nu**(1/3)*K['Pipsi']*GX[-1], K['Pi']*TX[-1])
    return h0, hS
def C2_Ub(pr, key, rho): return np.inf

# ---------------- rota 4: constantes conjuntas de 2 variaveis (Lema F6), tabela por terna ----------------
_MT = {}
def Mtab(N, J, G):
    key = (N, J, G)
    if key not in _MT:
        f = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Mtab_%g_%g_%g.json' % (N, J, G))
        _MT[key] = json.load(open(f)) if os.path.exists(f) else None
    return _MT[key]
def corrS(pr, key, wabs_in, wabs_out_kind):
    """|(funcional em [0,S]) - (funcional S=inf)| / |phi_f|  (Lema F6, parte (c)); wabs_in: |w| na grade x; fora: 1 ('one') ou s ('s')"""
    p = pr.pair[key]; k, f = p['k'], p['f']; bk, bf = pr.B[k], pr.B[f]; x = pr.x; S = pr.S; C = bk['C']
    kap = bk['nu']**(-1/3); a = K['a']; at = K['at']
    OmS = min(tp.F('OmPsi', bf['nu']**(-1/3)*S, rate_tail=a), np.exp(-a*bf['nu']**(-1/3)*S))
    u = np.arange(0, 60, 0.01)
    KP1 = np.minimum(tp.F('KPi1', kap*u, rate_tail=at), K['Pi1']*np.exp(-at*kap*u)); KP = np.minimum(tp.F('KPi', kap*u, rate_tail=at), K['Pi']*np.exp(-at*kap*u))
    tailKP1 = np.concatenate([np.cumsum((KP1[::-1][1:] + KP1[::-1][:-1])/2*0.01)[::-1], [0]])   # int_u^inf KP1
    ins = trap(wabs_in*(np.interp(S - x, u, KP)*OmS + kap*OmS*np.interp(S - x, u, tailKP1)), x)
    # fora: s em [S, S+40]
    so = S + np.arange(0, 40, 0.02); Omf = lambda t: np.minimum(tp.F('OmPsi', bf['nu']**(-1/3)*t, rate_tail=a), np.exp(-a*bf['nu']**(-1/3)*t))
    I1 = np.trapezoid(KP1, u); IOm = np.trapezoid(Omf(u), u)
    Uinf = np.interp(kap*so, kap*u, KP) + kap*(2*I1*Omf(so/2) + np.interp(kap*so/2, kap*u, KP1)*IOm)
    wout = np.ones_like(so) if wabs_out_kind == 'one' else so
    out = np.trapezoid(wout*Uinf, so)*MARG + (Uinf[-1]*10)        # cauda alem de S+40 (grosseira)
    return C*(ins + out)
def Uinf_mom(pr, key, bk, L=150.0, ds=0.01):
    """int_0^inf s e^{-bk s} U_inf(s) ds, U_inf = perfil da rota 1 com S=inf (sem termo OmP_f(S)):
       U_inf(s) = C_k[KPi_k(s) + nu_k^{-1/3} int_0^inf KPi1_k(|s-t|) OmP_f(t) dt]  (cota de |y_inf|/|phi_f|)"""
    p = pr.pair[key]
    if 'Uinf' not in p:
        k, f = p['k'], p['f']; bk_, bf = pr.B[k], pr.B[f]
        xs = np.arange(0, L + ds/2, ds)
        lk, lf = bk_['nu']**(-1/3)*xs, bf['nu']**(-1/3)*xs
        KPk = np.minimum(tp.F('KPi', lk, rate_tail=K['at']), K['Pi']*np.exp(-K['at']*lk))
        KP1k = np.minimum(tp.F('KPi1', lk, rate_tail=K['at']), K['Pi1']*np.exp(-K['at']*lk))
        Omf = np.minimum(tp.F('OmPsi', lf, rate_tail=K['a']), np.exp(-K['a']*lf))
        n = len(xs)
        conv = np.convolve(Omf, np.r_[KP1k[::-1], KP1k[1:]], mode='full')[n-1:2*n-1]*ds*MARG
        # cauda t>L da convolucao: <= Pi1 int_L^inf OmP_f
        conv += K['Pi1']*Omf[-1]/(K['a']*bf['nu']**(-1/3))
        U = bk_['C']*(KPk + bk_['nu']**(-1/3)*conv)
        p['Uinf'] = (xs, U)
    xs, U = p['Uinf']
    r = K['a']*min(pr.B['1']['nu'], pr.B['g']['nu'])**(-1/3)*0.5     # taxa conservadora para a cauda s>L
    tail = U[-1]*np.exp(0)*(xs[-1]/r + 1/r**2)
    return trap(xs*np.exp(-bk*xs)*U, xs) + tail
def route4(pr, key, beta, lam, mode, N, J, G):
    T4 = Mtab(N, J, G)
    if T4 is None: return np.inf
    idx = 0 if key == 'g1' else 1; p = pr.pair[key]; kb = pr.B[p['k']]; x = pr.x; S = pr.S
    if mode == 'm':
        nodes = sorted(float(b) for b in T4 if b != 's')
        if beta > nodes[-1]: return np.inf
        bk = max(b for b in nodes if b <= beta)
        M = T4['%g' % bk][idx]
        lip = (beta - bk)*Uinf_mom(pr, key, bk)          # v5: Lipschitz com U_inf em [0,inf) (correcao pedida pelo Arb v3/v4)
        return kb['C']*M + lip + corrS(pr, key, np.exp(-beta*x), 'one')
    if mode == 'd' and lam == 1 and beta*S < 30:
        Ms = T4['s'][idx]
        pref = np.exp(-beta*S)
        rem = trap((beta*x)**3/3*0.5*(np.exp(-beta*(S - x)) + np.exp(-beta*(S + x)))*p['U'], x)
        return pref*2*beta*(kb['C']*Ms + corrS(pr, key, x, 's')) + rem
    return np.inf

_RC = {}
PHIP = 0.8784*1.01            # sup_{H_delta}|F'/F| (numerico; a certificar)
def theta_star(nu1, g):
    k = (round(nu1, 6), round(g, 6))
    if k not in _RC:
        f = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'relcone_cache.json')
        dd = json.load(open(f)) if os.path.exists(f) else {}
        kk = '%g|%g' % k
        if kk not in dd:
            import relcone; dd[kk] = relcone.theta_star((nu1/g)**(1/3))[0]*1.01; json.dump(dd, open(f, 'w'))
        _RC[k] = dd[kk]
    return _RC[k]
_RT = {}
def rc_table(N, J, G):
    k = (N, J, G)
    if k not in _RT:
        f = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'relcone2_%g_%g_%g.json' % (N, J, G))
        _RT[k] = json.load(open(f)) if os.path.exists(f) else None
    return _RT[k]
_PT = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ptable.npz'))
def Pfine(Sp):
    """P(S') <= min(tabela fina (no' a esquerda), tabela/forma exponencial do G0)"""
    v = C2.Pbound(Sp); Sg, Pg = _PT['Sp'], _PT['P']
    if Sp >= Sg[0]: v = min(v, Pg[min(int(np.floor((Sp - Sg[0])/(Sg[1] - Sg[0]) + 1e-12)), len(Sg) - 1)])
    return v
_PR = {}
def best(N, J, G, sig, beta, lams=(0.0, 1.0), ths=(0.35, 0.5, 0.65, 0.8), variants=('i', 'ii', 'ir', 'iir'), coupling=1.0, overrides=None):
    if 2*sig**(1/3) > 90:
        return C2.best(N, J, G, sig, beta, coupling=coupling)
    key = (N, J, G, sig)
    if key not in _PR:
        if len(_PR) > 4: _PR.clear()
        _PR[key] = Prof(N, J, G, sig)
    pr = _PR[key]; bv = None
    for lam in lams:
        for th in ths:
            for v in variants:
                o = crit(N, J, G, sig, beta, lam, th, v, pr, coupling, None if overrides is None else overrides[lam])
                sc = (1 - o['leak'])*o['thK'] if o.get('ok') else -1
                if bv is None or sc > bv[0]: bv = (sc, o)
    return bv
if __name__ == '__main__':
    import sys, time
    N, J, G = map(float, sys.argv[1:4])
    for sig in [300, 600, 1000, 2000]:
        for beta in [1e-8, 0.05, 0.3, 1.0, 3.0]:
            t0 = time.time(); sc, o = best(N, J, G, sig, beta)
            print('sig=%5.0f beta=%6.1e ok=%-5s lam=%g v=%s P1=%.3f eW=%.3f eD=%.3f thmax=%.3f thK=%.3f leak=%.2g %s (%.1fs)' % (
                sig, beta, o.get('ok'), o['lam'], o.get('cone'), o.get('P1', np.nan), o.get('eW', np.nan), o.get('eD', np.nan),
                o.get('th_max', np.nan), o.get('thK', np.nan), o.get('leak', np.nan), o.get('why', ''), time.time() - t0), flush=True)
