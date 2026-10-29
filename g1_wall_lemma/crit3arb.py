"""Rigorous BOX version of crit3 (G1'' cotas finas) -- 'crit3a' = crit3 WITHOUT routes 3/4 (pluggable later).
Ingredients (all certified): T1 OmA, T2 OmPsi (tp_cert_0.05.json), T3 KPi, T4 KPi1, T5 Q, T6 Q* (T36_cert_0.05.json),
T7 fine P (ptable_0.05.json), and the crit2 constants.  Tables are NONINCREASING step majorants on l = 0:0.025:16,
evaluated at the left node <= argument; beyond 16: T(16) e^{-rate (l-16)} (submultiplicativity, Lemma F1).
Integrals in s: upper sums over cells of width ds (left-node values of nonincreasing factors), floating point sums of
NONNEGATIVE terms with a relative safety factor (1+1e-9) >> accumulated rounding; scalar combination in Arb (arbcrit).
Box: sigma in [s1,s2] (S in [S1,S2]), beta in [b1,b2].
Interface for driver.py (G1_CRIT=crit3arb):  CRIT(P, s1, s2, b1, b2, choice, coupling) ; cands(NJG, s, b, coupling) ; params."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, os, json, math
import numpy as np
sys.path.insert(0, _os.path.join(_R, 'g1_wall_lemma')); sys.path.insert(0, _os.path.join(_R, 'g1_wall_lemma/heuristic'))
from arbprim import *
import arbcrit as AC
from arbcrit import K, Weight, WeightRefl, route_integrals, cells as tcells
def params(N, J, G):
    P = AC.params(N, J, G); P['key'] = '%g_%g_%g' % (float(N), float(J), float(G)); return P
_JD = {}
def Jdata(P):
    k = P.get('key')
    if k not in _JD:
        f = CD + 'J_cert_%s.json' % k
        _JD[k] = json.load(open(f)) if (k and os.path.exists(f)) else None
    return _JD[k]
UP = 1 + 1e-9; DN = 1 - 1e-9
_T8 = {}
def Theta_bar(P, b1, b2, s1):
    """certified T8 (Lemma F8): sup of the relative-cone angle over beta in [b1,b2], valid for sigma >= sigma_min.
    Several tables may exist (T8_cert_<key>.json, T8_cert_<key>_s<smin>.json): each one valid with sigma_min <= s1
    gives an upper bound; the minimum is used."""
    import glob as _g
    k = P.get('key')
    if k not in _T8:
        fs = sorted(_g.glob(_os.path.join(_CERT, 'g0_constants', 'T8_cert_%s.json' % k)) + _g.glob(_os.path.join(_CERT, 'g0_constants', 'T8_cert_%s_s*.json' % k)))
        _T8[k] = [json.load(open(f)) for f in fs]
    best = None
    for d in _T8[k]:
        if s1 < d['sigma_min']: continue
        v = None; okc = True
        lo_cov = min(a for a, b, th in d['Theta_edge']); hi_cov = max(b for a, b, th in d['Theta_edge'])
        if b1 < lo_cov or b2 > hi_cov: continue
        for a, b, th in d['Theta_edge']:
            if a <= b2 and b >= b1:
                if not (th < 10): okc = False; break
                v = th if v is None else max(v, th)
        if okc and v is not None: best = v if best is None else min(best, v)
    return best
CD = _os.path.join(_CERT, 'g0_constants') + _os.sep
_t12 = json.load(open(CD + 'tp_cert_0.05.json')); _t36 = json.load(open(CD + 'T36_cert_0.05.json'))
DL = 0.025
TAB = {'OmA': np.array(_t12['OmA']), 'OmPsi': np.array(_t12['OmPsi']), 'KPi': np.array(_t36['KPi']),
       'KPi1': np.array(_t36['KPi1']), 'Q': np.array(_t36['Q']), 'Qs': np.array(_t36['Qs'])}
af, atf = 0.44452699, 0.35346270                      # lower bounds a, at
RATE = {'OmA': atf, 'OmPsi': af, 'KPi': atf, 'KPi1': atf, 'Q': atf, 'Qs': atf}
UNIF = {'OmA': (1.0, atf), 'OmPsi': (1.0, af), 'KPi': (float(K['Pi'].upper()), atf), 'KPi1': (float(K['Pi1'].upper()), atf),
        'Q': (float(K['Pi1psi'].upper()), af), 'Qs': (float(K['Pi1psi'].upper()), af)}
def T(name, l):
    """upper bound of table 'name' at arguments >= l (array, l >= 0)"""
    l = np.maximum(np.asarray(l, float), 0.0)
    v = TAB[name]; n = len(v)
    idx = np.floor(l*(1 - 1e-12)/DL).astype(np.int64)
    out = np.where(idx < n, v[np.minimum(idx, n - 1)], v[-1]*np.exp(-RATE[name]*(l - (n - 1)*DL)*DN))
    A0, r = UNIF[name]
    return np.minimum(out, A0*np.exp(-r*l*DN)*UP)*UP
CDJ = _os.path.join(_CERT, 'g0_constants') + _os.sep
_ptab = sorted((float(a), float(b)) for a, b in json.load(open(CD + 'ptable_0.05.json')))
def Pfine(Sp_lo):
    v = float(AC.Pbound(max(Sp_lo, 2.0)).upper()) if Sp_lo >= 2 else math.inf
    for s, p in _ptab:
        if Sp_lo >= s: v = min(v, p)
    return v

class Prof:
    """sigma-box profiles (beta independent).  s-cells [i ds, (i+1) ds], i < n, covering [0, S2]."""
    def __init__(s, P, s1, s2, ds=0.01):
        s.P = P; s.s1, s.s2 = s1, s2
        s.g = float(P['g'].mid()); s.nu1 = float(P['nu1'].mid()); s.c = P['c']
        s.S1 = float((2*arb(s1)**(ONE/3)).lower()); s.S2 = float((2*arb(s2)**(ONE/3)).upper())
        s.ds = ds; s.n = int(math.ceil(s.S2/ds)); s.lo = np.arange(s.n)*ds; s.hi = s.lo + ds
        s.B = {}
        for nm, nuA in (('1', P['nu1']), ('g', P['g'])):
            nu = float(nuA.mid())
            s.B[nm] = dict(nu=nu, nuA=nuA, kap=float((nuA**(-ONE/3)).lower()), kapU=float((nuA**(-ONE/3)).upper()),
                           C=float((2*PI*nuA**(-2*ONE/3)*arb(s1)**(-2*ONE/3)).upper()),
                           r=float((K['a']*nuA**(-ONE/3)).lower()), at=float((K['at']*nuA**(-ONE/3)).lower()))
        lo = s.lo
        s.KP = {k: T('KPi', b['kap']*lo) for k, b in s.B.items()}
        s.KP1 = {k: T('KPi1', b['kap']*lo) for k, b in s.B.items()}
        s.OmA = {k: T('OmA', b['kap']*lo) for k, b in s.B.items()}
        s.OmP = {k: T('OmPsi', b['kap']*lo) for k, b in s.B.items()}
        s.OmPS = {k: float(T('OmPsi', np.array([b['kap']*s.S1]))[0]) for k, b in s.B.items()}
        s.pair = {}
        for k, f in (('g', '1'), ('1', 'g')):
            bk, bf = s.B[k], s.B[f]; Ck = bk['C']
            y0 = Ck*np.sum(s.KP[k]*s.OmA[f])*ds*UP
            # conv(i) = int_0^inf KPi1(kap_k|s-t|) OmPsi_f(t) dt, s in cell i: |s-t| >= (|i-j|-1)_+ ds
            m = np.arange(s.n)
            kern = T('KPi1', bk['kap']*np.maximum(m - 1, 0)*ds)
            full = np.convolve(s.OmP[f], np.r_[kern[::-1], kern[1:]])          # sum_j Om_j kern[|i-j|]
            conv = full[s.n - 1: 2*s.n - 1]*ds
            # tail t >= n ds: OmPsi_f(t) <= OmPsi_f(n ds) e^{-a kap_f (t - n ds)}, KPi1 <= Pi1
            X = s.n*ds; tailc = UNIF['KPi1'][0]*float(T('OmPsi', np.array([bf['kap']*X]))[0])/(af*bf['kap'])
            conv = (conv + tailc)*UP
            KPS = T('KPi', bk['kap']*np.maximum(s.S1 - s.hi, 0.0))
            U = Ck*(s.KP[k] + s.OmPS[f]*KPS + bk['kapU']*conv)*UP
            s.pair[k + f] = dict(k=k, f=f, y0=y0, U=U, US=U[-1])
    def IK(s, k, rho):   # int_0^S KPi(kap s)(e^{rho s}+e^{-rho s}) ds
        lo, ds = s.lo, s.ds
        return np.sum(s.KP[k]*(np.exp(rho*(lo + ds)) + np.exp(-rho*lo)))*ds*UP
    def K0(s, k, rho):
        return np.sum(s.KP[k]*np.exp(-rho*s.lo))*s.ds*UP
    def L1K(s, k):        # int_0^S KPi(kap s) ds
        return np.sum(s.KP[k])*s.ds*UP

_PROF = {}
def prof(P, s1, s2):
    key = (id(P), s1, s2)
    if key not in _PROF:
        if len(_PROF) > 6: _PROF.clear()
        _PROF[key] = Prof(P, s1, s2)
    return _PROF[key]

def wup(pr, b1, b2, lam, mode):
    """upper bound of the weight on each s-cell (mode 'm': e^{-beta s}; 'd': lambda=1 normalised w_hat, or lambda=0 e^{-beta(S-s)})"""
    lo, hi = pr.lo, pr.hi
    if mode == 'm': return np.exp(-b1*lo)
    if lam == 1:   # w_hat(s) = 2 sinh(beta s)/(e^{beta S}-1) <= min(Kn(b1,S1) sinh(b2 s)/b2 ..., 1+e^{-b1 S1}); Kn = 2/(S psi(beta S))
        # w_hat is NONINCREASING in beta and in S for fixed s < S (d/dbeta log w_hat = s coth(beta s) - S/(1-e^{-beta S}) <= 0),
        # so w_hat <= w_hat(b1, S1; s) = Kn(b1,S1) sinh(b1 s)/b1
        Kn = float((2/(arb(pr.S1)*psi(arb(b1)*arb(pr.S1)))).upper())
        sh = hi*np.where(b1*hi > 1e-8, np.sinh(np.minimum(b1*hi, 700))/np.maximum(b1*hi, 1e-300), 1.0)
        return np.minimum(Kn*sh*UP, 1 + math.exp(-b1*pr.S1))
    return None    # lambda=0 handled by reflection
def int_wU(pr, U, b1, b2, lam, mode):
    if mode == 'd' and lam != 1:      # int_0^S e^{-beta(S-s)} U(s) ds = int_0^S e^{-beta v} U(S-v) dv
        v = pr.lo; n = pr.n
        # U(S - v) for v in cell, S in [S1,S2]: s in [S1 - v_hi, S2 - v_lo] -> max of U over those cells
        a = np.floor(np.maximum(pr.S1 - (v + pr.ds), 0)/pr.ds).astype(int); b = np.minimum(np.floor((pr.S2 - v)/pr.ds).astype(int), n - 1)
        cm = np.maximum.accumulate(U[::-1])[::-1]      # suffix max (U roughly decreasing): max over [a, n) >= max over [a,b]
        return np.sum(np.exp(-b1*v)*cm[np.minimum(a, n - 1)])*pr.ds*UP
    return np.sum(wup(pr, b1, b2, lam, mode)*U)*pr.ds*UP

def crit_box3(P, s1, s2, b1, b2, lam, th_rho, cone, coupling=1, route2=True, extra_routes=None):
    pr = prof(P, s1, s2)
    a, at, ct, mF, Pi = K['a'], K['at'], K['ct'], K['mF'], K['Pi']
    g, nu1, c = P['g'], P['nu1'], P['c']
    k12, k21 = float((coupling*P['k12']).upper()), float((coupling*P['k21']).upper())
    cf = float(c.upper()); c2 = float((2 - c).upper())
    B1, Bg = pr.B['1'], pr.B['g']
    rho = th_rho*B1['at']
    assert rho < Bg['at'] and rho < B1['at'] and rho <= B1['r'] and rho <= Bg['r']
    S1, S2 = pr.S1, pr.S2
    out = dict(lam=lam, th_rho=th_rho, cone=cone)
    def fail(w): return dict(out, margin=-math.inf, why=w)
    N1, Ng = B1['C']*pr.IK('1', rho), Bg['C']*pr.IK('g', rho)
    tau2 = k12*k21*N1*Ng*UP
    if tau2 >= 1: return fail('tau2')
    K1n, Kgn = 2*B1['C']*pr.L1K('1'), 2*Bg['C']*pr.L1K('g')          # ||K_nu||_{L1->L1} (no extra nu^{1/3})
    tau1 = k12*k21*K1n*Kgn*UP
    if tau1 >= 1: return fail('tau1')
    pu, pv = pr.pair['g1'], pr.pair['1g']
    eS = math.exp(-rho*S1)
    # column U_1
    u1 = Bg['C']*pr.IK('g', rho); B_1 = k21*u1/(1 - tau2)
    B20 = k12*k21*Bg['C']*pr.K0('g', rho)*N1*B_1; A0 = k12*B1['C']*pr.K0('1', rho)*B_1
    P1 = (A0 + c2*(k21*pu['y0'] + B20))*UP
    Ub = np.max(np.exp(rho*pr.hi)*pu['U'])*UP
    Bp = k21*Ub/(1 - tau2); Ap = k12*N1*Bp
    Q1 = min(eS*(Ap + c2*Bp), c2*k21*pu['US'] + eS*(Ap + c2*tau2*Bp))*UP
    L1u = np.sum(pu['U'])*pr.ds*UP; L1B = k21*L1u/(1 - tau1); L1A = k12*K1n*L1B
    # column U_g
    v1 = B1['C']*pr.IK('1', rho); A_g = k12*v1/(1 - tau2)
    A20 = k12*k21*B1['C']*pr.K0('1', rho)*Ng*A_g; B0 = k21*Bg['C']*pr.K0('g', rho)*A_g
    Pg = (k12*pv['y0'] + A20 + c2*B0)*UP
    Vb = np.max(np.exp(rho*pr.hi)*pv['U'])*UP
    Agp = k12*Vb/(1 - tau2); Bgp = k21*Ng*Agp
    Qg = min(eS*(Agp + c2*Bgp), k12*pv['US'] + eS*(tau2*Agp + c2*Bgp))*UP
    L1v = np.sum(pv['U'])*pr.ds*UP; L1Ag = k12*L1v/(1 - tau1); L1Bg = k21*Kgn*L1Ag
    if P1 >= 1: return fail('P1')
    # crit2 routes (Arb, uniform constants) for the same box
    Sa = iv(*[2*arb(x)**(ONE/3) for x in (s1, s2)]) if False else 2*iv(s1, s2)**(ONE/3)
    Ba = iv(b1, b2)
    s23 = iv(s1, s2)**(-2*ONE/3)
    KB1 = dict(nu=P['nu1'], r=a*P['nu1']**(-ONE/3), at=at*P['nu1']**(-ONE/3), C=2*PI*P['nu1']**(-2*ONE/3)*s23)
    KBg = dict(nu=P['g'], r=a*P['g']**(-ONE/3), at=at*P['g']**(-ONE/3), C=2*PI*P['g']**(-2*ONE/3)*s23)
    cl = tcells(Sa)
    Wm = Weight('m', Ba, Sa); Wd = Weight('d1n', Ba, Sa) if lam == 1 else WeightRefl(Ba, Sa)
    rhoA = arb(th_rho)*at*P['nu1']**(-ONE/3)
    Mo1, Mog = {}, {}
    for mode, Wt in (('m', Wm), ('d', Wd)):
        R1u, R2u = route_integrals(Wt, KBg, KB1['r'], route2 and b2 <= 1.0, cl)
        R1v, R2v = route_integrals(Wt, KB1, KBg['r'], route2 and b2 <= 1.0, cl)
        MU = min(int_wU(pr, pu['U'], b1, b2, lam, mode), float(R1u.upper()), float(R2u.upper()))
        MV = min(int_wU(pr, pv['U'], b1, b2, lam, mode), float(R1v.upper()), float(R2v.upper()))
        if extra_routes is not None:
            e = extra_routes(pr, P, b1, b2, lam, mode)
            MU, MV = min(MU, e[0]), min(MV, e[1])
        wr = float((Wt.ThR(Sa, rhoA) if getattr(Wt, 'refl', False) else Wt.ThR(ZERO, rhoA)).upper())
        wmax = float(Wt.wmax().upper())
        Mo1[mode] = (cf*k21*MU + min((Ap + cf*tau2*Bp)*wr, wmax*(L1A + cf*tau1*L1B)))*UP
        Mog[mode] = (k12*MV + min((tau2*Agp + cf*Bgp)*wr, wmax*(tau1*L1Ag + cf*L1Bg)))*UP
    # ---- Lemma 4.1 with OmPsi (F6): b_low, tails, Bm, Dl
    def Iom(k, b1_, b2_):   # beta int_0^inf e^{-beta s} OmPsi(kap s) ds  (upper), beta in [b1_,b2_]
        # = int_0^inf e^{-u} Om(kap u/beta) du is NONDECREASING in beta (Om step majorant nonincreasing) -> evaluate at b2_
        X = pr.n*pr.ds
        return (b2_*np.sum(np.exp(-b2_*pr.lo)*pr.OmP[k])*pr.ds + float(T('OmPsi', np.array([pr.B[k]['kap']*X]))[0])*math.exp(-b2_*X))*UP
    def blow(k):
        r = pr.B[k]['r']
        return max(r/(b2 + r)*DN, (1 - Iom(k, b1, b2))*DN)
    def tail(k): return 2*pr.OmPS[k]*UP
    def Bm(k): return (1 + math.exp(-b1*S1)*pr.OmPS[k] + Iom(k, b1, b2))*UP
    def Dl(k):   # e^{-2al} Delta(lam) / e2   (upper)
        kap = pr.B[k]['kap']; OS = pr.OmPS[k]
        if lam == 1:
            Kn = float((2/(arb(S1)*psi(arb(b1)*arb(S1)))).upper())
            ch = np.cosh(np.minimum(b2*pr.hi, 700))
            return ((1 + math.exp(-b1*S1))*OS + Kn*np.sum(ch*pr.OmP[k])*pr.ds)*UP
        v = pr.lo                                        # beta int_0^S e^{-beta v} OmPsi(kap (S - v)) dv
        Om_rev = T('OmPsi', kap*np.maximum(S1 - (v + pr.ds), 0.0))
        return (math.exp(-b1*S1) + OS + b2*np.sum(np.exp(-b1*v)*Om_rev)*pr.ds)*UP
    vals = dict(Mo1=Mo1, Mog=Mog, P1=P1, Pg=Pg, Q1=Q1, Qg=Qg, bl1=blow('1'), blg=blow('g'), tail1=tail('1'), tailg=tail('g'),
                D1=Dl('1'), Dg=Dl('g'), Bm1=Bm('1'),
                P_1=Pfine(float((Sa*P['nu1']**(-ONE/3)).lower())), P_g=Pfine(float((Sa*P['g']**(-ONE/3)).lower())))
    if isinstance(cone, str) and cone.endswith('r'):
        vals['Thbar'] = Theta_bar(P, b1, b2, s1)
        if vals['Thbar'] is None: return fail('noT8')
    return finish(P, vals, lam, cone, out)

def finish(P, v, lam, cone, out):
    """Props. C8-C9 in Arb from the (rigorous float) ingredient bounds."""
    A = lambda x: arb(x)
    ct, mF = K['ct'], K['mF']; g, nu1, c = P['g'], P['nu1'], P['c']; th = ONE/3
    def fail(w): return dict(out, margin=-math.inf, why=w)
    P1, Pg, Q1, Qg = A(v['P1']), A(v['Pg']), A(v['Q1']), A(v['Qg'])
    epx = ((2 - c)*P1 + Pg)/((1 - P1)*(2 - c))
    if not epx.upper() < 1: return fail('epx')
    ph1 = nu1**th/mF; xb = (2 - c)*(1 + epx)
    bl1, blg, t1, tg = A(v['bl1']), A(v['blg']), A(v['tail1']), A(v['tailg'])
    kmg = ct*blg - tg
    D1, Dg = A(v['D1']), A(v['Dg'])
    M1m, M1d, Mgm, Mgd = A(v['Mo1']['m']), A(v['Mo1']['d']), A(v['Mog']['m']), A(v['Mog']['d'])
    variant = cone if isinstance(cone, str) else ('i' if cone else 'ii')
    rel = variant.endswith('r'); base = variant[:-1] if rel else variant
    if rel:     # Lemma F8 (relative cone) with the certified T8 table
        e1, eg = 1 - bl1, 1 - blg
        if not (e1.upper() < 1 and eg.upper() < 1): return fail('rel')
        # theta_rel of Lemma F8 = Theta*(beta) (+ the sigma-shift, already inside the certified T8 cell values); the
        # beta-dependent table already contains b_nu, so no asin(1-b_low) terms (those belong to the beta=0 fallback)
        threl = A(v['Thbar']) + (epx.asin() if base == 'i' else ZERO)
        if not threl.upper() < PI.lower(): return fail('rel')
        crel = (threl/2).cos()
        kmg = crel*blg - tg
        if base == 'i':
            km1 = (1 - epx)*crel*bl1 - (1 + epx)*t1
        else:
            km1 = crel*bl1 - t1
        if not (km1.lower() > 0 and kmg.lower() > 0): return fail('km')
        if base == 'i':
            th_max = 1 - amax((1 + epx)*D1/km1, Dg/kmg)
            eW = amax(xb*M1m/((2 - c)*km1), Mgm/(c*kmg))
            eD = amax(xb*M1d/((2 - c)*km1), Mgd/(c*kmg))
        else:
            th_max = 1 - amax(D1/km1, Dg/kmg)
            dx = (2 - c)*epx
            eW = amax((dx*A(v['Bm1']) + xb*M1m)/((2 - c)*km1), Mgm/(c*kmg))
            eD = amax((dx*D1 + xb*M1d)/((2 - c)*km1), Mgd/(c*kmg))
    elif base == 'i':
        thp = ct.acos() + epx.asin()
        if not thp.upper() < (PI/2).lower(): return fail('cone')
        km1 = (1 - epx)*thp.cos()*bl1 - (1 + epx)*t1
        if not (km1.lower() > 0 and kmg.lower() > 0): return fail('km')
        th_max = 1 - amax((1 + epx)*D1/km1, Dg/kmg)
        eW = amax(xb*M1m/((2 - c)*km1), Mgm/(c*kmg))
        eD = amax(xb*M1d/((2 - c)*km1), Mgd/(c*kmg))
    else:
        km1 = ct*bl1 - t1
        if not (km1.lower() > 0 and kmg.lower() > 0): return fail('km')
        th_max = 1 - amax(D1/km1, Dg/kmg)
        dx = (2 - c)*epx
        eW = amax((dx*A(v['Bm1']) + xb*M1m)/((2 - c)*km1), Mgm/(c*kmg))
        eD = amax((dx*D1 + xb*M1d)/((2 - c)*km1), Mgd/(c*kmg))
    thK = th_max - eW - eD
    out.update(eW=eW, eD=eD, th_max=th_max, thK=thK)
    if not (thK.lower() > 0 and eW.upper() < 1): return fail('thK')
    P_1, P_g = A(v['P_1']), A(v['P_g'])
    a4 = xb*(nu1**(-th)*P_1 + Q1)/(2 - c); b4 = ((2 - c)*g**(-th)*P_g + Qg)/c
    E4 = amax(a4/km1, b4/kmg); e3 = ph1*(nu1**(-th)*P_1 + Q1)
    lamA = arb(lam)
    y_m = A(v['Bm1']) + M1m; y_l = D1 + M1d
    dR = (1 - th_max + eD)/(1 - eW)
    Y = y_m*dR + y_l; Z2 = 2*y_m + lamA*y_l + lamA*y_m*dR + dR*y_l; Z3 = 2*y_m**2 + 2*lamA*y_m*y_l + y_l**2
    q = 1 - P1
    leak = 2*E4*ph1*Y/(q*thK) + e3**2/q**2 + 2*E4*e3*ph1*Z2/(q**2*thK) + E4**2*ph1**2*Z3/(q**2*thK*(1 - eW))
    m = amin(thK, 1 - leak)
    return dict(out, margin=float(m.lower()), why='ok' if m.lower() > 0 else 'leak', leak=leak)

# ---------------- route 3 (Lemma F4) with certified J constants ----------------
def route3_extra(pr, P, b1, b2, lam, mode):
    JD = Jdata(P)
    if JD is None or (mode == 'd' and lam != 1): return (math.inf, math.inf)
    S1, S2 = pr.S1, pr.S2
    Sa = 2*iv(pr.s1, pr.s2)**(ONE/3); Ba = iv(b1, b2); s23 = iv(pr.s1, pr.s2)**(-2*ONE/3)
    Wt = Weight('m', Ba, Sa) if mode == 'm' else Weight('d1n', Ba, Sa)
    lo, hi, ds = pr.lo, pr.hi, pr.ds
    res = []
    for key in ('g1', '1g'):
        k, f = key[0], key[1]
        kb, fb = pr.B[k], pr.B[f]; kap = kb['kap']; Ck = kb['C']; jd = JD[key]
        nuA = kb['nuA']; KB = dict(nu=nuA, r=K['a']*nuA**(-ONE/3), at=K['at']*nuA**(-ONE/3), C=2*PI*nuA**(-2*ONE/3)*s23)
        r, at = KB['r'], KB['at']
        h0 = KB['C']*amin(nuA**(ONE/3)*K['Pipsi']*(Wt.val(ZERO) + Wt.wS()*(-r*Sa).exp() + Wt.G(ZERO, r)), K['Pi']*Wt.ThR(ZERO, at))
        hS = KB['C']*amin(nuA**(ONE/3)*K['Pipsi']*(Wt.val(Sa) + Wt.w0()*(-r*Sa).exp() + Wt.H(Sa, r)), K['Pi']*Wt.ThX(Sa, at))
        h0, hS = float(h0.upper()), float(hS.upper())
        Omf = pr.OmP[f]; OmS = pr.OmPS[f]
        n = pr.n; mm = np.arange(n)
        Qk = T('Q', kap*np.maximum(mm - 1, 0)*ds); Qsk = T('Qs', kap*np.maximum(mm - 1, 0)*ds)
        # Phi_Q(i) <= sum_{j>=i} Om_j Q(kap (j-i-1)_+ ds) ds ; Phi_Q*(i) <= sum_{j<=i} Om_j Q*(kap (i-j-1)_+ ds) ds
        PhiQ = np.correlate(Omf, Qk, mode='full')[n - 1:]*ds if False else np.convolve(Omf[::-1], Qk)[:n][::-1]*ds
        PhiQs = np.convolve(Omf, Qsk)[:n]*ds
        if mode == 'm':
            wS = math.exp(-b1*S1); wp = b2*np.exp(-b1*lo)
            bj = max(b for b in jd['beta'] if b <= b1); JG = jd['JG'][jd['beta'].index(bj)]
            main = JG + np.sum(Omf*(1 - np.exp(-b2*hi))*T('Q', kap*lo))*ds
        else:
            wS = 1 + math.exp(-b1*S1)
            Kn = float((2/(arb(S1)*psi(arb(b1)*arb(S1)))).upper())
            wp = Kn*np.cosh(np.minimum(b2*hi, 700))
            bo = float((ONE/(arb(S1)*psi(-arb(b2)*arb(S1)))).upper())      # beta/(1-e^{-beta S}) upper
            main = min(2*bo*jd['JsH'], wS*jd['JH'])
        tS = wS*np.sum(Omf*T('Qs', kap*np.maximum(S1 - hi, 0.0)))*ds
        tw = np.sum(wp*(PhiQ + PhiQs))*ds
        res.append((h0 + OmS*hS + Ck*(main + tS + tw))*UP)
    return tuple(res)

# ---------------- driver plug-in ----------------
import crit3 as C3F
C3F._PT = dict(C3F._PT)      # materialise the lazy NpzFile before fork (BadZipFile race between workers)
def CRIT(P, s1, s2, b1, b2, ch, coupling):
    return crit_box3(P, s1, s2, b1, b2, ch['lam'], ch['th_rho'], ch['cone'], coupling=coupling,
                     extra_routes=(route3_extra if os.environ.get('G1_ROUTE3', '1') == '1' else None))
_NOR = {}
def cands(NJG, s, b, coupling):
    """float heuristic: crit3 (routes 3/4 disabled, consistent with this port) best choice, then best with other lambda"""
    r3, r4 = C3F.route3, C3F.route4
    if os.environ.get('G1_ROUTE3', '1') != '1': C3F.route3 = lambda *a, **k: np.inf
    C3F.route4 = lambda *a, **k: np.inf
    try:
        sc, o = C3F.best(*NJG, s, max(b, 1e-9), coupling=coupling)
        yield dict(lam=o['lam'], th_rho=o['th_rho'], cone=o['cone'])
        sc2, o2 = C3F.best(*NJG, s, max(b, 1e-9), lams=(1.0 - o['lam'],), coupling=coupling)
        yield dict(lam=o2['lam'], th_rho=o2['th_rho'], cone=o2['cone'])
    finally:
        C3F.route3, C3F.route4 = r3, r4
