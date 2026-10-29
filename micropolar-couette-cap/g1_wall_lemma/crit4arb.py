"""Criterion of Theorem C in Arb with route 4 (Lemma F5, joint constants) added to crit3arb (crit3 + Lemma F8 + route 3).
For each box the moment bounds are the minimum over routes 1, 2, 3 and 4; crit3arb itself is unchanged.
Route 4 uses the certified joint constants M_u(beta_j), M_v(beta_j) (nodes beta_j in [0, 1]) and M^s_u, M^s_v
(certificates/g1_wall_lemma/route4/M_cert_1_1_0.5.json, produced by g1_wall_lemma/route4/mjoint_cert.py):
  weight m (w = e^{-beta s}, beta in [b1, b2], beta_j = largest node <= b1):
     |int_0^S w y|/|phi_f| <= C_k M(beta_j) + (b2 - beta_j) C_k int_0^inf s e^{-beta_j s} U~(s) ds + corr_S        (Lemma F5(c))
  weight d, lambda = 1 (normalised w^ = 2 sinh(beta s)/(e^{beta S} - 1), beta S < 30):
     |int w^ y|/|phi_f| <= Kn (C_k M^s + corr^(s)_S) + Kn b2^2/6 int_0^S s^3 cosh(b2 s) U(s) ds,  Kn >= 2 beta/(e^{beta S}-1)   (Lemma F5(d))
  with U_inf = C_k U~, U~(s) = K_Pi(kap_k s) + kap_k int_0^inf K_Pi1(kap_k |s-t|) Omega_f(t) dt, Omega_f(t) = Omega_psi(kap_f t),
  corr_S = C_k [int_0^S |w| (K_Pi(kap_k(S-s)) Omega_f(S) + Omega_f(S) int_{kap_k(S-s)}^inf K_Pi1) ds + int_S^inf |w| U~].
Plug-in for driver.py:  G1_CRIT=crit4arb."""
import os, sys, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'heuristic'))
import crit3arb as C3A
from crit3arb import T, UP, DN, crit_box3, route3_extra, Jdata, params, K, ONE, arb, psi
_R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
MFILE = os.environ.get('G1_MCERT', os.path.join(_R, 'certificates', 'g1_wall_lemma', 'route4', 'M_cert_1_1_0.5.json'))
_M = None
def Mdata():
    global _M
    if _M is None:
        d = json.load(open(MFILE))
        assert d['NJG'] == [1, 1, 0.5] and d['all_ok'] is True
        nodes = sorted((float(b), v) for b, v in d['M_exp'].items())
        _M = dict(nodes=[b for b, v in nodes], Mu=[v[0] for b, v in nodes], Mv=[v[1] for b, v in nodes], Ms=d['M_s'])
    return _M
AF = 0.44452699; ATF = 0.35346270

def _tailint_table(name, L):
    """G(x) >= int_x^inf T(name, l) dl on the grid x_i = i*0.01 (nonincreasing; evaluate at a lower argument)"""
    dx = 0.01; x = np.arange(0, L + dx/2, dx)
    v = T(name, x)                                   # left-node (upper) values on [x_i, x_i + dx]
    c = ATF*DN
    tail = float(T(name, np.array([x[-1]]))[0])/c   # beyond L: T <= T(L) e^{-at (l - L)}
    G = (np.cumsum(v[::-1])[::-1]*dx + tail)*UP
    return x, G
_TK1 = _tailint_table('KPi1', 120.0)
def TKP1int(x):
    x = np.maximum(np.asarray(x, float), 0.0); idx = np.floor(x/0.01*(1 - 1e-12)).astype(int)
    xs, G = _TK1
    return np.where(idx < len(G), G[np.minimum(idx, len(G) - 1)], G[-1]*np.exp(-ATF*DN*(x - xs[-1])))

_UT = {}
def Utilde(kap_lo, kap_hi, kapf_lo, Lmax=200.0, ds=0.01):
    """cells [i ds, (i+1) ds] of [0, Lmax]: upper bound of U~(s) (without C_k)"""
    key = (kap_lo, kapf_lo)
    if key not in _UT:
        n = int(round(Lmax/ds)); lo = np.arange(n)*ds
        Om = T('OmPsi', kapf_lo*lo)
        m = np.arange(n); kern = T('KPi1', kap_lo*np.maximum(m - 1, 0)*ds)
        conv = np.convolve(Om, np.r_[kern[::-1], kern[1:]])[n - 1:2*n - 1]*ds
        conv += float(K['Pi1'].upper())*float(T('OmPsi', np.array([kapf_lo*Lmax]))[0])/(AF*kapf_lo)
        U = (T('KPi', kap_lo*lo) + kap_hi*conv*UP)*UP
        # tail s > Lmax: U~(s) <= Pi e^{-at kap s} + kap (K_Pi1(kap s/2) I_Om + Pi1 Omega_psi(kapf s/2)/(a kapf)) <= Cs e^{-c s/2}
        IOm = (np.sum(Om)*ds + float(T('OmPsi', np.array([kapf_lo*Lmax]))[0])/(AF*kapf_lo))*UP
        Cs = float(K['Pi'].upper()) + kap_hi*(float(K['Pi1'].upper())*IOm + float(K['Pi1'].upper())/(AF*kapf_lo))
        c = min(ATF*kap_lo, AF*kapf_lo)*DN
        _UT[key] = dict(lo=lo, hi=lo + ds, U=U, ds=ds, L=Lmax, Cs=Cs, c=c)
    return _UT[key]
def int_tail(u, b, a0, moment):
    """upper bound of int_{a0}^inf s^moment e^{-b s} U~(s) ds"""
    lo, hi, U, ds = u['lo'], u['hi'], u['U'], u['ds']
    sel = hi > a0
    lo_c = np.maximum(lo[sel], a0)
    w = np.exp(-b*lo_c)*(hi[sel]**moment if moment else 1.0)
    acc = np.sum(w*U[sel]*(hi[sel] - lo_c))
    L, Cs, c = u['L'], u['Cs'], u['c'] + b
    L0 = max(L, a0); cc = u['c']/2 + b
    acc += Cs*math.exp(-cc*L0)*((L0/cc + 1/cc**2) if moment else 1/cc)
    return acc*UP

# extra relative-cone table (Lemma F8) with sigma_min = 228 on beta in [0, 0.45] (certificates/g1_wall_lemma/route4/); it is
# used in addition to the tables read by crit3arb (each valid table gives an upper bound; the minimum is taken)
T8X = os.environ.get('G1_T8X', os.path.join(_R, 'certificates', 'g1_wall_lemma', 'route4', 'T8_cert_1_1_0.5_s228.json'))
_T8X = [json.load(open(T8X))] if os.path.exists(T8X) else []
_theta3 = C3A.Theta_bar
def Theta_bar(P, b1, b2, s1):
    best = _theta3(P, b1, b2, s1)
    for d in _T8X:
        if [str(x) for x in d['NJG']] != ['1', '1', '0.5'] or s1 < d['sigma_min']: continue
        lo_cov = min(a for a, b, th in d['Theta_edge']); hi_cov = max(b for a, b, th in d['Theta_edge'])
        if b1 < lo_cov or b2 > hi_cov: continue
        v = None; okc = True
        for a, b, th in d['Theta_edge']:
            if a <= b2 and b >= b1:
                if not (th < 10): okc = False; break
                v = th if v is None else max(v, th)
        if okc and v is not None: best = v if best is None else min(best, v)
    return best
C3A.Theta_bar = Theta_bar          # crit_box3 looks the function up in its module at call time

def route4_extra(pr, P, b1, b2, lam, mode):
    MD = Mdata()
    if mode == 'd' and lam != 1: return (math.inf, math.inf)
    S1, S2 = pr.S1, pr.S2
    res = []
    for key, idx in (('g1', 0), ('1g', 1)):
        k, f = key[0], key[1]
        kb, fb = pr.B[k], pr.B[f]
        Ck = kb['C']; kap_lo, kap_hi = kb['kap'], kb['kapU']; kapf_lo = fb['kap']
        u = Utilde(kap_lo, kap_hi, kapf_lo)
        OmS = pr.OmPS[f]
        lo, hi, ds = pr.lo, pr.hi, pr.ds
        sep = kap_lo*np.maximum(S1 - hi, 0.0)
        inner = T('KPi', sep)*OmS + OmS*TKP1int(sep)
        if mode == 'm':
            nodes = MD['nodes']
            if b1 > nodes[-1]: res.append(math.inf); continue
            j = max(i for i, b in enumerate(nodes) if b <= b1); bj = nodes[j]
            M = (MD['Mu'] if idx == 0 else MD['Mv'])[j]
            lip = (b2 - bj)*Ck*int_tail(u, bj, 0.0, 1)
            corr = Ck*np.sum(np.exp(-b1*lo)*inner)*ds + Ck*int_tail(u, b1, S1, 0)
            res.append((Ck*M + lip + corr)*UP)
        else:
            if not b2*S2 < 30: res.append(math.inf); continue
            Kn = float((2/(arb(S1)*psi(arb(b1)*arb(S1)))).upper())
            Ms = MD['Ms'][idx]
            corrs = Ck*np.sum(hi*inner)*ds + Ck*int_tail(u, 0.0, S1, 1)
            Up = pr.pair[key]['U']
            rem = Kn*b2**2/6*np.sum(hi**3*np.cosh(np.minimum(b2*hi, 700))*Up)*ds
            res.append((Kn*(Ck*Ms + corrs) + rem)*UP)
    return tuple(res)

def routes34(pr, P, b1, b2, lam, mode):
    r3 = route3_extra(pr, P, b1, b2, lam, mode)
    r4 = route4_extra(pr, P, b1, b2, lam, mode)
    return (min(r3[0], r4[0]), min(r3[1], r4[1]))

import crit3 as C3F
C3F._PT = dict(C3F._PT)
def CRIT(P, s1, s2, b1, b2, ch, coupling):
    return crit_box3(P, s1, s2, b1, b2, ch['lam'], ch['th_rho'], ch['cone'], coupling=coupling, extra_routes=routes34)
def cands(NJG, s, b, coupling):
    """float heuristic: crit3 with routes 3 and 4 (numerical joint constants of heuristic/Mtab_*.json): best choice, then
    best with the other lambda.  Only the choice of the free parameters depends on it."""
    sc, o = C3F.best(*NJG, s, max(b, 1e-9), coupling=coupling)
    yield dict(lam=o['lam'], th_rho=o['th_rho'], cone=o['cone'])
    sc2, o2 = C3F.best(*NJG, s, max(b, 1e-9), lams=(1.0 - o['lam'],), coupling=coupling)
    yield dict(lam=o2['lam'], th_rho=o2['th_rho'], cone=o2['cone'])
