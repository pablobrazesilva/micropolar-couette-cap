"""Certification of the crit3 route-3 constants J_G(beta), J_H, J_sH (Lemma F4, eta=0.05), one triple.
   J_G(beta) = sup_{t0} int_0^inf Om_f(s) e^{-beta s} |G(z(t0); kap_k s)| ds,  G(Z0;l) = H(Z0+le) + w Ai'(Z0+le) psi(wb Z0)
   J_H = sup_{t0} int Om_f(s)|H(z(t0)+kap_k s e)| ds,  J_sH: weight s Om_f(s);  H = wb Ai'(wb Z)psi(Z) - w Ai'(Z)psi(wb Z).
   Sup over H_eta = sup over the line (Lemma F2).  On the line: z(t)+le = z(t+l), psi(wb z(t)) = conj psi(z(-t)).
   METHOD: Arb point values (Ai, Ai', Ai(wb z), Ai'(wb z), psi = Ai/F with F from the G0 integral representation) on the
   half-grid t_k = k h/2, |t| <= TG.  For an (a-cell, s-cell) pair with p = a + kap s, a Taylor model at the midpoint
   (a_m, p_m) (half-grid points):  G = G_m + Phi_a da + Phi_p dp + R,  |da|<=h/2, |dp|<=h, |R| <= (Maa h^2/4 + 2 Map h^2/2 + Mpp h^2)/2,
     Phi_a = -e Ai'(z(p)) Ai(wb z(a)),  Phi_p = e[z H0(z) - e^{i pi/6}/(2 pi)](p) + w e z(p) Ai(z(p)) psi(wb z(a)),
     Maa <= rhoA(p) K_Pi1(l),  Map <= |z(p)| K_Pi(l),  Mpp <= 2 Pi_psi + |z|(Q(0)+Q*(0)) + Om_A(l) Pi_psi + |z| Q(l)   (l = kap s_lo)
   (H' = Z H0 - e^{i pi/6}/(2pi), H0 = Ai(wbZ)psi(Z) - w Ai(Z)psi(wbZ), |H0| <= 2 Pi_psi, H'' = H0 + Z H).  Sup of the linear model
   over the box at its corners.  rhoA = |Ai'/Ai| on [p_m-h, p_m+h] by the Riccati a-priori bound.  Tails: explicit (see code).
   usage: cert_J2.py N J G"""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import sys, json, math
sys.path.insert(0, _os.path.join(_R, 'g0_constants'))
from cl import *
from flint import acb, arb, ctx
import numpy as np
from multiprocessing import Pool
ctx.prec = 64
eta = '0.05'; ETA = 0.05
h = 0.025; hh = h/2; TG = 52.0
EMc = acb(arb(3).sqrt()/2, -arb(1)/2); E3c = acb(arb(1)/2, arb(3).sqrt()/2); WB = acb(-arb(1)/2, -arb(3).sqrt()/2)
af = 0.44452699
CD = _os.path.join(_CERT, 'g0_constants') + _os.sep; CDW = _os.path.dirname(_wout('g0_constants', 'x')) + _os.sep
T2 = np.array(json.load(open(CD + 'tp_cert_0.05.json'))['OmPsi']); T1 = np.array(json.load(open(CD + 'tp_cert_0.05.json'))['OmA'])
T36 = json.load(open(CD + 'T36_cert_0.05.json'))
cp = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pi_0.05.json'))); cq = json.load(open(_os.path.join(_CERT, 'g1_wall_lemma', 'cert_pipsi.json')))
Pipsi = cq['Pi_psi']; Q0 = T36['Q'][0]; Qs0 = T36['Qs'][0]
def tab(v, l, rate):
    l = np.maximum(np.asarray(l, float), 0); idx = np.floor(l*(1 - 1e-12)/0.025).astype(int)
    return np.where(idx < len(v), v[np.minimum(idx, len(v) - 1)], v[-1]*np.exp(-rate*(l - 16)*(1 - 1e-9)))*(1 + 1e-9)
OmPsi = lambda l: np.minimum(tab(T2, l, af), np.exp(-af*np.asarray(l)*(1 - 1e-9))*(1 + 1e-9))
OmA = lambda l: tab(T1, l, 0.35346270)
KPi = lambda l: tab(np.array(T36['KPi']), l, 0.35346270)
KPi1 = lambda l: tab(np.array(T36['KPi1']), l, 0.35346270)
Qt = lambda l: tab(np.array(T36['Q']), l, 0.35346270)
def cr(x):
    c = complex(float(x.real.mid()), float(x.imag.mid())); r = ub(abs(x - acb(c.real, c.imag)))
    return c.real, c.imag, r
def pdata(k):
    t = -TG + k*hh
    z = acb(arb(repr(t)))*EMc - arb(eta)*E3c
    F, u, ut = point(repr(t), eta, 64)
    ai, aip, _, _ = z.airy(); bi, bip, _, _ = (WB*z).airy()
    return (k,) + cr(ai) + cr(aip) + cr(bi) + cr(bip) + cr(ai/F) + cr(z)
if __name__ == '__main__':
    N, J, G = sys.argv[1:4]
    Nf, Jf, Gf = float(N), float(J), float(G); g = Gf/Jf; nu1 = 1 + Nf
    nK = int(round(2*TG/hh)) + 1
    import os, pickle
    cache = CDW + 'Jpoints_h%g_T%g.pkl' % (h, TG)
    if os.path.exists(cache):
        D = pickle.load(open(cache, 'rb'))
    else:
        with Pool(2) as P: res = P.map(pdata, range(nK), chunksize=64)
        res.sort(); D = np.array([r[1:] for r in res]); pickle.dump(D, open(cache, 'wb'))
    def C(col): return D[:, col] + 1j*D[:, col + 1], D[:, col + 2]
    Ai, rAi = C(0); Aip, rAip = C(3); Bi, rBi = C(6); Bip, rBip = C(9); Psi, rPsi = C(12); Z, rZ = C(15)
    tk = -TG + np.arange(nK)*hh
    # psi(wb z(t)) = conj psi(z(-t)) : index nK-1-k
    Pw, rPw = np.conj(Psi[::-1]), rPsi[::-1]
    w = complex(-0.5, math.sqrt(3)/2); wb = np.conj(w); e = complex(math.sqrt(3)/2, -0.5); eip6 = complex(math.sqrt(3)/2, 0.5)
    rel = 1e-13
    def mul(x, rx, y, ry): return x*y, np.abs(x)*ry + rx*np.abs(y) + rx*ry + rel*np.abs(x*y)
    def add(x, rx, y, ry): return x + y, rx + ry + rel*(np.abs(x) + np.abs(y))
    t1, r1 = mul(Bip, rBip, Psi, rPsi); t2, r2 = mul(Aip, rAip, Pw, rPw)
    H, rH = add(wb*t1, r1, -w*t2, r2)
    u1, s1 = mul(Bi, rBi, Psi, rPsi); u2, s2 = mul(Ai, rAi, Pw, rPw)
    H0, rH0 = add(u1, s1, -w*u2, s2)
    zH0, rzH0 = mul(Z, rZ, H0, rH0)
    Hp, rHp = e*(zH0 - eip6/(2*math.pi)), rzH0 + rel                     # dH/dt along the line
    absZ = np.abs(Z) + rZ
    rho0 = np.abs(Aip)/np.maximum(np.abs(Ai) - rAi, 1e-300) + (rAip/np.maximum(np.abs(Ai) - rAi, 1e-300))*1.0
    def rho_bound(k, span):   # |Ai'/Ai| on [t_k - span, t_k + span] (Riccati a priori bound)
        Zm = absZ[k] + span; r0 = rho0[k]; disc = 1 - 4*span*(r0 + span*Zm)
        return np.where(disc > 0, (1 - np.sqrt(np.maximum(disc, 0)))/(2*span), np.inf)
    kap = {'1': nu1**(-1/3), 'g': g**(-1/3)}
    BG = [0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3, 5, 8, 12, 20, 40]
    out = {}
    def tails(T, side):
        th = tail_side(T, eta, side); tr = tail_side(T, eta, 'left' if side == 'right' else 'right')
        eB = max(th[3], tr[3]); eD = max(th[4], tr[4]); x = math.sqrt(T)
        P0 = 1/(2*math.pi*(2*x - 2*eD))
        return 2*(1 + eD/x)*(1 + eB)*P0*(1 + 1e-9), P0, eB, eD, th[1]
    for k_, f_ in (('g', '1'), ('1', 'g')):
        kk, kf = kap[k_]*(1 + 1e-12), kap[f_]*(1 - 1e-12)
        smax = 8.0/kf; ds = h/kk; m = int(math.ceil(smax/ds)); smax = m*ds
        slo = np.arange(m)*ds; shi = slo + ds; lo_l = kk*slo*(1 - 1e-12)
        Om = OmPsi(kf*slo); Wexp = np.array([np.exp(-b*slo) for b in BG])
        KP, KP1, OA, QQ = KPi(lo_l), KPi1(lo_l), OmA(lo_l), Qt(lo_l)
        # a-cells: midpoints a_m = odd half-grid index ia; p_m = a_m + kap(j+1/2)ds = index ia + 1 + 2j
        best = np.zeros(len(BG)); bH = 0.0; bsH = 0.0
        ia_lo = int(round((-28.0 + TG)/hh)) | 1
        ia_hi = nK - 1 - 2*m - 4
        for ia in range(ia_lo, ia_hi, 2):
            ip = ia + 1 + 2*np.arange(m)
            # G_m and derivatives at (a_m, p_m)
            gx, rg = mul(Aip[ip], rAip[ip], np.full(m, Pw[ia]), np.full(m, rPw[ia]))
            Gm, rGm = add(H[ip], rH[ip], w*gx, rg)
            fa, rfa = mul(Aip[ip], rAip[ip], np.full(m, Bi[ia]), np.full(m, rBi[ia])); Phia = -e*fa
            zA, rzA = mul(Z[ip], rZ[ip], Ai[ip], rAi[ip]); q, rq = mul(zA, rzA, np.full(m, Pw[ia]), np.full(m, rPw[ia]))
            Phip, rPhip = add(Hp[ip], rHp[ip], w*e*q, rq)
            # corners: da in {+-h/2}, ds in {+-ds/2} -> dp = da + kap ds in {+-h/2 +- h/2}
            lin = np.zeros(m)
            for da in (-hh, hh):
                for dp_ in (-hh, hh):
                    lin = np.maximum(lin, np.abs(Gm + Phia*da + Phip*(da + dp_)))
            lin = lin + rGm + rfa*hh + rPhip*h
            Zmax = absZ[ip] + h
            Maa = rho_bound(ip, h)*KP1; Map = Zmax*KP
            Mpp = 2*Pipsi + Zmax*(Q0 + Qs0) + OA*Pipsi + Zmax*QQ
            R = 0.5*(Maa*hh**2 + 2*Map*hh*h + Mpp*h**2)
            Gsup = lin + R
            best = np.maximum(best, (Wexp*(Om*Gsup)).sum(axis=1)*ds)
            # H only (1D in p): |H(p)| <= max |H_m +- H'_m h| + |H''| h^2/2, |H''| <= 2Pipsi + |z|(Q0+Qs0)
            Hsup = np.maximum(np.abs(H[ip] + Hp[ip]*h), np.abs(H[ip] - Hp[ip]*h)) + rH[ip] + rHp[ip]*h + 0.5*h**2*(2*Pipsi + Zmax*(Q0 + Qs0))
            bH = max(bH, np.sum(Om*Hsup)*ds); bsH = max(bsH, np.sum(shi*Om*Hsup)*ds)
        a_lo = tk[ia_lo] - hh; a_hi = tk[ia_hi] + hh            # certified t0 range [a_lo, a_hi] (cells around midpoints)
        # ---- tails in t0 ----
        # right: t0 >= a_hi, all p >= a_hi
        HcR, P0R, eBR, eDR, rateR = tails(a_hi, 'right')
        secR = (1 + np.sqrt(kk*shi/a_hi))*(1 + eDR/math.sqrt(a_hi))*(1 + eBR)*P0R*np.exp(-kk*slo*rateR*(1 - 1e-9))
        TR = ((Wexp*(Om*(HcR + secR))).sum(axis=1)*ds, np.sum(Om*HcR)*ds, np.sum(shi*Om*HcR)*ds)
        # left: t0 <= a_lo; p <= -T* (T* = 15) where H crude-left; p > -T* only for kap s >= |a_lo| - T*
        Ts = 15.0
        HcL, P0L, eBL, eDL, rateL = tails(Ts, 'left')
        _, P0L2, eBL2, eDL2, _ = tails(abs(a_lo), 'left')
        secL = (1 + eDL2/math.sqrt(abs(a_lo)))*(1 + eBL2)*P0L2*np.exp(-kk*slo*rateL*(1 - 1e-9))
        far = (kk*shi >= abs(a_lo) - Ts)                          # p may be > -T*: crude uniform |G| <= 2Q0+Qs0
        GL = np.where(far, 2*Q0 + Qs0, HcL + secL); HL = np.where(far, Q0 + Qs0, HcL)
        TLt = ((Wexp*(Om*GL)).sum(axis=1)*ds, np.sum(Om*HL)*ds, np.sum(shi*Om*HL)*ds)
        # s-tail beyond smax: |G| <= 2Q0+Qs0, |H| <= Q0+Qs0
        Oms = float(OmPsi(np.array([kf*smax]))[0])
        stG = np.array([(2*Q0 + Qs0)*Oms*math.exp(-b*smax)/(af*kf + b) for b in BG])
        stH = (Q0 + Qs0)*Oms/(af*kf); stsH = (Q0 + Qs0)*Oms*(smax/(af*kf) + 1/(af*kf)**2)
        UP = 1 + 1e-9
        JG = np.maximum(np.maximum(best, TR[0]), TLt[0])*UP + stG
        JH = max(bH, TR[1], TLt[1])*UP + stH; JsH = max(bsH, TR[2], TLt[2])*UP + stsH
        out[k_ + f_] = dict(beta=BG, JG=[float(x) for x in JG], JH=float(JH), JsH=float(JsH), t0_range=[float(a_lo), float(a_hi)],
                            compact=dict(JG=[float(x) for x in best], JH=float(bH), JsH=float(bsH)),
                            tails=dict(right=[float(TR[0][0]), float(TR[1]), float(TR[2])], left=[float(TLt[0][0]), float(TLt[1]), float(TLt[2])]))
        print((N, J, G), k_ + f_, 'JG(0,1,5,40)=', ['%.4f' % JG[BG.index(b)] for b in (0, 1, 5, 40)], 'JH=%.4f JsH=%.4f' % (JH, JsH),
              '| compact JG0=%.4f JH=%.4f JsH=%.4f | tailsR %.4f/%.4f/%.4f tailsL %.4f/%.4f/%.4f' % (best[0], bH, bsH, TR[0][0], TR[1], TR[2], TLt[0][0], TLt[1], TLt[2]), flush=True)
    json.dump(out, open(CDW + 'J_cert_%s_%s_%s.json' % (N, J, G), 'w'), indent=1)
