"""Funcoes de dois pontos e funcoes na reta de bordo (eta=0.05) para o crit3.  [NUMERICO; a certificar em Arb]
   Reta: z(t) = t e^{-i pi/6} - eta e^{i pi/3};  e := e^{-i pi/6};  z(t)+l e = z(t+l).
   Tabelas (majorantes nao-crescentes em l, grade dl=0.025 em [0,16]):
     OmA(l)   = sup_H |Ai(Z+le)/Ai(Z)|          OmPsi(l) = sup_H |psi(Z+le)/psi(Z)|
     KPi(l)   = sup_H |Ai(Z+le) Ai(wb Z)|       KPi1(l)  = sup_H max(|Ai'(Z+le)Ai(wbZ)|, |Ai(Z+le)Ai'(wbZ)|)
     Q(l)     = sup_H |Ai'(Z+le) psi(wb Z)|     Qs(l)    = sup_H |Ai'(wb Z) psi(Z+le)|
   Funcoes na reta: |H|, G(t0,l) = H(z(t0+l)) + w Ai'(z(t0+l)) psi(wb z(t0)),  H = wb Ai'(wb Z)psi(Z) - w Ai'(Z)psi(wb Z)."""
import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
import numpy as np, os, mpmath as mp
from scipy.special import airy
from scipy.integrate import cumulative_trapezoid
HERE = os.path.dirname(os.path.abspath(__file__))
eta = 0.05; EM = np.exp(-1j*np.pi/6); E3 = np.exp(1j*np.pi/3); wb = np.exp(-2j*np.pi/3); w = np.exp(2j*np.pi/3)
DT = 0.0025; T = 32.0
def line():
    t = np.arange(-T, T + DT/2, DT); Z = t*EM - eta*E3
    A, Ap, _, _ = airy(Z); B, Bp, _, _ = airy(wb*Z)
    mp.mp.dps = 30; psiR = complex(mp.quad(lambda u: mp.airyai(mp.mpc(Z[-1]) + u), [0, 2, 8, mp.inf]))
    # integracao de psi' = -Ai*EM da direita para a esquerda (Simpson composto via cumtrapz fino + correcao Richardson)
    f = (A*EM)[::-1]; tt = t[::-1]
    I1 = cumulative_trapezoid(f, tt, initial=0)
    I2 = np.zeros_like(I1); I2[::2] = cumulative_trapezoid(f[::2], tt[::2], initial=0)
    I = I1.copy(); I[::2] = I1[::2] + (I1[::2] - I2[::2])/3
    psi = (psiR - I)[::-1]
    psi[1::2] = 0.5*(psi[:-1:2] + psi[2::2]) + (psi[1::2] - 0.5*(psi[:-1:2] + psi[2::2]))  # (mantem)
    psiwb = np.conj(psi[::-1])
    H = np.conj(w)*Bp*psi - w*Ap*psiwb
    return dict(t=t, Z=Z, A=A, Ap=Ap, B=B, Bp=Bp, psi=psi, psiwb=psiwb, H=H)
def build():
    d = line(); t = d['t']; n = len(t)
    ut = np.real(-EM*d['Ap']/d['A']); uF = np.real(EM*d['A']/d['psi'])
    cu = np.concatenate([[0], np.cumsum((ut[1:] + ut[:-1])/2*DT)]); cF = np.concatenate([[0], np.cumsum((uF[1:] + uF[:-1])/2*DT)])
    L = np.arange(0, 16.0001, 0.025); out = {k: [] for k in ['OmA', 'OmPsi', 'KPi', 'KPi1', 'Q', 'Qs']}
    aA, aAp, aB, aBp, apsi, apw = [abs(d[k]) for k in ['A', 'Ap', 'B', 'Bp', 'psi', 'psiwb']]
    for l in L:
        k = int(round(l/DT)); m = n - k
        out['OmA'].append(np.exp(-(cu[k:] - cu[:m]).min())); out['OmPsi'].append(np.exp(-(cF[k:] - cF[:m]).min()))
        out['KPi'].append((aA[k:]*aB[:m]).max())
        out['KPi1'].append(max((aAp[k:]*aB[:m]).max(), (aA[k:]*aBp[:m]).max()))
        out['Q'].append((aAp[k:]*apw[:m]).max()); out['Qs'].append((aBp[:m]*apsi[k:]).max())
    for key in out:
        v = np.array(out[key])*1.002            # margem numerica 0,2%
        out[key] = np.maximum.accumulate(v[::-1])[::-1]
    np.savez(os.path.join(HERE, 'tp_tables2.npz'), L=L, ut_min=ut.min(), uF_min=uF.min(),
             HL1=np.trapezoid(abs(d['H']), t), Hsup=abs(d['H']).max(), **out)
    return out
_tab = None
def tab():
    global _tab
    if _tab is None:
        f = os.path.join(HERE, 'tp_tables2.npz')
        if not os.path.exists(f): build()
        _tab = dict(np.load(f))
        # tabelas CERTIFICADAS (agente Arb, v2) para OmA e OmPsi: mesma grade L
        import json
        cf = _os.path.join(_CERT, 'g0_constants', 'tp_cert_0.05.json')
        if os.path.exists(cf):
            dc = json.load(open(cf))
            assert np.allclose(np.array(dc['L']), _tab['L'])
            _tab['OmA'] = np.array(dc['OmA']); _tab['OmPsi'] = np.array(dc['OmPsi'])
    return _tab
def F(name, l, rate_tail=None):
    """majorante em degrau (valor no no' a esquerda) de uma tabela nao-crescente; cauda: submultiplicatividade/taxa."""
    tb = tab(); L = tb['L']; v = tb[name]; l = np.asarray(l, float)
    idx = np.clip(np.floor(l/(L[1] - L[0]) + 1e-12).astype(int), 0, len(L) - 1)
    out = v[idx]
    if rate_tail is not None:
        out = np.where(l > L[-1], v[-1]*np.exp(-rate_tail*(l - L[-1])), out)
    return out
_lin = None
def lin():
    global _lin
    if _lin is None: _lin = line()
    return _lin
def J(nu_k, nu_f, weight='1', smax=40.0, ds=0.01, t0s=None, beta=0.0):
    """sup_{t0} int_0^smax rho(s) |F(z(t0); nu_k^{-1/3} s)| ds, com rho(s) = OmPsi(nu_f^{-1/3} s) (weight '1') ou s*OmPsi(...) ('s');
       F = G (two-point, weight '1G') ou H (one-point)."""
    d = lin(); t = d['t']; aH = abs(d['H'])
    s = np.arange(0, smax, ds); rho = F('OmPsi', nu_f**(-1/3)*s, rate_tail=0.4445)
    if weight.startswith('s'): rho = rho*s
    rho = rho*np.exp(-beta*s)
    useG = weight.endswith('G')
    if t0s is None: t0s = np.arange(-24, 14, 0.02)
    best = 0.0; sc = nu_k**(-1/3)
    for t0 in t0s:
        ip = np.clip(np.round((t0 + sc*s + T)/DT).astype(int), 0, len(t) - 1)
        if useG:
            i0 = int(round((t0 + T)/DT))
            val = abs(d['H'][ip] + w*d['Ap'][ip]*d['psiwb'][i0])
        else:
            val = aH[ip]
        best = max(best, np.trapezoid(rho*val, s))
    return best*1.01            # margem 1% (discretizacao)
if __name__ == '__main__':
    out = build(); tb = tab()
    print('min ut=%.5f min uF=%.5f  ||H||_L1=%.4f sup|H|=%.4f' % (tb['ut_min'], tb['uF_min'], tb['HL1'], tb['Hsup']))
    for k in out: print('%-6s' % k, ' '.join('%.4g' % out[k][i] for i in [0, 20, 40, 80, 120, 160, 240, 320, 480, 640]), ' int=%.4f' % np.trapezoid(out[k], tb['L']))
    for (n1, g) in [(1.7, 1/3), (2, 0.5)]:
        print('nu1=%.2f g=%.3f: J_G(g<-nu1)=%.4f J_H=%.4f J_sH=%.4f | J_G(nu1<-g)=%.4f J_H=%.4f J_sH=%.4f' % (
            n1, g, J(g, n1, '1G'), J(g, n1, '1'), J(g, n1, 's'), J(n1, g, '1G'), J(n1, g, '1'), J(n1, g, 's')))
def J1(nu_k, nu_f, beta, smax=30.0, ds=0.02, t0s=None):
    """Rota 1 conjunta (Lema F5): sup_{t0} int_0^inf e^{-beta s}[ |Ai(z(t0+ks))Ai(wb z(t0))|
         + k int_0^inf OmPsi(nu_f^{-1/3} t) |dK(s,t)| dt ] ds,  k = nu_k^{-1/3},
       |dK(s,t)| = |Ai(z(t0+ks)) Ai'(wb z(t0+kt))| (t<s) ou |Ai'(z(t0+kt)) Ai(wb z(t0+ks))| (t>s)."""
    d = lin(); n = len(d['t']); aA, aAp, aB, aBp = [abs(d[q]) for q in ('A', 'Ap', 'B', 'Bp')]
    s = np.arange(0, smax, ds); k = nu_k**(-1/3)
    Om = F('OmPsi', nu_f**(-1/3)*s, rate_tail=0.4445); Om = np.minimum(Om, np.exp(-0.4445*nu_f**(-1/3)*s))
    e = np.exp(-beta*s)
    if t0s is None: t0s = np.arange(-24, 14, 0.05)
    best = 0.0
    for t0 in t0s:
        ip = np.clip(np.round((t0 + k*s + T)/DT).astype(int), 0, n - 1); i0 = int(round((t0 + T)/DT))
        A_s, X_s, Ap_t, Xp_t = aA[ip], aB[ip], aAp[ip], aBp[ip]
        term0 = np.sum(e*A_s)*ds*aB[i0]
        cumX = np.cumsum(Om*Xp_t)*ds                    # int_{t<s} Om Xp
        part1 = np.sum(e*A_s*cumX)*ds
        cumE = np.cumsum(e*X_s)*ds                      # int_{s<t} e X
        part2 = np.sum(Om*Ap_t*cumE)*ds
        best = max(best, term0 + k*(part1 + part2))
    return best*1.02
def Mjoint(nu_k, nu_f, beta=0.0, wkind='exp', smax=36.0, ds=0.01, tk=None, tf=None, verbose=False):
    """Constante conjunta de 2 variaveis (Lema F6): sup_{t_k, t_f} |nu_f^{-1/3} int_0^inf h(t) Ai(z(t_f + nu_f^{-1/3} t)) dt / psi(z(t_f))|,
       h(t) = R(t) A_X(t) + X(t) A_R(t)  (nucleo do bloco k com inicio z(t_k), S = inf),
       A_X(t) = int_0^t w X, A_R(t) = int_t^inf w R, w(s) = e^{-beta s} (wkind 'exp') ou s (wkind 's').
       Pelo PL iterado (holomorfa e limitada em cada variavel), e' o sup sobre H_delta x H_delta."""
    d = lin(); n = len(d['t']); A, B, psi = d['A'], d['B'], d['psi']
    s = np.arange(0, smax, ds); kk = nu_k**(-1/3); kf = nu_f**(-1/3)
    w = np.exp(-beta*s) if wkind == 'exp' else s
    if tk is None: tk = np.arange(-18, 10, 0.1)
    if tf is None: tf = np.arange(-24, 14, 0.1)
    jf = np.clip(np.round((tf[:, None] + kf*s[None, :] + T)/DT).astype(int), 0, n - 1)
    AiF = A[jf]                                   # (n_tf, n_s)
    psiF = psi[np.clip(np.round((tf + T)/DT).astype(int), 0, n - 1)]
    best = 0.0; arg = None
    for t0 in tk:
        ip = np.clip(np.round((t0 + kk*s + T)/DT).astype(int), 0, n - 1)
        R = A[ip]; X = B[ip]
        AX = np.concatenate([[0], np.cumsum((w*X)[1:] + (w*X)[:-1])*ds/2])
        tot = np.concatenate([[0], np.cumsum((w*R)[1:] + (w*R)[:-1])*ds/2]); AR = tot[-1] - tot
        h = R*AX + X*AR
        val = np.abs(kf*np.trapezoid(AiF*h[None, :], s, axis=1)/psiF)
        i = val.argmax()
        if val[i] > best: best = val[i]; arg = (t0, tf[i])
    if verbose: print('  argmax (t_k, t_f) =', arg)
    return best*1.02
