"""Certification of the joint constants of route 4 (Lemma F5), triple (N,J,G) = (1,1,1/2), eta = delta = 0.05.

  f(t_k, t_f) = kf int_0^inf h(s; t_k) Ai(z(t_f + kf s)) ds / psi(z(t_f)),   kf = nu_f^{-1/3},
  h(s; t_k)   = R(s) A_X(s) + X(s) A_R(s),  R(s) = Ai(z(t_k + k s)),  X(s) = Ai(wb z(t_k + k s)),  k = nu_k^{-1/3},
  A_X(s) = int_0^s w X,  A_R(s) = int_s^inf w R,  w(s) = e^{-beta s} ('exp') or s ('s');
  M = sup_{(t_k,t_f) in R^2} |f|  (Lemma F5(b)); pairs 'u' (k = g, f = nu1) and 'v' (k = nu1, f = g).

The plane is covered by four regions, each with a rigorous upper bound (floating point with explicit error terms for
quadrature, interpolation and rounding; Airy values and F = Ai/psi from Arb):
  (1) strip t_k in [K1, K2], t_f in [-20, 30]: values, gradients and quadrature errors at the grid points
      (t_k, t_f) in (Z/128)^2 plus a second-order Taylor bound on each grid cell of half-width 1/256 (sups of the
      second derivatives per block of width 1 in t_f);
  (2) strip t_k, t_f <= -20 or t_f >= 30: integration by parts against psi, |f| <= |h(0)| + int |h'| |psi(z(t_f+kf s))/psi(z(t_f))|,
      with the certified lower bounds of u_F = Re(e F) (G0 boxes; Lemma B tails);
  (3) |t_k| <= TB outside the strip: the same IBP bound uniform in t_f, |f| <= |h(0)| + int |h'| Omega_psi(kf s) ds;
  (4) |t_k| >= TB: quasi-static representation h = h0 + (2pi/k) int G h0'' + (2pi/k) R W0(h0), h0 = k w/(2 pi Q), with the
      certified two-point bounds K_Pi, K_Pi1; the bound is nonincreasing in |t_k| and evaluated at |t_k| = TB.
usage: mjoint_cert.py pair(u|v) weight(exp|s) beta thr line_tab.npz Ftab128.npz   (prints RESULT {json})"""
import os, sys, json, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R0 = os.environ.get('CAP_REPO', os.path.abspath(os.path.join(HERE, '..', '..')))
CERT = os.path.join(R0, 'certificates')
WORK = os.environ.get('CAP_WORK', os.path.join(R0, 'work'))
sys.path.insert(0, os.path.join(R0, 'g0_constants')); sys.path.insert(0, os.path.join(R0, 'g0_sector_lemma'))

DT = 1.0/128; STEP = 1            # tau-grid, (t_k, t_f) grid and sigma grid: 1/128
UP = 1 + 1e-12                    # relative safety factor for floating-point sums of nonnegative terms (>> n u)
ETA = 0.05
E = np.exp(-1j*np.pi/6); WB = np.exp(-2j*np.pi/3); E3 = np.exp(1j*np.pi/3)

# ---------------- certified tables (step majorants on l = 0:0.025:16, exponential tails), as in crit3arb ----------------
_t12 = json.load(open(os.path.join(CERT, 'g0_constants', 'tp_cert_0.05.json')))
_t36 = json.load(open(os.path.join(CERT, 'g0_constants', 'T36_cert_0.05.json')))
_cp = json.load(open(os.path.join(CERT, 'g1_wall_lemma', 'cert_pi_0.05.json')))
_cq = json.load(open(os.path.join(CERT, 'g1_wall_lemma', 'cert_pipsi.json')))
AF, ATF = 0.44452699, 0.35346270          # certified lower bounds of a(0.05), a~(0.05)
PI_, PI1 = _cp['Pi'], _cp['Pi1']
DL = 0.025
TAB = {'OmA': np.array(_t12['OmA']), 'OmPsi': np.array(_t12['OmPsi']), 'KPi': np.array(_t36['KPi']), 'KPi1': np.array(_t36['KPi1'])}
RATE = {'OmA': ATF, 'OmPsi': AF, 'KPi': ATF, 'KPi1': ATF}
UNIF = {'OmA': 1.0, 'OmPsi': 1.0, 'KPi': PI_, 'KPi1': PI1}
DN = 1 - 1e-9
def T(name, l):
    """upper bound of the nonincreasing certified function 'name' at every argument >= l"""
    l = np.maximum(np.asarray(l, float), 0.0); v = TAB[name]; n = len(v)
    idx = np.floor(l*(1 - 1e-12)/DL).astype(np.int64)
    out = np.where(idx < n, v[np.minimum(idx, n - 1)], v[-1]*np.exp(-RATE[name]*(l - (n - 1)*DL)*DN))
    return np.minimum(out, UNIF[name]*np.exp(-RATE[name]*l*DN))*(1 + 1e-9)
def Tint(name, x, moment=0):
    """upper bound of int_x^inf l^moment T(l) dl (x scalar >= 0)"""
    v = TAB[name]; n = len(v); L = (n - 1)*DL
    acc = 0.0
    if x < L:
        i0 = int(math.floor(x/DL)); ls = np.arange(i0, n - 1)*DL
        lo = np.maximum(ls, x); hi = ls + DL
        vals = T(name, lo)
        acc = np.sum(vals*(hi**(moment + 1) - lo**(moment + 1))/(moment + 1))
    c = RATE[name]*DN; X = max(x, L); v0 = float(T(name, X))
    acc += v0*(1/c if moment == 0 else (X/c + 1/c**2))
    return acc*(1 + 1e-9)
IKPI, IKPI1 = Tint('KPi', 0), Tint('KPi1', 0)
IKPIm, IKPI1m = Tint('KPi', 0, 1), Tint('KPi1', 0, 1)

# ---------------- line table (Arb) and F table (Arb) ----------------
def load_tables(linef, ftabf):
    d = np.load(linef); f = np.load(ftabf)
    L = {k: d[k] for k in d.files}; L['imax'] = int(L['imax'])
    a = L['a']; L['absa'] = abs(a)
    Amin = np.minimum(abs(a[:-1]) - L['ra'][:-1], abs(a[1:]) - L['ra'][1:]) - DT/2*L['A1m'][:-1]
    L['Amin'] = np.append(np.maximum(Amin, 0.0), 0.0)
    L['rho'] = np.where(L['Amin'] > 0, L['A1m']/np.maximum(L['Amin'], 1e-300), np.inf)      # sup |Ai'/Ai| per cell
    L['rel_a'] = float(np.max(L['ra']/abs(a)))
    Ft = dict(k=f['k'], F=f['F'], rF=f['rF'])
    return L, Ft

# ---------------- G0 boxes: running-min lower bound of u_F = Re(e F) on [-20, 30] ----------------
def uF_profile():
    """v(tau) on tau-cells of width 1/256 from -20: nonincreasing lower bound of u_F, floored by the tail bound 1.6423
    (u_F >= 1.6423 for tau <= -20, tail_side(20, left)); used for t_f <= -20 (shift argument).  Returns (edges, v)."""
    d = json.load(open(os.path.join(CERT, 'g0_constants', 'line_0.05.json')))
    b = sorted(d['boxes']); lo = np.array([x[0] for x in b]); ulo = np.array([x[4] for x in b])
    v = np.minimum.accumulate(np.minimum(ulo, 1.6423))
    return lo, v, float(d['Fmax_up'])
U_RIGHT30 = 3.8313    # u_F >= 3.8313 for tau >= 30 (tail_side(30, right))

class Pair:
    def __init__(s, name):
        nu1, g = 2.0, 0.5
        nk, nf = (g, nu1) if name == 'u' else (nu1, g)
        s.name = name; s.kap = nk**(-1/3); s.kapf = nf**(-1/3)
        s.kap_lo, s.kapf_lo = s.kap*(1 - 1e-15), s.kapf*(1 - 1e-15)
        s.ds = DT/s.kap

# ---------------- profile along s for one t_k (tau-index i0) ----------------
def profile(L, P, i0, wk, beta, J=None):
    im = L['imax']; a, a1, ra, ra1 = L['a'], L['a1'], L['ra'], L['ra1']
    Am, A1m, Zm = L['Am'], L['A1m'], L['Zm']
    if J is None: J = im - 1 - i0
    kap, ds = P.kap, P.ds
    j = np.arange(J + 1); s = j*ds
    iR = im + i0 + j; iX = im - i0 - j
    R = a[iR]; X = np.conj(a[iX]); rR = ra[iR]; rX = ra[iX]
    Rp = kap*E*a1[iR]; Xp = kap*WB*E*np.conj(a1[iX]); rRp = kap*ra1[iR]; rXp = kap*ra1[iX]
    m = np.arange(J); cR = im + i0 + m; cX = im - i0 - m - 1
    RA, RA1 = Am[cR], A1m[cR]; XA, XA1 = Am[cX], A1m[cX]; Zc = np.maximum(Zm[cR], Zm[cX])
    if wk == 'exp':
        w = np.exp(-beta*s); wmax = w[:-1]; w1 = beta; w2 = beta*beta
        eX = ds**3/12*wmax*(w2*XA + 2*w1*kap*XA1 + kap**2*Zc*XA)
        eR = ds**3/12*wmax*(w2*RA + 2*w1*kap*RA1 + kap**2*Zc*RA)
    else:
        w = s.copy(); wmax = s[1:]
        eX = ds**3/12*(2*kap*XA1 + wmax*kap**2*Zc*XA)
        eR = ds**3/12*(2*kap*RA1 + wmax*kap**2*Zc*RA)
    tX = (w[:-1]*X[:-1] + w[1:]*X[1:])*ds/2; tR = (w[:-1]*R[:-1] + w[1:]*R[1:])*ds/2
    radX = (w[:-1]*rX[:-1] + w[1:]*rX[1:])*ds/2; radR = (w[:-1]*rR[:-1] + w[1:]*rR[1:])*ds/2
    AX = np.concatenate([[0], np.cumsum(tX)])
    EX = np.concatenate([[0], np.cumsum(eX + radX + 1e-13*abs(tX))])
    EX = EX + 1e-13*np.concatenate([[0], np.cumsum(abs(tX))])
    rev = np.cumsum(tR[::-1])[::-1]; AR = np.append(rev, 0.0)
    # tail of A_R beyond s_J: |Ai(z(tau))| <= |Ai(z(tau_J))| e^{-ut (tau - tau_J)}, ut = utlo(tail_side(tau_J, right)) >= 7
    utJ = UTJ[0] if UTJ else ut_far(J + i0)
    RJ = abs(R[-1]) + rR[-1]
    tail = RJ*(w[-1]/(utJ*kap) if wk == 'exp' else (s[-1]/(utJ*kap) + 1/(utJ*kap)**2))
    ER = np.append(np.cumsum((eR + radR + 2e-13*abs(tR))[::-1])[::-1], 0.0) + tail
    absR, absX = abs(R) + rR, abs(X) + rX
    h = R*AX + X*AR
    Eh = absR*EX + rR*abs(AX) + absX*ER + rX*abs(AR) + 1e-13*(abs(R*AX) + abs(X*AR))
    hp = Rp*AX + Xp*AR
    Ehp = (abs(Rp) + rRp)*EX + rRp*abs(AX) + (abs(Xp) + rXp)*ER + rXp*abs(AR) + 1e-13*(abs(Rp*AX) + abs(Xp*AR))
    # cell suprema of |h|, |h'|, |h''| (h'' = Q h - kap w/(2 pi), |Q| = kap^2 |z|)
    m0 = np.maximum(abs(h[:-1]) + Eh[:-1], abs(h[1:]) + Eh[1:])
    m1 = np.maximum(abs(hp[:-1]) + Ehp[:-1], abs(hp[1:]) + Ehp[1:])
    c = kap*wmax/(2*np.pi); dd = ds/2
    den = 1 - dd*dd*kap**2*Zc
    assert np.all(den > 0.5)
    H0 = (m0 + dd*m1 + dd*dd*c)/den*UP
    H2 = (kap**2*Zc*H0 + c)*UP
    H1 = (m1 + dd*H2)*UP
    return dict(s=s, w=w, R=R, X=X, Rp=Rp, Xp=Xp, h=h, hp=hp, Eh=Eh, Ehp=Ehp, H0=H0, H1=H1, H2=H2, Zc=Zc, RA=RA, RA1=RA1,
                XA=XA, XA1=XA1, wmax=wmax, J=J, X0=X[0], absX0=absX[0], absX0p=abs(Xp[0]) + rXp[0])

UTJ = []
def ut_far(iJ):
    from cl import tail_side
    v = tail_side(repr(float(int(iJ))/128), '0.05', 'right')[1]
    UTJ.append(v); return v

def unif(P, wk, beta, s):
    """uniform bounds (all t_k) of |h(s)| and |h'(s)| from |G| <= K_Pi, |d_s G| <= kap K_Pi1"""
    kap = P.kap
    if wk == 'exp':
        return 2*IKPI/kap*np.ones_like(s), 2*IKPI1*np.ones_like(s)
    return (2*s*IKPI + 2*IKPIm/kap)/kap, 2*s*IKPI1 + 2*IKPI1m/kap

def dk_bounds(P, wk, beta, pr, d):
    """per s-cell: DK0 >= sup |d_{t_k} h|, DK1 >= sup |d_{t_k} h'| over the t_k-cell of half-width d around the grid point,
    and the cell-sups of |h|, |h'| over (t_k-cell x s-cell).  Uses |X0 R(s)| <= K_Pi(kap s), |X0 R'| , |X0' R| <= kap K_Pi1."""
    kap = P.kap; s = pr['s'][:-1]
    H0u, H1u = unif(P, wk, beta, pr['s'][1:])
    KP, KP1 = T('KPi', kap*s), T('KPi1', kap*s)
    if wk == 'exp':
        DK0 = (H1u + beta*H0u + KP)/kap
        H0c = pr['H0'] + d*DK0
        H2c = kap**2*(pr['Zc'] + d)*H0c + kap*pr['wmax']/(2*np.pi)
        DK1 = (H2c + beta*H1u + kap*KP1)/kap
    else:     # d_t h = (h' - h1)/kap, h1 = profile with w = 1 (|h1| <= 2 IKPI/kap, |h1'| <= 2 IKPI1)
        DK0 = (H1u + 2*IKPI/kap)/kap
        H0c = pr['H0'] + d*DK0
        H2c = kap**2*(pr['Zc'] + d)*H0c + kap*pr['wmax']/(2*np.pi)
        DK1 = (H2c + 2*IKPI1)/kap
    H1c = pr['H1'] + d*DK1
    return DK0*UP, DK1*UP, H0c*UP, H1c*UP, H2c*UP

def ibp_bound(P, wk, beta, pr, d, kern):
    """sup over the t_k-cell of |h(0)| + int_0^inf |h'(s)| kern(kf s) ds  (kern nonincreasing, given on s-cell left nodes)"""
    DK0, DK1, H0c, H1c, H2c = dk_bounds(P, wk, beta, pr, d)
    s = pr['s']; ds = P.ds
    h0 = abs(pr['h'][0]) + pr['Eh'][0] + d*DK0[0]
    body = np.sum(H1c*kern(P.kapf_lo*s[:-1]))*ds
    # beyond s_J: |h'| <= uniform bound, kern <= Omega_psi <= e^{-a l}
    sJ = s[-1]; c = AF*P.kapf_lo
    H0u, H1u = unif(P, wk, beta, np.array([sJ]))
    tailv = (2*IKPI1*math.exp(-c*sJ)/c) if wk == 'exp' else ((2*IKPI1*(sJ/c + 1/c**2) + 2*IKPI1m/P.kap/c)*math.exp(-c*sJ))
    return (h0 + body + tailv)*UP

# ---------------- region (1): strip, second-order cells ----------------
def strip(L, Ft, P, wk, beta, kidx, sigmax=60.0, blk=1.0, verbose=False):
    """kidx: tau-indices of the strip rows (grid 1/128).  t_f grid = F-table grid (1/128) restricted to [-20, 30].
    Each grid point is the centre of a cell of half-width d = 1/256 in both variables."""
    im = L['imax']; ds = P.ds; kap, kapf = P.kap, P.kapf
    d = STEP*DT/2
    dsig = STEP*DT; Ns = int(round(sigmax/dsig)); sig = np.arange(Ns + 1)*dsig
    tw = np.full(Ns + 1, dsig); tw[0] = tw[-1] = dsig/2
    fk = Ft['k']; scale = int(round(1/(STEP*DT)))           # F table at t = k/scale
    sel = (fk >= -20*scale) & (fk <= 30*scale)
    tf_idx = (fk[sel]*(128//scale)).astype(int); Fv = Ft['F'][sel]; rF = Ft['rF'][sel]
    af = L['a'][im + tf_idx]; psi = af/Fv
    relpsi = L['rel_a'] + rF/abs(Fv) + 1e-15
    rows = tf_idx[None, :] + STEP*np.arange(Ns + 1)[:, None]
    OmA = T('OmA', sig)
    rhoc = L['rho']; Zm = L['Zm']
    from numpy.lib.stride_tricks import sliding_window_view as swv
    def wmax(arr, lo_tau, width):
        """sup of arr over tau-cells meeting [lo_tau + sigma_i, lo_tau + width + sigma_{i+1}] for every i"""
        n0 = im + int(math.floor(lo_tau*128)) - 1; W = int(math.ceil((width + dsig)*128)) + 2
        win = swv(arr[n0:n0 + STEP*Ns + W], W)[::STEP][:Ns + 1]
        return win.max(1)
    # global windows (all t_f in [-20,30]) for the error terms; local windows per t_f block for the second derivatives
    rho_w = wmax(rhoc, -20.0, 50.0); Z_w = wmax(Zm, -20.0, 50.0)
    lo_box, vprof, Fs_glob = uF_profile()
    bx = json.load(open(os.path.join(CERT, 'g0_constants', 'line_0.05.json')))['boxes']
    bxl = np.array([x[0] for x in bx]); bxF = np.array([x[3] for x in bx])
    blocks = np.arange(-20.0, 30.0, blk)
    Fs_b = np.array([bxF[(bxl >= b0 - d - 1/256) & (bxl < b0 + blk + d)].max() for b0 in blocks])
    rho_b = np.array([wmax(rhoc, b0 - d, blk + 2*d) for b0 in blocks])          # (nb, Ns+1)
    Z_b = np.array([wmax(Zm, b0 - d, blk + 2*d) for b0 in blocks])
    Kq1_b = rho_b + Fs_b[:, None]; Kq2_b = Kq1_b**2 + Z_b + rho_b**2 + Fs_b[:, None]*Kq1_b
    nk = len(kidx)
    Hm = np.zeros((nk, Ns + 1), complex); Km = np.zeros_like(Hm)
    Erow = np.zeros(nk); ErK = np.zeros(nk); ErF = np.zeros(nk)
    D2r = np.zeros((nk, Ns + 1)); DS0r = np.zeros((nk, Ns + 1)); HS0r = np.zeros((nk, Ns + 1))
    tails = -1; t0 = time.time(); c = ATF*DN
    for r, i0 in enumerate(kidx):
        pr = profile(L, P, i0, wk, beta)
        if wk == 's': pr1 = profile(L, P, i0, 'exp', 0.0)
        sst = sig/kapf; mcell = np.floor(sst/ds).astype(int); th = sst/ds - mcell
        assert mcell[-1] + 1 + int(math.ceil(dsig/kapf/ds)) + 1 <= pr['J'] - 1, 'profile too short'
        h, hp, Eh, Ehp, H0, H1, H2 = pr['h'], pr['hp'], pr['Eh'], pr['Ehp'], pr['H0'], pr['H1'], pr['H2']
        Hv = (1 - th)*h[mcell] + th*h[mcell + 1]
        EHv = np.maximum(Eh[mcell], Eh[mcell + 1]) + th*(1 - th)/2*ds**2*H2[mcell] + 1e-12*H1[mcell]*ds
        if wk == 'exp':
            dkh = (hp + beta*h - pr['X0']*pr['R'])/kap
            Edkh = (Ehp + beta*Eh)/kap + 1e-12*abs(pr['X0']*pr['R'])
            h3 = kap**3*H0 + kap**2*pr['Zc']*H1 + kap*beta*pr['wmax']/(2*np.pi)
            dkh2 = (h3 + beta*H2 + abs(pr['X0'])*kap**2*pr['Zc']*pr['RA'])/kap
        else:
            dkh = (hp - pr1['h'])/kap
            Edkh = (Ehp + pr1['Eh'])/kap
            h3 = kap**3*H0 + kap**2*pr['Zc']*H1 + kap/(2*np.pi)
            dkh2 = (h3 + pr1['H2'])/kap
        Kv = (1 - th)*dkh[mcell] + th*dkh[mcell + 1]
        EKv = np.maximum(Edkh[mcell], Edkh[mcell + 1]) + th*(1 - th)/2*ds**2*dkh2[mcell]
        Hm[r] = tw*Hv; Km[r] = tw*Kv
        nsc = int(math.ceil(dsig/kapf/ds)) + 1
        idxc = np.minimum(mcell[:, None] + np.arange(nsc)[None, :], pr['J'] - 1)
        H0s, H1s, H2s = H0[idxc].max(1), H1[idxc].max(1), H2[idxc].max(1)
        trap = dsig**3/12*OmA*(H2s/kapf**2 + 2*H1s/kapf*rho_w + H0s*Z_w)
        H0u, H1u = unif(P, wk, beta, np.array([sigmax/kapf]))
        tailH = (H0u[0]*math.exp(-c*sigmax)/c) if wk == 'exp' else \
                ((2*IKPI/kap/kapf*(sigmax/c + 1/c**2) + 2*IKPIm/kap**2/c)*math.exp(-c*sigmax))
        Erow[r] = (np.sum(tw*EHv*OmA) + np.sum(trap[:-1]) + tailH + L['rel_a']*np.sum(tw*abs(Hv)*OmA))*UP
        DK0, DK1, H0c, H1c, H2c = dk_bounds(P, wk, beta, pr, 0.0)
        dK0s = DK0[idxc].max(1); dK1s = DK1[idxc].max(1); dK2s = dkh2[idxc].max(1)
        trapK = dsig**3/12*OmA*(dK2s/kapf**2 + 2*dK1s/kapf*rho_w + dK0s*Z_w)
        tailK = (2*IKPI1 + 2*IKPI/kap + PI_ + 2*IKPI1m)/kap*(1 + sigmax/kapf + 1/(c*kapf))*math.exp(-c*sigmax)/c
        ErK[r] = (np.sum(tw*EKv*OmA) + np.sum(trapK[:-1]) + tailK)*UP
        trapF = dsig**3/12*OmA*(H2s/kapf**2*rho_w + 2*H1s/kapf*Z_w + H0s*(1 + Z_w*rho_w))
        ErF[r] = (np.sum(tw*EHv*OmA*rho_w) + np.sum(trapF[:-1]) + tailH*rho_w.max())*UP
        DKd0, DKd1, H0cd, H1cd, H2cd = dk_bounds(P, wk, beta, pr, d)
        HS0 = H0cd[idxc].max(1); HS1 = H1cd[idxc].max(1); HS2 = H2cd[idxc].max(1); DS0 = DKd0[idxc].max(1)
        sK = P.kap_lo*sig
        if wk == 'exp':
            D2 = (HS2 + 2*beta*HS1 + beta**2*HS0 + 3*kap*T('KPi1', sK) + beta*T('KPi', sK))/kap**2
        else:
            D2 = (HS2 + 2*2*IKPI1 + T('KPi', sK))/kap**2
        D2r[r] = dsig*OmA*D2; DS0r[r] = dsig*OmA*DS0; HS0r[r] = dsig*OmA*HS0
        def kernL(l):
            cum = np.concatenate([[0], np.cumsum(vprof)*(1/256)])
            k = np.minimum(np.floor(l*256).astype(int), len(vprof))
            extra = np.where(k >= len(vprof), (l - len(vprof)/256)*vprof[-1], 0.0)
            return np.exp(-(cum[k] + extra)*DN)*(1 + 1e-12)
        tl = ibp_bound(P, wk, beta, pr, d, kernL)
        tr = ibp_bound(P, wk, beta, pr, d, lambda l: np.exp(-U_RIGHT30*l*DN))
        tails = max(tails, tl, tr)
        if verbose and r % 100 == 0: print('   row %d/%d t_k=%.4f  Erow %.2e  tails %.4f  (%.0fs)' % (r, nk, i0*DT, Erow[r], max(tl, tr), time.time() - t0), flush=True)
    # second-derivative suprema per (row, t_f block); the sigma-tail (sigma > sigmax) adds a tiny explicit term
    tailS = 2e4*(sigmax/c + 1/c**2)*math.exp(-c*sigmax)
    Skk = (D2r @ np.broadcast_to(Fs_b, (Ns + 1, len(blocks))) + tailS)*UP
    Skf = (DS0r @ (Fs_b[:, None]*Kq1_b).T + tailS)*UP
    Sff = (HS0r @ (Fs_b[:, None]*Kq2_b).T + tailS)*UP
    SO = 0.5*d*d*(Skk + 2*Skf + Sff)                                    # (nk, nblocks)
    bidx = np.minimum(np.floor((tf_idx*DT - blocks[0] + 1e-12)/blk).astype(int), len(blocks) - 1)
    # cells at a block edge overlap the neighbouring block: use the max of both
    bidx2 = np.minimum(np.floor((tf_idx*DT + d - blocks[0])/blk).astype(int), len(blocks) - 1)
    bidx1 = np.maximum(np.floor((tf_idx*DT - d - blocks[0])/blk).astype(int), 0)
    SOc = np.maximum(SO[:, bidx1], SO[:, bidx2])
    B = 800; Mmax = -1; arg = None; fmax = -1; farg = None; LOW = [-1, None]
    for c0 in range(0, len(tf_idx), B):
        cs = slice(c0, c0 + B)
        rws = rows[:, cs]
        Aq = L['a'][im + rws]; A1q = L['a1'][im + rws]
        I0 = Hm @ Aq; I1 = Hm @ A1q; IK = Km @ Aq
        ps = psi[cs]; Fc = Fv[cs]; aF = abs(Fc) + rF[cs]
        f = I0/ps; fK = IK/ps; fF = E*(I1 + Fc*I0)/ps
        Ef = aF[None, :]*Erow[:, None] + abs(f)*relpsi[cs][None, :]
        EK = aF[None, :]*ErK[:, None] + abs(fK)*relpsi[cs][None, :]
        EF = aF[None, :]*(ErF[:, None] + aF[None, :]*Erow[:, None]) + abs(fF)*relpsi[cs][None, :]
        cell = (abs(f) + Ef + d*(abs(fK) + EK) + d*(abs(fF) + EF) + SOc[:, cs])*UP
        i, jj = np.unravel_index(np.argmax(cell), cell.shape)
        if cell[i, jj] > Mmax:
            Mmax = float(cell[i, jj]); arg = dict(t_k=kidx[i]*DT, t_f=tf_idx[cs][jj]*DT, absf=float(abs(f[i, jj])), Ef=float(Ef[i, jj]),
                                                  grad=float(d*(abs(fK[i, jj]) + abs(fF[i, jj]))), second=float(SOc[i, c0 + jj]))
        lowb = abs(f) - Ef
        i3, j3 = np.unravel_index(np.argmax(lowb), f.shape)
        if lowb[i3, j3] > LOW[0]: LOW[0] = float(lowb[i3, j3]); LOW[1] = (float(kidx[i3]*DT), float(tf_idx[cs][j3]*DT))
        i2, j2 = np.unravel_index(np.argmax(abs(f)), f.shape)
        if abs(f[i2, j2]) > fmax: fmax = float(abs(f[i2, j2])); farg = (kidx[i2]*DT, tf_idx[cs][j2]*DT)
    return dict(M_cells=Mmax, argmax=arg, fmax_grid=fmax, fmax_at=farg, M_lower=LOW[0], M_lower_at=LOW[1], second_order_max=float(SO.max()),
                tails_tf=float(tails), Erow_max=float(Erow.max()), rows=int(nk), cols=int(len(tf_idx)))

# ---------------- region (3) ----------------
def bregion(L, P, wk, beta, kidx, d):
    best = -1; arg = None; vals = []
    for i0 in kidx:
        pr = profile(L, P, i0, wk, beta)
        b = ibp_bound(P, wk, beta, pr, d, lambda l: T('OmPsi', l))
        vals.append(b)
        if b > best: best, arg = b, i0*DT
    return best, arg, np.array(vals)

# ---------------- region (4): |t_k| >= TB, quasi-static representation ----------------
def farfield(P, wk, beta, TB, side, T1=None, ds=0.01, U=400.0):
    """upper bound of |h(0)| + int |h'| Omega_psi(kf s) ds for all t_k with t_k >= TB (side 'right') or t_k <= -TB ('left').
    h = h0 + (1/c0) int G h0'' du + (1/c0) R W0(h0), c0 = kap/(2 pi), h0 = c0 w/Q, Q = kap^2 e^2 z,
    |G(s,u)| <= K_Pi(kap|s-u|), |d_s G| <= kap K_Pi1(kap|s-u|), |R'(s)X'(0)| <= kap^2 |rho(tau_s)| K_Pi1(kap s),
    |rho| <= r + eps_D (Lemma D), |z(tau)|^2 = tau^2 + eta^2.  Evaluated at |t_k| = TB (nonincreasing in |t_k|)."""
    from cl import tail_side
    kap, kapf = P.kap, P.kapf_lo
    c0 = kap/(2*np.pi)
    s = np.arange(0, U, ds); slo, shi = s, s + ds
    tau_lo = (TB + kap*slo) if side == 'right' else (-TB + kap*slo)
    tau_hi = tau_lo + kap*ds
    # min |z| on each s-cell
    tmin = np.where((tau_lo <= 0) & (tau_hi >= 0), 0.0, np.minimum(abs(tau_lo), abs(tau_hi)))
    zmin = np.sqrt(tmin**2 + ETA**2)*DN
    if wk == 'exp':
        wmax = np.exp(-beta*slo); w1 = beta*wmax; w2 = beta*beta*wmax; w0, w10 = 1.0, beta
    else:
        wmax = shi; w1 = np.ones_like(s); w2 = 0*s; w0, w10 = 0.0, 1.0
    h0pp = (w2/zmin + 2*kap*w1/zmin**2 + 2*kap**2*wmax/zmin**3)/(2*np.pi*kap)
    h0p = w1/(2*np.pi*kap*zmin) + wmax/(2*np.pi*zmin**2)          # |w'|/(2 pi kap |z|) + w/(2 pi |z|^2)
    z0 = math.sqrt(TB**2 + ETA**2)
    h00 = w0/(2*np.pi*kap*z0); h0p0 = w10/(2*np.pi*kap*z0) + w0/(2*np.pi*z0**2)
    # rho bound on the s-cells (Lemma D): |Ai'/Ai| <= r + eps_D for |tau| >= 20; inside, the table bound would be needed:
    # here every s-cell used with the boundary term has |tau| >= T1 >= 20 (left) or tau >= TB (right).
    eps = max(tail_side(20, '0.05', 'left')[4], tail_side(20, '0.05', 'right')[4])
    rhob = np.sqrt(np.sqrt(tau_hi**2 + ETA**2) if side == 'right' else np.sqrt(tau_lo**2 + ETA**2)) + eps
    # convolution int kap K_Pi1(kap|s-u|) |h0''(u)| du over s-cells (lower separation)
    n = len(s); m = np.arange(n)
    KP1 = T('KPi1', P.kap_lo*np.maximum(m - 1, 0)*ds); KP = T('KPi', P.kap_lo*np.maximum(m - 1, 0)*ds)
    conv1 = np.convolve(h0pp, np.r_[KP1[::-1], KP1[1:]])[n - 1:2*n - 1]*ds*kap
    # tail u > U: |h0''| <= its value at U (decreasing for |tau| large): crude
    conv1 += kap*PI1*h0pp[-1]*s[-1]          # u > U: |h0''| <= h0''(U) (U/u)^2
    hp_s = h0p + conv1/c0 + (kap*T('KPi1', P.kap_lo*slo)*h0p0 + kap**2*rhob*T('KPi1', P.kap_lo*slo)*h00)/c0
    if side == 'left':
        if T1 is None: T1 = TB/2
        s1 = (TB - T1)/kap
        H0u, H1u = unif(P, wk, beta, shi)
        hp_s = np.where(shi <= s1, hp_s, H1u)
    body = np.sum(hp_s*T('OmPsi', kapf*slo))*ds
    cA = AF*kapf; Ue = s[-1] + ds
    H0u, H1u = unif(P, wk, beta, np.array([Ue]))
    tailv = (2*IKPI1*math.exp(-cA*Ue)/cA) if wk == 'exp' else ((2*IKPI1*(Ue/cA + 1/cA**2) + 2*IKPI1m/kap/cA)*math.exp(-cA*Ue))
    hzero = h00 + (np.sum(T('KPi', P.kap_lo*slo)*h0pp)*ds + PI_*h0pp[-1]*s[-1])/c0 + (PI_*h0p0 + kap*PI1*h00)/c0
    return float((hzero + body + tailv)*UP)

def certify(pair, wk, beta, linef, ftabf, thr, TB=80.0, verbose=True):
    """thr: rows whose region-3 bound exceeds thr are treated as strip rows (any thr is valid; it only affects sharpness)."""
    t0 = time.time()
    L, Ft = load_tables(linef, ftabf); P = Pair(pair)
    # region (3): coarse grid 1/16 (half-width 1/32) on [-TB, TB]; cells above 0.8 thr are recomputed on the grid 1/64
    kc = np.arange(int(-TB*128), int(TB*128) + 1, 8)
    bc, _, vc = bregion(L, P, wk, beta, kc, 4*DT)
    fine = sorted(set(int(i) for c in kc[vc > 0.8*thr] for i in range(c - 6, c + 7, 2) if abs(i) <= TB*128))
    if fine:
        bf, _, vf = bregion(L, P, wk, beta, np.array(fine), DT)
    else:
        vf = np.array([])
    keep_c = vc <= 0.8*thr
    b3 = max(vc[keep_c].max() if keep_c.any() else -1, vf[vf <= thr].max() if (vf <= thr).any() else -1)
    sel = np.array(fine)[vf > thr] if len(fine) else np.array([], int)
    if verbose: print('  region 3: %d coarse + %d fine rows, max outside strip %.5f, strip rows %d (%.0fs)'
                      % (len(kc), len(fine), b3, len(sel), time.time() - t0), flush=True)
    rows = np.array(sorted(set(int(i) + o for i in sel for o in (-1, 0, 1))))
    st = strip(L, Ft, P, wk, beta, rows, verbose=False) if len(rows) else dict(M_cells=-1, tails_tf=-1)
    if verbose: print('  strip: %s (%.0fs)' % ({k: st[k] for k in ('M_cells', 'fmax_grid', 'fmax_at', 'tails_tf', 'second_order_max')}, time.time() - t0), flush=True)
    fl = farfield(P, wk, beta, TB, 'left'); fr = farfield(P, wk, beta, TB, 'right')
    M = max(st['M_cells'], st['tails_tf'], b3, fl, fr)
    res = dict(pair=pair, weight=wk, beta=beta, M=M, regions=dict(strip_cells=st['M_cells'], strip_tf_tails=st['tails_tf'],
               region3=float(b3), far_left=fl, far_right=fr), fmax_grid=st.get('fmax_grid'), fmax_at=st.get('fmax_at'),
               argmax=st.get('argmax'), M_lower=st.get('M_lower'), M_lower_at=st.get('M_lower_at'), strip_tk=[float(rows.min()*DT), float(rows.max()*DT)] if len(rows) else None,
               strip_rows=int(len(rows)), TB=TB, thr=thr, seconds=time.time() - t0)
    return res
if __name__ == '__main__':
    pair, wk, beta = sys.argv[1], sys.argv[2], float(sys.argv[3]); thr = float(sys.argv[4])
    linef, ftabf = sys.argv[5], sys.argv[6]
    r = certify(pair, wk, beta, linef, ftabf, thr)
    print('RESULT', json.dumps(r, default=float), flush=True)
