"""Tail lemma (sigma >= sigma_h): uniform-in-beta bounds for crit2 with lambda=1, variant (ii), route 1, hi_sup.
For S >= S_h = 2 sigma_h^{1/3} and 0 <= beta <= beta_A(S) = sqrt(S/(2 c_A)), every 'bad' quantity is bounded by an
expression that is (a) free of beta and (b) nonincreasing in S on [S_h, inf); we evaluate it at S = S_h.
See report, section 4.  usage: tail.py N J G cA sigma_h"""
import sys
from arbprim import *
from arbcrit import K, params, Pbound

Y2SINH = arb(6).sqrt()/2          # sup_{y>=0} y^2/sinh(y) <= sup y/(1+y^2/6) = sqrt(6)/2

def tail_margin(P, cA, sigma_h, th_rho='0.5', coupling=1):
    a, at, ct, mF, Pi, Pi1 = K['a'], K['at'], K['ct'], K['mF'], K['Pi'], K['Pi1']
    g, nu1, c, k12, k21 = P['g'], P['nu1'], P['c'], coupling*P['k12'], coupling*P['k21']
    cA = arb(cA); sig = arb(sigma_h); S = 2*sig**(ONE/3); th = ONE/3
    bA = (S/(2*cA)).sqrt()
    r1, rg, a1, ag = a*nu1**(-th), a*g**(-th), at*nu1**(-th), at*g**(-th)
    C1, Cg = 2*PI*nu1**(-2*th)*sig**(-2*th), 2*PI*g**(-2*th)*sig**(-2*th)
    # monotonicity conditions used (decreasing in S on [S_h, inf)):
    for r in (r1, rg):
        assert (S*r).lower() > 2      # (bA(S)+r) e^{-rS/2}, (bA+r) S e^{-rS}, S e^{-rS} decreasing
    rho = arb(th_rho)*a1
    assert (rho - r1).upper() < 0 and (rho - rg).upper() < 0
    I = lambda av: 1/(av - rho) + 1/(av + rho)
    k0 = {1: C1*Pi, 'g': Cg*Pi}; k1 = {1: C1*Pi1*nu1**(-th), 'g': Cg*Pi1*g**(-th)}
    N1, Ng = k0[1]*I(a1), k0['g']*I(ag)
    tau2 = k12*k21*N1*Ng; assert tau2.upper() < 1
    u1, u0 = k0['g']*I(ag), k0['g']/(ag + a1)
    B_1 = k21*u1/(1 - tau2); B20 = k12*k21*k0['g']/(ag + rho)*N1*B_1; A0 = k12*k0[1]/(a1 + rho)*B_1
    P1 = A0 + (2 - c)*(k21*u0 + B20)
    Ub = k0['g']*(1 + ((rho - r1)*S).exp()) + k1['g']*(1/(ag - rho) + 1/(ag + r1))
    Bp = k21*Ub/(1 - tau2); Ap = k12*N1*Bp; Q1 = (-rho*S).exp()*(Ap + (2 - c)*Bp)
    v1, v0 = k0[1]*I(a1), k0[1]/(a1 + ag)
    A_g = k12*v1/(1 - tau2); A20 = k12*k21*k0[1]/(a1 + rho)*Ng*A_g; B0 = k21*k0['g']/(ag + rho)*A_g
    Pg = k12*v0 + A20 + (2 - c)*B0
    Vb = k0[1]*(1 + ((rho - rg)*S).exp()) + k1[1]*(1/(a1 - rho) + 1/(a1 + rg))
    Agp = k12*Vb/(1 - tau2); Bgp = k21*Ng*Agp; Qg = (-rho*S).exp()*(Agp + (2 - c)*Bgp)
    assert P1.upper() < 1
    epx = ((2 - c)*P1 + Pg)/((1 - P1)*(2 - c)); assert epx.upper() < 1
    xb1 = 1 + epx
    hi1 = Ap + c*tau2*Bp          # higher-order coefficient for column 1 (times int w e^{-rho s})
    hig = tau2*Agp + c*Bgp
    mx = lambda x, y: amax(x, y)
    # (beta+r) * Mo(w_m)   [uniform in beta]
    J0 = lambda rf, ak: 1/ak + 1/(rf + ak)                  # (beta+rf) * J(beta; rf, ak) <= this
    BMo1m = c*k21*(Cg*Pi*(mx(ONE, r1/ag) + (bA + r1)*(-r1*S).exp()/ag) + Cg*Pi1*g**(-th)*J0(r1, ag)) + hi1*mx(ONE, r1/rho)
    BMogm = k12*(C1*Pi*(mx(ONE, rg/a1) + (bA + rg)*(-rg*S).exp()/a1) + C1*Pi1*nu1**(-th)*J0(rg, a1)) + hig*mx(ONE, rg/rho)
    # Mo(w_m) itself (for y_m), beta-uniform: beta=0 is the worst case of each term
    Jb0 = lambda rf, ak: 1/(rf*ak) + 1/(rf*(rf + ak))
    Mo1m = c*k21*(Cg*Pi*(1/ag + (-r1*S).exp()/ag) + Cg*Pi1*g**(-th)*Jb0(r1, ag)) + hi1/rho
    # Mo(w_hat) (normalised lambda=1 weight, w_hat <= 2s/S and <= 2)
    J1 = lambda rf, ak: 2/(ak*rf**2) + 1/(ak**2*rf)
    Mo1h = c*k21*(Cg*Pi*(2/(S*ag**2) + 2*(-r1*S).exp()/ag) + Cg*Pi1*g**(-th)*(2/S)*J1(r1, ag)) + hi1*2/(S*rho**2)
    Mogh = k12*(C1*Pi*(2/(S*a1**2) + 2*(-rg*S).exp()/a1) + C1*Pi1*nu1**(-th)*(2/S)*J1(rg, a1)) + hig*2/(S*rho**2)
    # (beta+r) Dl/e2  and denominators (beta+r) kappa_m
    BDl = lambda r: 2/S + 2*Y2SINH/(r*S**2) + (bA + r)*(2*(-r*S).exp() + (-r*S/2).exp())
    den = lambda r: ct*r - 2*(bA + r)*(-r*S).exp()
    d1, dg = den(r1), den(rg); assert d1.lower() > 0 and dg.lower() > 0
    th_max = 1 - mx(BDl(r1)/d1, BDl(rg)/dg)
    eW = mx(((bA + r1)*epx*2 + xb1*BMo1m)/d1, BMogm/(c*dg))
    eD = mx((epx*BDl(r1) + xb1*(bA + r1)*Mo1h)/d1, (bA + rg)*Mogh/(c*dg))
    thK = th_max - eW - eD
    assert thK.lower() > 0 and eW.upper() < 1
    # leakage
    Sp1 = float((S*nu1**(-th)).lower()); Spg = float((S*g**(-th)).lower())
    P_1, P_g = Pbound(Sp1), Pbound(Spg)
    ph1 = nu1**th/mF
    E4 = mx(xb1*(nu1**(-th)*P_1 + Q1)*(bA + r1)/d1, ((2 - c)*g**(-th)*P_g + Qg)/c*(bA + rg)/dg)
    e3 = ph1*(nu1**(-th)*P_1 + Q1)
    y_m = 2 + Mo1m; y_l = BDl(r1)/r1 + Mo1h
    dR = (1 - th_max + eD)/(1 - eW)
    Y = y_m*dR + y_l; Z2 = 2*y_m + y_l + y_m*dR + dR*y_l; Z3 = 2*y_m**2 + 2*y_m*y_l + y_l**2
    q = 1 - P1
    leak = 2*E4*ph1*Y/(q*thK) + e3**2/q**2 + 2*E4*e3*ph1*Z2/(q**2*thK) + E4**2*ph1**2*Z3/(q**2*thK*(1 - eW))
    m = amin(thK, 1 - leak)
    return dict(margin=float(m.lower()), th_max=th_max, eW=eW, eD=eD, leak=leak, S=S, bA=bA)

if __name__ == '__main__':
    N, J, G, cA, sh = sys.argv[1:6]; cp = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    o = tail_margin(params(N, J, G), cA, sh, coupling=cp)
    print('(%s,%s,%s) tail sigma>=%s (S_h=%s, bA(S_h)=%s): margin >= %.4f  th_max>=%s eW<=%s eD<=%s leak<=%s'
          % (N, J, G, sh, o['S'].str(5), o['bA'].str(4), o['margin'], lo(o['th_max']).str(5), hi(o['eW']).str(4), hi(o['eD']).str(4), hi(o['leak']).str(3)))
