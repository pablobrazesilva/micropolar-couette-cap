"""Interval-safe primitives (Arb balls) for the G1' coupling criterion.
Every primitive returns an enclosure valid for ALL parameter values in the input balls.
Monotone primitives are evaluated at corner POINTS (no dependency blow-up)."""
from flint import arb, ctx
ctx.prec = 80
ZERO = arb(0); ONE = arb(1); TWO = arb(2); PI = arb.pi()

def lo(x):  return arb(x.lower())
def hi(x):  return arb(x.upper())
def iv(a, b): return arb(a).union(arb(b))
def amax(x, y):  # interval max
    return iv(max(lo(x), lo(y)), max(hi(x), hi(y))) if True else None
def amin(x, y):
    return iv(min(lo(x), lo(y)), min(hi(x), hi(y)))
def pos(x):      # clamp to >= 0 (used where the true quantity is >= 0)
    return iv(max(lo(x), ZERO), max(hi(x), ZERO))

def _hp(f):
    """evaluate f at higher precision (point evaluations with cancellation)"""
    old = ctx.prec; ctx.prec = 256
    try: return f()
    finally: ctx.prec = old

# ---------- one-argument monotone functions ----------
def _psi_pt(y):  # (e^y - 1)/y, entire, increasing
    if y == 0: return ONE
    return _hp(lambda: y.expm1()/y)
def psi(y): return _psi_pt(lo(y)).union(_psi_pt(hi(y)))
def _shc_pt(y):  # sinh(y)/y, even, increasing in |y|
    if y == 0: return ONE
    return _hp(lambda: y.sinh()/y)
def shc(y):
    a, b = abs(lo(y)), abs(hi(y))
    m = ZERO if (lo(y) <= 0 <= hi(y)) else min(a, b)
    return _shc_pt(m).union(_shc_pt(max(a, b)))

# ---------- E(q,L) = int_0^L e^{q u} du : increasing in q and in L (L>=0) ----------
def _E_pt(q, L):
    if L <= 0: return ZERO
    return L*_psi_pt(q*L)
def E(q, L):
    L = pos(L)
    return _E_pt(lo(q), lo(L)).union(_E_pt(hi(q), hi(L)))

# ---------- E1(q,L) = int_0^L v e^{q v} dv : increasing in q, L ----------
def _E1_pt(q, L):
    if L <= 0: return ZERO
    def f():
        y = q*L
        if abs(y) < arb('1e-3'):
            s = ZERO; fact = ONE
            for n in range(0, 10):
                if n > 0: fact = fact*n
                s += q**n*L**(n + 2)/(fact*(n + 2))
            # remainder: sum_{n>=10} |y|^n L^2/(n!(n+2)) <= L^2 |y|^10/10! * e^{|y|}
            rem = L*L*abs(y)**10/arb(3628800)*abs(y).exp()
            return s + arb(0, rem.upper())
        return (y.exp()*(y - 1) + 1)/(q*q)
    return _hp(f)
def E1(q, L):
    L = pos(L)
    return _E1_pt(lo(q), lo(L)).union(_E1_pt(hi(q), hi(L)))

# ---------- SHE(b,q,T,L) = int_0^L sh(b,T+v) e^{qv} dv,  sh(b,x)=sinh(bx)/b ----------
# increasing in b>=0, q, T>=0, L>=0
def _SHE0(q, T, L): return T*_E_pt(q, L) + _E1_pt(q, L)
def _SHE_pt(b, q, T, L):
    if L <= 0: return ZERO
    if b == 0: return _SHE0(q, T, L)
    if b*(T + L) <= arb('1e-3'):
        base = _SHE0(q, T, L)
        return base.union(base*_shc_pt(b*(T + L)))
    return _hp(lambda: ((b*T).exp()*_E_pt(q + b, L) - (-b*T).exp()*_E_pt(q - b, L))/(2*b))
def SHE(b, q, T, L):
    b = pos(b); T = pos(T); L = pos(L)
    return _SHE_pt(lo(b), lo(q), lo(T), lo(L)).union(_SHE_pt(hi(b), hi(q), hi(T), hi(L)))

# ---------- SHX(b,q,T) = int_0^T sh(b,T-v) e^{qv} dv : increasing in b, q, T ----------
def _SHX0(q, T): return T*_E_pt(q, T) - _E1_pt(q, T)
def _SHX_pt(b, q, T):
    if T <= 0: return ZERO
    if b == 0: return _SHX0(q, T)
    if b*T <= arb('1e-3'):
        base = _SHX0(q, T)
        return base.union(base*_shc_pt(b*T))
    return _hp(lambda: ((b*T).exp()*_E_pt(q - b, T) - (-b*T).exp()*_E_pt(q + b, T))/(2*b))
def SHX(b, q, T):
    b = pos(b); T = pos(T)
    return _SHX_pt(lo(b), lo(q), lo(T)).union(_SHX_pt(hi(b), hi(q), hi(T)))

def expi(x): return x.exp()
