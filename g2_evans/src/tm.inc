/* tm.inc -- Taylor models in the box variables t = (t1,t2,t3) in [-1,1]^3 (total degree <= TK) for the Evans function,
 * s-formulation (s in [0,S]):   Y_s = A(s) Y,  A affine in (Lam, b2 = beta^2, e = sigma^{-2/3});  S = 2 sigma^{1/3} = 2/sqrt(e).
 *   Lam = Lam_c + u t1  (u complex: half of the contour segment),   b2 = b2_c + rb t2,   S = S_c + rS t3,  e = 4/S^2.
 * A TM is a polynomial P(t) (acb coefficients) plus a remainder disk radius r:  f(t) in P(t) + D(0,r) for all t in the box.
 * Arithmetic: products truncated at degree TK (truncated monomials bounded by |c| since |t^k| <= 1).
 * Integration: Taylor steps on [0,S_c] (point step length h), then ONE step of TM length delta = rS t3 to reach S.
 * Because the dependence on the parameters is carried symbolically, the enclosure does not suffer the "wrapping"
 * of ball arithmetic in (Lam, alpha, sigma); only the degree > TK terms and the step tails go to the remainder. */
static int TK = 3, NM = 0;
static int MI[128][3];
static int MULI[128][128];
static int MONO_T1 = -1, MONO_T2 = -1, MONO_T3 = -1;
static int TKP = 99;
static double MSHIFT = 0.0;   /* growth normalisation (experimental; loses precision for large beta, see report) */      /* max degree in the parameter variables (t2,t3); total degree <= TK */
static int mono_in(int a, int b, int c) { return a + b + c <= TK && b + c <= TKP; }
static void tm_setup2(int K, int Kp)
{
    TK = K; TKP = Kp; NM = 0;
    for (int d = 0; d <= K; d++) for (int a = d; a >= 0; a--) for (int b = d - a; b >= 0; b--) { int c = d - a - b;
        if (!mono_in(a, b, c)) continue; MI[NM][0] = a; MI[NM][1] = b; MI[NM][2] = c; NM++; }
    for (int i = 0; i < NM; i++) for (int j = 0; j < NM; j++) {
        int a = MI[i][0] + MI[j][0], b = MI[i][1] + MI[j][1], c = MI[i][2] + MI[j][2]; MULI[i][j] = -1;
        if (mono_in(a, b, c)) for (int k = 0; k < NM; k++) if (MI[k][0] == a && MI[k][1] == b && MI[k][2] == c) { MULI[i][j] = k; break; }
    }
    for (int k = 0; k < NM; k++) { if (MI[k][0] == 1 && MI[k][1] == 0 && MI[k][2] == 0) MONO_T1 = k;
        if (MI[k][0] == 0 && MI[k][1] == 1 && MI[k][2] == 0) MONO_T2 = k; if (MI[k][0] == 0 && MI[k][1] == 0 && MI[k][2] == 1) MONO_T3 = k; }
}
static void tm_setup(int K) { tm_setup2(K, K); }
static double dabs_ub(const acb_t z) { arf_t u; arf_init(u); acb_get_abs_ubound_arf(u, z, 64); double d = arf_get_d(u, ARF_RND_UP); arf_clear(u); return d; }
/* sup over the box of |P| (crude, upper) */
static double tm_norm(acb_srcptr c, double r) { double s = r; for (int i = 0; i < NM; i++) if (!acb_is_zero(c + i)) s += dabs_ub(c + i); return s * (1 + 1e-12); }
static int tm_isconst(acb_srcptr c, double r) { if (r != 0) return 0; for (int i = 1; i < NM; i++) if (!acb_is_zero(c + i)) return 0; return 1; }
/* z += x*y */
static void tm_muladd(acb_ptr z, double *rz, acb_srcptr x, double rx, int xconst, acb_srcptr y, double ry)
{
    acb_t t; acb_init(t);
    if (xconst) {
        for (int j = 0; j < NM; j++) if (!acb_is_zero(y + j)) { acb_mul(t, x, y + j, PREC); acb_add(z + j, z + j, t, PREC); }
        *rz += dabs_ub(x) * ry * (1 + 1e-12);
    } else {
        double tr = 0;
        for (int i = 0; i < NM; i++) { if (acb_is_zero(x + i)) continue;
            for (int j = 0; j < NM; j++) { if (acb_is_zero(y + j)) continue;
                int k = MULI[i][j];
                if (k >= 0) { acb_mul(t, x + i, y + j, PREC); acb_add(z + k, z + k, t, PREC); }
                else tr += dabs_ub(x + i) * dabs_ub(y + j);
            } }
        double nx = tm_norm(x, 0), ny = tm_norm(y, 0);
        *rz += (tr + nx * ry + ny * rx + rx * ry) * (1 + 1e-12);
    }
    acb_clear(t);
}
/* sparse 20x20 matrix with TM entries */
typedef struct { int n; int row[400], col[400]; acb_ptr val; double r[400]; int cst[400]; } tsmat;
static void tsmat_init(tsmat *M) { M->n = 0; M->val = _acb_vec_init(400 * 128); }
static void tsmat_clear(tsmat *M) { _acb_vec_clear(M->val, 400 * 128); }
/* compound of a 6x6 TM matrix A (A[36*NM] coefficients, rA[36]) */
static void tm_compound3(tsmat *C, acb_srcptr A, const double *rA)
{
    acb_ptr D = _acb_vec_init(NW * NW * NM); double *rD = calloc(NW * NW, sizeof(double));
    for (int J = 0; J < NW; J++)
        for (int pos = 0; pos < 3; pos++)
            for (int m = 0; m < 6; m++) {
                int ent = m * 6 + TRIP[J][pos];
                int nz = rA[ent] != 0; for (int i = 0; i < NM && !nz; i++) if (!acb_is_zero(A + ent * NM + i)) nz = 1;
                if (!nz) continue;
                int t[3] = {TRIP[J][0], TRIP[J][1], TRIP[J][2]}; t[pos] = m;
                if (t[0] == t[1] || t[0] == t[2] || t[1] == t[2]) continue;
                int s = 1;
                for (int i = 0; i < 3; i++) for (int j = 0; j < 2 - i; j++) if (t[j] > t[j + 1]) { int x = t[j]; t[j] = t[j + 1]; t[j + 1] = x; s = -s; }
                int I = tidx(t[0], t[1], t[2]);
                for (int i = 0; i < NM; i++) { if (s > 0) acb_add(D + (I * NW + J) * NM + i, D + (I * NW + J) * NM + i, A + ent * NM + i, PREC);
                    else acb_sub(D + (I * NW + J) * NM + i, D + (I * NW + J) * NM + i, A + ent * NM + i, PREC); }
                rD[I * NW + J] += rA[ent];
            }
    C->n = 0;
    for (int I = 0; I < NW; I++) for (int J = 0; J < NW; J++) {
        int nz = rD[I * NW + J] != 0; for (int i = 0; i < NM && !nz; i++) if (!acb_is_zero(D + (I * NW + J) * NM + i)) nz = 1;
        if (!nz) continue;
        C->row[C->n] = I; C->col[C->n] = J; _acb_vec_set(C->val + C->n * NM, D + (I * NW + J) * NM, NM); C->r[C->n] = rD[I * NW + J];
        C->cst[C->n] = tm_isconst(C->val + C->n * NM, C->r[C->n]); C->n++; }
    _acb_vec_clear(D, NW * NW * NM); free(rD);
}
static double tsmat_norm(const tsmat *M)
{
    double rs[NW] = {0};
    for (int e = 0; e < M->n; e++) rs[M->row[e]] += tm_norm(M->val + e * NM, M->r[e]);
    double m = 0; for (int i = 0; i < NW; i++) if (rs[i] > m) m = rs[i]; return m * (1 + 1e-12);
}
/* out (20 TMs) += M * v */
static void tsmat_mulvec_add(acb_ptr out, double *rout, const tsmat *M, acb_srcptr v, const double *rv)
{
    for (int e = 0; e < M->n; e++)
        tm_muladd(out + M->row[e] * NM, rout + M->row[e], M->val + e * NM, M->r[e], M->cst[e], v + M->col[e] * NM, rv[M->col[e]]);
}
static double tmvec_maxnorm(acb_srcptr v, const double *rv) { double m = 0; for (int i = 0; i < NW; i++) { double d = tm_norm(v + i * NM, rv[i]); if (d > m) m = d; } return m; }

static long TM_STEPS = 0;
static double LAST_RPOLY = 0, LAST_RLOH = 0;
static acb_ptr LAST_DT = NULL; static double LAST_R = 0;
/* one Taylor step from Y (TM vector) with matrix B0 (TM), C1 (const TM matrix), step length h (point arb) or the TM step
 * delta = rS t3 (if tmstep).  Returns 0 if the majorant condition fails. */
static int tm_step(acb_ptr Y, double *rY, const tsmat *B0, const tsmat *C1, const arb_t h, int tmstep, double rS)
{
    int n20 = NW * NM;
    acb_ptr Ycur = _acb_vec_init(n20), Yprev = _acb_vec_init(n20), Ynew = _acb_vec_init(n20), acc = _acb_vec_init(n20);
    double rcur[NW], rprev[NW], rnew[NW], racc[NW];
    double hd = tmstep ? rS : arf_get_d(arb_midref(h), ARF_RND_UP) * (1 + 1e-15);
    double M0 = tsmat_norm(B0), M1 = tsmat_norm(C1);
    double th = (M0 * hd + M1 * hd * hd) / (DEG + 1);
    if (th >= 0.5) { _acb_vec_clear(Ycur, n20); _acb_vec_clear(Yprev, n20); _acb_vec_clear(Ynew, n20); _acb_vec_clear(acc, n20); return 0; }
    _acb_vec_set(Ycur, Y, n20); _acb_vec_zero(Yprev, n20); _acb_vec_set(acc, Y, n20);
    for (int i = 0; i < NW; i++) { rcur[i] = rY[i]; rprev[i] = 0; racc[i] = rY[i]; }
    arb_t hp; arb_init(hp); arb_one(hp); acb_t t; acb_init(t);
    acb_ptr dm = _acb_vec_init(NM);      /* delta^m as TM (monomial rS^m t3^m or remainder) */
    for (int m = 0; m < DEG; m++) {
        _acb_vec_zero(Ynew, n20); for (int i = 0; i < NW; i++) rnew[i] = 0;
        tsmat_mulvec_add(Ynew, rnew, B0, Ycur, rcur);
        if (m > 0) tsmat_mulvec_add(Ynew, rnew, C1, Yprev, rprev);
        _acb_vec_scalar_div_ui(Ynew, Ynew, n20, m + 1, PREC);
        for (int i = 0; i < NW; i++) rnew[i] = rnew[i] / (m + 1) * (1 + 1e-12);
        if (!tmstep) {
            arb_mul(hp, hp, h, PREC);
            for (int i = 0; i < n20; i++) { acb_mul_arb(t, Ynew + i, hp, PREC); acb_add(acc + i, acc + i, t, PREC); }
            double hpd = arf_get_d(arb_midref(hp), ARF_RND_UP) * (1 + 1e-12) + mag_get_d(arb_radref(hp));
            for (int i = 0; i < NW; i++) racc[i] += rnew[i] * hpd;
        } else {
            /* acc += Ynew * delta^{m+1},  delta = rS t3 */
            int p = m + 1; double rp = pow(rS, p) * (1 + 1e-12);
            for (int i = 0; i < NW; i++) {
                acb_srcptr y = Ynew + i * NM;
                for (int k = 0; k < NM; k++) { if (acb_is_zero(y + k)) continue;
                    int a = MI[k][0], b = MI[k][1], c = MI[k][2] + p, kk = -1;
                    if (mono_in(a, b, c)) for (int q = 0; q < NM; q++) if (MI[q][0] == a && MI[q][1] == b && MI[q][2] == c) { kk = q; break; }
                    if (kk >= 0) { acb_mul_2exp_si(t, y + k, 0); arb_t rr; arb_init(rr); arb_set_d(rr, rS); arb_pow_ui(rr, rr, p, PREC);
                        acb_mul_arb(t, t, rr, PREC); acb_add(acc + i * NM + kk, acc + i * NM + kk, t, PREC); arb_clear(rr); }
                    else racc[i] += dabs_ub(y + k) * rp;
                }
                racc[i] += rnew[i] * rp;
            }
        }
        acb_ptr sw = Yprev; Yprev = Ycur; Ycur = Ynew; Ynew = sw;
        for (int i = 0; i < NW; i++) { rprev[i] = rcur[i]; rcur[i] = rnew[i]; }
    }
    double cN = tmvec_maxnorm(Ycur, rcur) * pow(hd, DEG), cN1 = tmvec_maxnorm(Yprev, rprev) * pow(hd, DEG - 1);
    double tail = fmax(cN, cN1) * 2 * th / (1 - th) * (1 + 1e-10);
    for (int i = 0; i < NW; i++) racc[i] += tail;
    _acb_vec_set(Y, acc, n20); for (int i = 0; i < NW; i++) rY[i] = racc[i];
    _acb_vec_clear(Ycur, n20); _acb_vec_clear(Yprev, n20); _acb_vec_clear(Ynew, n20); _acb_vec_clear(acc, n20); _acb_vec_clear(dm, NM);
    arb_clear(hp); acb_clear(t); TM_STEPS++;
    return 1;
}

#include "lohner.inc"

/* TM enclosure of D over: Lam = Lc + u t1, b2 in [b2lo,b2hi], S in [Slo,Shi]  ->  acb rectangle E (incl. remainder) */
static int evans_tm(acb_t E, const params_t *P, const acb_t Lc, const acb_t u, double b2lo, double b2hi, double Slo, double Shi)
{
    double Sc = 0.5 * (Slo + Shi), rS = 0.5 * (Shi - Slo) * (1 + 1e-12) + 1e-300;
    double b2c = 0.5 * (b2lo + b2hi), rb = 0.5 * (b2hi - b2lo) * (1 + 1e-12);
    /* TMs of Lam, b2, e */
    acb_ptr TL = _acb_vec_init(NM), TB = _acb_vec_init(NM), TE = _acb_vec_init(NM); double rE = 0;
    acb_set(TL + 0, Lc); acb_set(TL + MONO_T1, u);
    arb_set_d(acb_realref(TB + 0), b2c); arb_set_d(acb_realref(TB + MONO_T2), rb);
    { /* e = 4/(Sc + rS t3)^2 = (4/Sc^2) sum_j (j+1)(-rho t3)^j ,  rho = rS/Sc */
      arb_t e0, rho, c; arb_init(e0); arb_init(rho); arb_init(c);
      arb_set_d(e0, Sc); arb_mul(e0, e0, e0, PREC); arb_inv(e0, e0, PREC); arb_mul_ui(e0, e0, 4, PREC);
      arb_set_d(rho, rS); arb_div_si(rho, rho, 1, PREC); { arb_t sc; arb_init(sc); arb_set_d(sc, Sc); arb_div(rho, rho, sc, PREC); arb_clear(sc); }
      int JE = TK < TKP ? TK : TKP;
      for (int j = 0; j <= JE; j++) {
          int k = -1; for (int q = 0; q < NM; q++) if (MI[q][0] == 0 && MI[q][1] == 0 && MI[q][2] == j) { k = q; break; }
          arb_pow_ui(c, rho, j, PREC); arb_mul_ui(c, c, j + 1, PREC); if (j & 1) arb_neg(c, c); arb_mul(c, c, e0, PREC);
          acb_set_arb(TE + k, c);
      }
      double rh = arf_get_d(arb_midref(rho), ARF_RND_UP) * (1 + 1e-12); double e0d = arf_get_d(arb_midref(e0), ARF_RND_UP) * (1 + 1e-12);
      if (rh >= 0.5) return 0;
      rE = e0d * pow(rh, JE + 1) * ((JE + 2) - (JE + 1) * rh) / ((1 - rh) * (1 - rh)) * (1 + 1e-10);
      arb_clear(e0); arb_clear(rho); arb_clear(c); }
    /* A(0) entries (TM) and A1 = dA/ds */
    acb_ptr A = _acb_vec_init(36 * NM), A1 = _acb_vec_init(36 * NM); double rA[36] = {0}, rA1[36] = {0};
    acb_t z; acb_init(z); arb_t t; arb_init(t);
#define AE(i, j) (A + ((i) * 6 + (j)) * NM)
    acb_one(AE(0, 1)); acb_one(AE(2, 3)); acb_one(AE(4, 5)); acb_one(AE(5, 0)); acb_set_arb(AE(5, 2), P->c);
    _acb_vec_set(AE(5, 4), TB, NM);
    /* q1(0) = (Lam + nu1 b2 + s1 e)/nu1 */
    for (int k = 0; k < NM; k++) {
        acb_set(z, TL + k); acb_mul_arb(AE(1, 0) + k, TB + k, P->nu1, PREC); acb_add(AE(1, 0) + k, AE(1, 0) + k, z, PREC);
        acb_mul_arb(z, TE + k, P->s1, PREC); acb_add(AE(1, 0) + k, AE(1, 0) + k, z, PREC); acb_div_arb(AE(1, 0) + k, AE(1, 0) + k, P->nu1, PREC);
        acb_set(z, TL + k); acb_mul_arb(AE(3, 2) + k, TB + k, P->g, PREC); acb_add(AE(3, 2) + k, AE(3, 2) + k, z, PREC);
        acb_mul_arb(z, TE + k, P->sg, PREC); acb_add(AE(3, 2) + k, AE(3, 2) + k, z, PREC); acb_div_arb(AE(3, 2) + k, AE(3, 2) + k, P->g, PREC);
        arb_div(t, P->k12, P->nu1, PREC); acb_mul_arb(AE(1, 2) + k, TE + k, t, PREC);
        arb_div(t, P->k21, P->g, PREC); acb_mul_arb(AE(3, 0) + k, TE + k, t, PREC);
    }
    arb_inv(t, P->nu1, PREC); rA[1 * 6 + 0] = rE * arf_get_d(arb_midref(P->s1), ARF_RND_UP) * arf_get_d(arb_midref(t), ARF_RND_UP) * 1.0001;
    arb_inv(t, P->g, PREC);   rA[3 * 6 + 2] = rE * fabs(arf_get_d(arb_midref(P->sg), ARF_RND_UP)) * arf_get_d(arb_midref(t), ARF_RND_UP) * 1.0001;
    arb_div(t, P->k12, P->nu1, PREC); rA[1 * 6 + 2] = rE * fabs(arf_get_d(arb_midref(t), ARF_RND_UP)) * 1.0001;
    arb_div(t, P->k21, P->g, PREC); rA[3 * 6 + 0] = rE * fabs(arf_get_d(arb_midref(t), ARF_RND_UP)) * 1.0001;
    arb_inv(t, P->nu1, PREC); arb_neg(t, t); acb_set_arb(z, t); acb_mul_onei(A1 + (1 * 6 + 0) * NM, z);
    arb_inv(t, P->g, PREC); arb_neg(t, t); acb_set_arb(z, t); acb_mul_onei(A1 + (3 * 6 + 2) * NM, z);
#undef AE
    tsmat C0, C1, B0; tsmat_init(&C0); tsmat_init(&C1); tsmat_init(&B0);
    tm_compound3(&C0, A, rA); tm_compound3(&C1, A1, rA1);
    /* growth normalisation: W~ = exp(-m s) W with m = MSHIFT*sqrt(beta^2) (TM in t2), i.e. C0 -> C0 - m I.
       D~ = b.W~(S) = exp(-m S) D: same zeros, same winding in Lambda (m does not depend on Lambda).  It removes the
       dominant parameter sensitivity of D ~ exp(3 beta S) for large beta. */
    /* C0L: matrix used for the Lohner BALL propagator.  With the growth normalisation W~ = e^{-m(p) s} W, m(p) = MSHIFT
       sqrt(beta^2) (a TM in t2), the exact step propagator is  e^{-(m(p)-m_c) h} * Phi_{m_c}(p),  Phi_{m_c} = propagator of
       A(p) - m_c I with the POINT constant m_c.  So the ball encloses Phi_{m_c} (no radius from m on the diagonal) and the
       scalar factor is bounded by e^{dm |h|}, dm = sup|m(p) - m_c| (applied to rho, see lohner_scale). */
    tsmat C0L; tsmat_init(&C0L); C0L.n = 0; double DM = 0;
    acb_t MC; acb_init(MC);
    if (MSHIFT > 0 && b2c > 2.5 * rb && b2c > 0) {
        acb_ptr TMm = _acb_vec_init(NM); double rM = 0;
        arb_t c0, x, bin; arb_init(c0); arb_init(x); arb_init(bin);
        arb_set_d(c0, b2c); arb_sqrt(c0, c0, PREC); arb_mul_si(c0, c0, 1, PREC);
        { arb_t ms; arb_init(ms); arb_set_d(ms, MSHIFT); arb_mul(c0, c0, ms, PREC); arb_clear(ms); }
        double rho2 = rb / b2c; int J = TKP < TK ? TKP : TK;
        arb_one(bin);                                            /* binom(1/2, j) */
        for (int j = 0; j <= J; j++) {
            int k = -1; for (int q = 0; q < NM; q++) if (MI[q][0] == 0 && MI[q][1] == j && MI[q][2] == 0) { k = q; break; }
            arb_set_d(x, rho2); arb_pow_ui(x, x, j, PREC); arb_mul(x, x, bin, PREC); arb_mul(x, x, c0, PREC);
            acb_set_arb(TMm + k, x);
            /* bin <- bin * (1/2 - j)/(j+1) */
            arb_set_d(x, 0.5 - j); arb_mul(bin, bin, x, PREC); arb_div_ui(bin, bin, j + 1, PREC);
        }
        rM = fabs(arf_get_d(arb_midref(bin), ARF_RND_UP)) * pow(rho2, J + 1) / (1 - rho2) * arf_get_d(arb_midref(c0), ARF_RND_UP) * 1.0001;
        for (int I0 = 0; I0 < NW; I0++) {
            int f = -1; for (int e = 0; e < C0.n; e++) if (C0.row[e] == I0 && C0.col[e] == I0) { f = e; break; }
            if (f < 0) { f = C0.n; C0.row[f] = I0; C0.col[f] = I0; _acb_vec_zero(C0.val + f * NM, NM); C0.r[f] = 0; C0.n++; }
            _acb_vec_sub(C0.val + f * NM, C0.val + f * NM, TMm, NM, PREC); C0.r[f] += rM;
            C0.cst[f] = tm_isconst(C0.val + f * NM, C0.r[f]);
        }
        acb_set(MC, TMm + 0); DM = (tm_norm(TMm, rM) - dabs_ub(TMm + 0) + 2 * mag_get_d(arb_radref(acb_realref(TMm + 0)))) * (1 + 1e-12);
        _acb_vec_clear(TMm, NM); arb_clear(c0); arb_clear(x); arb_clear(bin);
    }
    system_t Sy; boundary_init(&Sy, P);
    int n20 = NW * NM; acb_ptr Y = _acb_vec_init(n20); double rY[NW] = {0};
    for (int i = 0; i < NW; i++) acb_set(Y + i * NM, Sy.W0 + i);
    double Mn = tsmat_norm(&C0) + tsmat_norm(&C1) * Sc;
    int nsteps = (int) ceil(Mn * Sc / HSCALE); if (nsteps < 4) nsteps = 4;
    arb_t h, sk; arb_init(h); arb_init(sk); arb_set_d(h, Sc); arb_div_si(h, h, nsteps, PREC);
    int ok = 1;
    acb_mat_t Q; acb_mat_init(Q, NW, NW); acb_mat_one(Q); double rho[NW] = {0};
    smat Bb, C1b; smat_init(&Bb); smat_init(&C1b); acb_ptr Phi = _acb_vec_init(NW * NW);
    int useL = (DM > 0 || !acb_is_zero(MC));
    tsmat B0L; tsmat_init(&B0L);
    if (useL) {
        tm_compound3(&C0L, A, rA);
        for (int I0 = 0; I0 < NW; I0++) {
            int f = -1; for (int e = 0; e < C0L.n; e++) if (C0L.row[e] == I0 && C0L.col[e] == I0) { f = e; break; }
            if (f < 0) { f = C0L.n; C0L.row[f] = I0; C0L.col[f] = I0; _acb_vec_zero(C0L.val + f * NM, NM); C0L.r[f] = 0; C0L.n++; }
            acb_sub(C0L.val + f * NM, C0L.val + f * NM, MC, PREC); C0L.cst[f] = tm_isconst(C0L.val + f * NM, C0L.r[f]);
        }
    }
    for (int k = 0; k <= nsteps && ok; k++) {
        /* B0 = C0 + s_k C1,  s_k = k h ; C1 const */
        arb_mul_si(sk, h, k, PREC);
        B0.n = 0;
        for (int e = 0; e < C0.n; e++) { B0.row[B0.n] = C0.row[e]; B0.col[B0.n] = C0.col[e]; _acb_vec_set(B0.val + B0.n * NM, C0.val + e * NM, NM); B0.r[B0.n] = C0.r[e]; B0.n++; }
        for (int e = 0; e < C1.n; e++) {
            int f = -1; for (int q = 0; q < B0.n; q++) if (B0.row[q] == C1.row[e] && B0.col[q] == C1.col[e]) { f = q; break; }
            acb_mul_arb(z, C1.val + e * NM, sk, PREC);
            if (f < 0) { f = B0.n; B0.row[f] = C1.row[e]; B0.col[f] = C1.col[e]; _acb_vec_zero(B0.val + f * NM, NM); B0.r[f] = 0; B0.n++; }
            acb_add(B0.val + f * NM, B0.val + f * NM, z, PREC);
        }
        for (int e = 0; e < B0.n; e++) B0.cst[e] = tm_isconst(B0.val + e * NM, B0.r[e]);
        if (LOHNER) for (int i = 0; i < NW; i++) rY[i] = 0;        /* rY := local errors tau of this step */
        if (k < nsteps) ok = tm_step(Y, rY, &B0, &C1, h, 0, 0);
        else ok = tm_step(Y, rY, &B0, &C1, h, 1, rS);           /* final TM step of length delta = rS t3 */
        if (ok && LOHNER) {
            if (useL) {      /* B0L = C0L + s_k C1 */
                B0L.n = 0;
                for (int e = 0; e < C0L.n; e++) { B0L.row[B0L.n] = C0L.row[e]; B0L.col[B0L.n] = C0L.col[e]; _acb_vec_set(B0L.val + B0L.n * NM, C0L.val + e * NM, NM); B0L.r[B0L.n] = C0L.r[e]; B0L.n++; }
                for (int e = 0; e < C1.n; e++) {
                    int f = -1; for (int q = 0; q < B0L.n; q++) if (B0L.row[q] == C1.row[e] && B0L.col[q] == C1.col[e]) { f = q; break; }
                    acb_mul_arb(z, C1.val + e * NM, sk, PREC);
                    if (f < 0) { f = B0L.n; B0L.row[f] = C1.row[e]; B0L.col[f] = C1.col[e]; _acb_vec_zero(B0L.val + f * NM, NM); B0L.r[f] = 0; B0L.n++; }
                    acb_add(B0L.val + f * NM, B0L.val + f * NM, z, PREC);
                }
                tsmat_to_ball(&Bb, &B0L);
            } else tsmat_to_ball(&Bb, &B0);
            tsmat_to_ball(&C1b, &C1);
            ok = ball_prop(Phi, &Bb, &C1b, h, k == nsteps, rS);
            if (ok) {
                if (useL) {   /* scalar factor |e^{-(m-m_c) h}| <= e^{DM |h|}: scale Phi's bound through rho and tau */
                    double hh = (k == nsteps) ? rS : arf_get_d(arb_midref(h), ARF_RND_UP) * (1 + 1e-15);
                    double fac = exp(DM * hh) * (1 + 1e-12);
                    for (int i = 0; i < NW; i++) rho[i] *= fac;     /* rho' = fac |Q'^-1 Phi_mc Q| rho + |Q'^-1| tau */
                }
                lohner_update(Q, rho, Phi, rY);
            }
            for (int i = 0; i < NW; i++) rY[i] = 0;
        }
    }
    double rLoh = 0;
    if (ok && LOHNER) {     /* remainder of D:  sum_j |(b^T Q)_j| rho_j */
        acb_t bq, tt; acb_init(bq); acb_init(tt);
        for (int j = 0; j < NW; j++) {
            acb_zero(bq);
            for (int i = 0; i < NW; i++) { acb_mul(tt, Sy.b + i, acb_mat_entry(Q, i, j), PREC); acb_add(bq, bq, tt, PREC); }
            rLoh += acb_ub(bq) * rho[j];
        }
        rLoh *= (1 + 1e-12); acb_clear(bq); acb_clear(tt);
    }
    acb_mat_clear(Q); smat_clear(&Bb); smat_clear(&C1b); _acb_vec_clear(Phi, NW * NW);
    if (ok) {
        acb_ptr Dt = _acb_vec_init(NM); double rD = 0;
        for (int i = 0; i < NW; i++) tm_muladd(Dt, &rD, Sy.b + i, 0, 1, Y + i * NM, rY[i]);
        /* range of the polynomial over [-1,1]^3 */
        acb_zero(E); arb_t rng; arb_init(rng); acb_t w; acb_init(w);
        for (int k = 0; k < NM; k++) {
            if (acb_is_zero(Dt + k)) continue;
            int odd = (MI[k][0] & 1) || (MI[k][1] & 1) || (MI[k][2] & 1), deg = MI[k][0] + MI[k][1] + MI[k][2];
            if (deg == 0) arb_one(rng); else if (odd) { arb_zero(rng); mag_one(arb_radref(rng)); }
            else { arb_set_d(rng, 0.5); mag_set_d(arb_radref(rng), 0.5); }
            acb_mul_arb(w, Dt + k, rng, PREC); acb_add(E, E, w, PREC);
        }
        arf_t ee; arf_init(ee); arf_set_d(ee, (rD + rLoh) * (1 + 1e-10)); acb_add_error_arf(E, ee); arf_clear(ee);
        LAST_RPOLY = rD; LAST_RLOH = rLoh;
        if (!LAST_DT) LAST_DT = _acb_vec_init(128);
        _acb_vec_set(LAST_DT, Dt, NM); LAST_R = (rD + rLoh) * (1 + 1e-10);
        arb_clear(rng); acb_clear(w); _acb_vec_clear(Dt, NM);
    }
    _acb_vec_clear(TL, NM); _acb_vec_clear(TB, NM); _acb_vec_clear(TE, NM); _acb_vec_clear(A, 36 * NM); _acb_vec_clear(A1, 36 * NM);
    tsmat_clear(&C0L); tsmat_clear(&B0L); acb_clear(MC);
    tsmat_clear(&C0); tsmat_clear(&C1); tsmat_clear(&B0); _acb_vec_clear(Y, n20); _acb_vec_clear(Sy.W0, NW); _acb_vec_clear(Sy.b, NW);
    acb_clear(z); arb_clear(t); arb_clear(h); arb_clear(sk);
    return ok;
}

/* range of the last D-polynomial over the sub-box prod [a_i,b_i] (a_i,b_i in [-1,1]) plus the remainder LAST_R */
static int poly_excl(const double *a, const double *b)
{
    arb_t T[3], pw, m; acb_t E, w; for (int i = 0; i < 3; i++) arb_init(T[i]); arb_init(pw); arb_init(m); acb_init(E); acb_init(w);
    for (int i = 0; i < 3; i++) { arb_set_d(T[i], 0.5 * (a[i] + b[i])); mag_set_d(arb_radref(T[i]), 0.5 * (b[i] - a[i]) * (1 + 1e-12) + 1e-300); }
    acb_zero(E);
    for (int k = 0; k < NM; k++) {
        if (acb_is_zero(LAST_DT + k)) continue;
        arb_one(m);
        for (int i = 0; i < 3; i++) if (MI[k][i]) { arb_pow_ui(pw, T[i], MI[k][i], PREC); arb_mul(m, m, pw, PREC); }
        acb_mul_arb(w, LAST_DT + k, m, PREC); acb_add(E, E, w, PREC);
    }
    arf_t e; arf_init(e); arf_set_d(e, LAST_R); acb_add_error_arf(E, e); arf_clear(e);
    int ok = acb_is_finite(E) && !acb_contains_zero(E);
    for (int i = 0; i < 3; i++) arb_clear(T[i]); arb_clear(pw); arb_clear(m); acb_clear(E); acb_clear(w);
    return ok;
}
/* exclusion on the t1-piece [a1,b1] for ALL (t2,t3) in [-1,1]^2 by bisection of the parameter square */
static int poly_excl_params(double a1, double b1, double a2, double b2, double a3, double b3, int depth)
{
    double a[3] = {a1, a2, a3}, b[3] = {b1, b2, b3};
    if (poly_excl(a, b)) return 1;
    if (depth <= 0) return 0;
    if (b2 - a2 >= b3 - a3) { double m = 0.5 * (a2 + b2); return poly_excl_params(a1, b1, a2, m, a3, b3, depth - 1) && poly_excl_params(a1, b1, m, b2, a3, b3, depth - 1); }
    double m = 0.5 * (a3 + b3); return poly_excl_params(a1, b1, a2, b2, a3, m, depth - 1) && poly_excl_params(a1, b1, a2, b2, m, b3, depth - 1);
}
