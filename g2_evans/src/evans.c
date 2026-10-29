/* evans.c -- validated Evans function (FLINT/Arb, C) for the micropolar Couette problem, case (a).
 *
 * Port of arb_evans.ArbEvansC (python-flint, G2 v2 report, sec. 3.1-3.2) plus:
 *   - derivative chain V_k = (1/k!) d^k W / dLam^k (k = 0..K+1), propagated with the same Taylor/majorant integrator;
 *   - rigorous segment tests for the argument principle:
 *       (B)  direct ball enclosure of D over the segment (= ball_winding of the Python code);
 *       (T)  Taylor model in Lambda:  D(c+e) = sum_{k<=K} D_k(c) e^k + R,  |R| <= sup_ball |D_{K+1}| r^{K+1}
 *            (D_k at the point c, D_{K+1} enclosed over the ball);  this removes most of the "wrapping" in Lambda;
 *   - parameter boxes (sigma, alpha) as Arb balls (the count is constant on the box when D != 0 on contour x box).
 *
 * 6th-order system (x = s/S in [0,1], Y = (zeta, zeta', w, w', Phi, Phi')):
 *   zeta'' = q1 zeta + k1 w,  w'' = k2 zeta + qg w,  Phi'' = zeta + c w + beta^2 Phi   (derivatives in s)
 *   q_nu(s) = (Lam + nu beta^2 + s_nu sigma^{-2/3} - i s)/nu, k1 = k12 sigma^{-2/3}/nu1, k2 = k21 sigma^{-2/3}/g
 *   BC at s=0 and s=S: zeta-(2-c)w = Phi = Phi' = 0.   D(Lam) = b . W(1),  W = y1^y2^y3 in Lambda^3 C^6 (dim 20),
 *   W' = A^{(3)} W (additive compound), b = 3x3 minors of B_S (Cauchy-Binet).
 * Taylor step (A = B0 + A1 t on [x_k, x_k+h]): (m+1) Y_{m+1} = B0 Y_m + A1 Y_{m-1} (+ Cd Y^{(k-1)}_m for the chain);
 *   tail: c_m = |Y_m|_inf h^m, theta = ((|B0|+|Cd|) h + |A1| h^2)/(n+1) < 1/2  =>  sum_{m>n} <= max(c_n,c_{n-1}) 2 theta/(1-theta).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arf.h>

#define NW 20
static slong PREC = 128;
static int DEG = 50;
static double HSCALE = 8.0;

typedef struct { arb_t g, kap, nu1, c, s1, sg, k12, k21, delta_p; } params_t;

static void params_init(params_t *P, const char *Ns, const char *Js, const char *Gs)
{
    arb_t N, J, G, t, u;
    arb_init(N); arb_init(J); arb_init(G); arb_init(t); arb_init(u);
    arb_set_str(N, Ns, PREC); arb_set_str(J, Js, PREC); arb_set_str(G, Gs, PREC);
    arb_init(P->g); arb_init(P->kap); arb_init(P->nu1); arb_init(P->c); arb_init(P->s1); arb_init(P->sg);
    arb_init(P->k12); arb_init(P->k21); arb_init(P->delta_p);
    arb_div(P->g, G, J, PREC);
    arb_mul_ui(P->kap, N, 4, PREC); arb_div(P->kap, P->kap, J, PREC);
    arb_add_ui(P->nu1, N, 1, PREC);
    arb_sub(t, P->nu1, P->g, PREC);                       /* t = 1+N-g */
    arb_mul_ui(P->c, N, 2, PREC); arb_div(P->c, P->c, t, PREC);
    arb_mul(P->s1, P->kap, N, PREC); arb_div(P->s1, P->s1, t, PREC);
    arb_sub_ui(u, P->g, 1, PREC); arb_neg(u, u);          /* u = 1-g */
    arb_mul(P->sg, P->kap, u, PREC); arb_div(P->sg, P->sg, t, PREC);
    arb_mul(P->k12, N, P->kap, PREC); arb_mul(P->k12, P->k12, u, PREC); arb_mul_si(P->k12, P->k12, -2, PREC);
    arb_div(P->k12, P->k12, t, PREC); arb_div(P->k12, P->k12, t, PREC);
    arb_mul_2exp_si(P->k21, P->kap, -1); arb_neg(P->k21, P->k21);
    /* delta' = 0.05 g^{1/3} */
    arb_set_str(t, "0.05", PREC); arb_root_ui(u, P->g, 3, PREC); arb_mul(P->delta_p, t, u, PREC);
    arb_clear(N); arb_clear(J); arb_clear(G); arb_clear(t); arb_clear(u);
}

/* ---------------- sparse 20x20 compound matrices ---------------- */
typedef struct { int n; int row[400], col[400]; acb_ptr val; } smat;
static void smat_init(smat *M) { M->n = 0; M->val = _acb_vec_init(400); }
static void smat_clear(smat *M) { _acb_vec_clear(M->val, 400); }
static int TRIP[NW][3];
static int tidx(int a, int b, int c) { for (int i = 0; i < NW; i++) if (TRIP[i][0] == a && TRIP[i][1] == b && TRIP[i][2] == c) return i; return -1; }
static void trip_init(void) { int n = 0; for (int i = 0; i < 6; i++) for (int j = i + 1; j < 6; j++) for (int k = j + 1; k < 6; k++) { TRIP[n][0] = i; TRIP[n][1] = j; TRIP[n][2] = k; n++; } }
/* additive compound of a 6x6 (dense acb array A[36]) */
static void compound3(smat *C, acb_srcptr A)
{
    acb_ptr D = _acb_vec_init(NW * NW);
    for (int J = 0; J < NW; J++)
        for (int pos = 0; pos < 3; pos++)
            for (int m = 0; m < 6; m++) {
                const acb_struct *a = A + m * 6 + TRIP[J][pos];
                if (acb_is_zero(a)) continue;
                int t[3] = {TRIP[J][0], TRIP[J][1], TRIP[J][2]}; t[pos] = m;
                if (t[0] == t[1] || t[0] == t[2] || t[1] == t[2]) continue;
                int s = 1;   /* bubble sort with sign */
                for (int i = 0; i < 3; i++) for (int j = 0; j < 2 - i; j++) if (t[j] > t[j + 1]) { int x = t[j]; t[j] = t[j + 1]; t[j + 1] = x; s = -s; }
                int I = tidx(t[0], t[1], t[2]);
                if (s > 0) acb_add(D + I * NW + J, D + I * NW + J, a, PREC); else acb_sub(D + I * NW + J, D + I * NW + J, a, PREC);
            }
    C->n = 0;
    for (int I = 0; I < NW; I++) for (int J = 0; J < NW; J++) if (!acb_is_zero(D + I * NW + J)) {
        C->row[C->n] = I; C->col[C->n] = J; acb_set(C->val + C->n, D + I * NW + J); C->n++; }
    _acb_vec_clear(D, NW * NW);
}
static void smat_mulvec_add(acb_ptr out, const smat *M, acb_srcptr v)
{
    acb_t t; acb_init(t);
    for (int e = 0; e < M->n; e++) { acb_mul(t, M->val + e, v + M->col[e], PREC); acb_add(out + M->row[e], out + M->row[e], t, PREC); }
    acb_clear(t);
}
static double smat_norm(const smat *M)   /* induced inf-norm, upper bound */
{
    double rs[NW] = {0}; arf_t u; arf_init(u);
    for (int e = 0; e < M->n; e++) { acb_get_abs_ubound_arf(u, M->val + e, 64); rs[M->row[e]] += arf_get_d(u, ARF_RND_UP); }
    arf_clear(u); double m = 0; for (int i = 0; i < NW; i++) if (rs[i] > m) m = rs[i];
    return m * (1 + 1e-12);
}
static double vec_maxabs(acb_srcptr v, int n)
{
    arf_t u; arf_init(u); double m = 0;
    for (int i = 0; i < n; i++) { acb_get_abs_ubound_arf(u, v + i, 64); double d = arf_get_d(u, ARF_RND_UP); if (d > m) m = d; }
    arf_clear(u); return m;
}

/* ---------------- the system ---------------- */
typedef struct { smat C0, C1, Cd; acb_ptr W0, b; } system_t;
static void build_system(system_t *Sy, const params_t *P, const acb_t Lam, const arb_t sigma, const arb_t alpha)
{
    arb_t S, beta, s23, t, third; acb_t L1, Lg, k1, k2, z;
    arb_init(S); arb_init(beta); arb_init(s23); arb_init(t); arb_init(third);
    acb_init(L1); acb_init(Lg); acb_init(k1); acb_init(k2); acb_init(z);
    arb_root_ui(t, sigma, 3, PREC);                                   /* sigma^{1/3} */
    arb_mul_2exp_si(S, t, 1); arb_div(beta, alpha, t, PREC); arb_mul(s23, t, t, PREC); arb_inv(s23, s23, PREC);
    acb_ptr A0 = _acb_vec_init(36), A1 = _acb_vec_init(36), dA = _acb_vec_init(36);
    /* L1 = (Lam + nu1 beta^2 + s1 s23)/nu1 */
    arb_mul(t, beta, beta, PREC); arb_mul(t, t, P->nu1, PREC); acb_add_arb(L1, Lam, t, PREC);
    arb_mul(t, P->s1, s23, PREC); acb_add_arb(L1, L1, t, PREC); acb_div_arb(L1, L1, P->nu1, PREC);
    arb_mul(t, beta, beta, PREC); arb_mul(t, t, P->g, PREC); acb_add_arb(Lg, Lam, t, PREC);
    arb_mul(t, P->sg, s23, PREC); acb_add_arb(Lg, Lg, t, PREC); acb_div_arb(Lg, Lg, P->g, PREC);
    arb_mul(t, P->k12, s23, PREC); arb_div(t, t, P->nu1, PREC); acb_set_arb(k1, t);
    arb_mul(t, P->k21, s23, PREC); arb_div(t, t, P->g, PREC); acb_set_arb(k2, t);
#define E(M, i, j) ((M) + (i) * 6 + (j))
    acb_set_arb(E(A0, 0, 1), S); acb_mul_arb(E(A0, 1, 0), L1, S, PREC); acb_mul_arb(E(A0, 1, 2), k1, S, PREC);
    acb_set_arb(E(A0, 2, 3), S); acb_mul_arb(E(A0, 3, 0), k2, S, PREC); acb_mul_arb(E(A0, 3, 2), Lg, S, PREC);
    acb_set_arb(E(A0, 4, 5), S); acb_set_arb(E(A0, 5, 0), S); arb_mul(t, S, P->c, PREC); acb_set_arb(E(A0, 5, 2), t);
    arb_mul(t, beta, beta, PREC); arb_mul(t, t, S, PREC); acb_set_arb(E(A0, 5, 4), t);
    if (acb_is_zero(E(A0, 5, 4))) { /* keep structural entry when beta == 0 exactly: harmless */ }
    arb_mul(t, S, S, PREC); arb_div(t, t, P->nu1, PREC); arb_neg(t, t); acb_set_arb(z, t); acb_mul_onei(E(A1, 1, 0), z);
    arb_mul(t, S, S, PREC); arb_div(t, t, P->g, PREC); arb_neg(t, t); acb_set_arb(z, t); acb_mul_onei(E(A1, 3, 2), z);
    arb_div(t, S, P->nu1, PREC); acb_set_arb(E(dA, 1, 0), t); arb_div(t, S, P->g, PREC); acb_set_arb(E(dA, 3, 2), t);
#undef E
    compound3(&Sy->C0, A0); compound3(&Sy->C1, A1); compound3(&Sy->Cd, dA);
    _acb_vec_clear(A0, 36); _acb_vec_clear(A1, 36); _acb_vec_clear(dA, 36);
    arb_clear(S); arb_clear(beta); arb_clear(s23); arb_clear(t); arb_clear(third);
    acb_clear(L1); acb_clear(Lg); acb_clear(k1); acb_clear(k2); acb_clear(z);
}
/* W0 = Pluecker of N0 (6x3), b = minors of B_S (3x6) */
static void det3(acb_t d, acb_srcptr M)   /* M row-major 3x3 */
{
    acb_t a, b; acb_init(a); acb_init(b); acb_zero(d);
    int p[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{1,0,2},{2,1,0}}; int sg[6] = {1,1,1,-1,-1,-1};
    for (int k = 0; k < 6; k++) { acb_mul(a, M + 0 * 3 + p[k][0], M + 1 * 3 + p[k][1], PREC); acb_mul(a, a, M + 2 * 3 + p[k][2], PREC);
        if (sg[k] > 0) acb_add(d, d, a, PREC); else acb_sub(d, d, a, PREC); }
    acb_clear(a); acb_clear(b);
}
static void boundary_init(system_t *Sy, const params_t *P)
{
    acb_ptr N0 = _acb_vec_init(18), BS = _acb_vec_init(18), M = _acb_vec_init(9);
    arb_t t; arb_init(t);
    arb_sub_ui(t, P->c, 2, PREC); arb_neg(t, t);                       /* 2-c */
    acb_set_arb(N0 + 0 * 3 + 0, t); acb_one(N0 + 2 * 3 + 0); acb_one(N0 + 1 * 3 + 1); acb_one(N0 + 3 * 3 + 2);
    acb_one(BS + 0 * 6 + 0); arb_neg(t, t); acb_set_arb(BS + 0 * 6 + 2, t); acb_one(BS + 1 * 6 + 4); acb_one(BS + 2 * 6 + 5);
    Sy->W0 = _acb_vec_init(NW); Sy->b = _acb_vec_init(NW);
    for (int I = 0; I < NW; I++) {
        for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) acb_set(M + r * 3 + c, N0 + TRIP[I][r] * 3 + c);
        det3(Sy->W0 + I, M);
        for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) acb_set(M + r * 3 + c, BS + r * 6 + TRIP[I][c]);
        det3(Sy->b + I, M);
    }
    _acb_vec_clear(N0, 18); _acb_vec_clear(BS, 18); _acb_vec_clear(M, 9); arb_clear(t);
}

/* ---------------- Taylor propagation of the chain V_0..V_L (L = levels-1) ---------------- */
static long NEVAL_VEC = 0;
static int propagate(acb_ptr V /* levels*NW, in/out */, int levels, const system_t *Sy, int nsteps)
{
    double M1 = smat_norm(&Sy->C1), Md = levels > 1 ? smat_norm(&Sy->Cd) : 0.0;
    double h = 1.0 / nsteps;
    acb_ptr Ycur = _acb_vec_init(levels * NW), Yprev = _acb_vec_init(levels * NW), Ynew = _acb_vec_init(levels * NW), acc = _acb_vec_init(levels * NW);
    smat B0; smat_init(&B0);
    arb_t hp, xk, inv; arb_init(hp); arb_init(xk); arb_init(inv); acb_t t; acb_init(t);
    arf_t err; arf_init(err);
    int ok = 1;
    for (int k = 0; k < nsteps; k++) {
        /* B0 = C0 + C1 * x_k (same sparsity as C0 union C1) */
        arb_set_si(xk, k); arb_div_si(xk, xk, nsteps, PREC);
        B0.n = 0;
        for (int e = 0; e < Sy->C0.n; e++) { B0.row[B0.n] = Sy->C0.row[e]; B0.col[B0.n] = Sy->C0.col[e]; acb_set(B0.val + B0.n, Sy->C0.val + e); B0.n++; }
        for (int e = 0; e < Sy->C1.n; e++) {
            int f = -1; for (int q = 0; q < B0.n; q++) if (B0.row[q] == Sy->C1.row[e] && B0.col[q] == Sy->C1.col[e]) { f = q; break; }
            acb_mul_arb(t, Sy->C1.val + e, xk, PREC);
            if (f < 0) { B0.row[B0.n] = Sy->C1.row[e]; B0.col[B0.n] = Sy->C1.col[e]; acb_set(B0.val + B0.n, t); B0.n++; }
            else acb_add(B0.val + f, B0.val + f, t, PREC);
        }
        double M0 = smat_norm(&B0);
        double th = ((M0 + Md) * h + M1 * h * h) / (DEG + 1);
        if (th >= 0.5) { ok = 0; break; }
        _acb_vec_set(Ycur, V, levels * NW); _acb_vec_zero(Yprev, levels * NW); _acb_vec_set(acc, V, levels * NW);
        arb_one(hp);
        arb_set_si(inv, 1); arb_div_si(inv, inv, nsteps, PREC);       /* exact h as arb */
        for (int m = 0; m < DEG; m++) {
            _acb_vec_zero(Ynew, levels * NW);
            for (int l = 0; l < levels; l++) {
                smat_mulvec_add(Ynew + l * NW, &B0, Ycur + l * NW);
                if (m > 0) smat_mulvec_add(Ynew + l * NW, &Sy->C1, Yprev + l * NW);
                if (l > 0) smat_mulvec_add(Ynew + l * NW, &Sy->Cd, Ycur + (l - 1) * NW);
            }
            _acb_vec_scalar_div_ui(Ynew, Ynew, levels * NW, m + 1, PREC);
            arb_mul(hp, hp, inv, PREC);
            for (int i = 0; i < levels * NW; i++) { acb_mul_arb(t, Ynew + i, hp, PREC); acb_add(acc + i, acc + i, t, PREC); }
            acb_ptr sw = Yprev; Yprev = Ycur; Ycur = Ynew; Ynew = sw;
        }
        /* tail: Ycur = Y_n, Yprev = Y_{n-1} */
        double cN = vec_maxabs(Ycur, levels * NW) * pow(h, DEG), cN1 = vec_maxabs(Yprev, levels * NW) * pow(h, DEG - 1);
        double tail = fmax(cN, cN1) * 2 * th / (1 - th) * (1 + 1e-10);
        arf_set_d(err, tail);
        for (int i = 0; i < levels * NW; i++) acb_add_error_arf(acc + i, err);
        _acb_vec_set(V, acc, levels * NW);
        NEVAL_VEC += levels;
    }
    _acb_vec_clear(Ycur, levels * NW); _acb_vec_clear(Yprev, levels * NW); _acb_vec_clear(Ynew, levels * NW); _acb_vec_clear(acc, levels * NW);
    smat_clear(&B0); arb_clear(hp); arb_clear(xk); arb_clear(inv); acb_clear(t); arf_clear(err);
    return ok;
}

/* D_k (k=0..levels-1), k-th Taylor coefficient in Lambda: D_k = b . V_k ;  Lam, sigma, alpha may be balls */
static int evans_chain(acb_ptr Dk, int levels, const params_t *P, const acb_t Lam, const arb_t sigma, const arb_t alpha)
{
    system_t Sy; smat_init(&Sy.C0); smat_init(&Sy.C1); smat_init(&Sy.Cd);
    build_system(&Sy, P, Lam, sigma, alpha); boundary_init(&Sy, P);
    double Mn = smat_norm(&Sy.C0) + smat_norm(&Sy.C1) + (levels > 1 ? smat_norm(&Sy.Cd) : 0.0);
    int nsteps = (int) ceil(Mn / HSCALE); if (nsteps < 4) nsteps = 4;
    acb_ptr V = _acb_vec_init(levels * NW);
    _acb_vec_set(V, Sy.W0, NW);
    int ok = propagate(V, levels, &Sy, nsteps);
    for (int l = 0; l < levels; l++) acb_dot(Dk + l, NULL, 0, Sy.b, 1, V + l * NW, 1, NW, PREC);
    _acb_vec_clear(V, levels * NW); _acb_vec_clear(Sy.W0, NW); _acb_vec_clear(Sy.b, NW);
    smat_clear(&Sy.C0); smat_clear(&Sy.C1); smat_clear(&Sy.Cd);
    return ok;
}

static double now(void);
#include "tm.inc"

/* ---------------- segment tests and winding ---------------- */
typedef struct { const params_t *P; arb_t sig, alp;      /* parameter box (balls) */
                 arb_t sigc, alpc;                         /* centre (points) for the endpoint values */
                 int K; double h0; int maxdepth; long nseg, nball, ntay, nfail;
                 int use_tm; double b2lo, b2hi, Slo, Shi; long ntm, npieces; int open; } wctx;

/* enclosure of D over the axis-parallel segment [a,b]; returns 1 if 0 not in enclosure */
static int seg_ok(wctx *W, double ax, double ay, double bx, double by, int *how)
{
    acb_t L, E; acb_init(L); acb_init(E);
    double cx = 0.5 * (ax + bx), cy = 0.5 * (ay + by), hx = fabs(bx - ax) / 2, hy = fabs(by - ay) / 2;
    arb_set_d(acb_realref(L), cx); arb_set_d(acb_imagref(L), cy);
    mag_set_d(arb_radref(acb_realref(L)), hx * (1 + 1e-12) + 1e-300); mag_set_d(arb_radref(acb_imagref(L)), hy * (1 + 1e-12) + 1e-300);
    int ok = 0;
    if (W->use_tm) {      /* Taylor model in (Lam along the segment, beta^2, S) */
        acb_t Lc, u; acb_init(Lc); acb_init(u);
        arb_set_d(acb_realref(Lc), cx); arb_set_d(acb_imagref(Lc), cy);
        arb_set_d(acb_realref(u), 0.5 * (bx - ax)); arb_set_d(acb_imagref(u), 0.5 * (by - ay));
        /* outward: the segment endpoints are exactly representable; u, Lc computed exactly in double for dyadic inputs is not
           guaranteed -> add tiny radius */
        mag_set_d(arb_radref(acb_realref(Lc)), 1e-15 * (1 + fabs(cx))); mag_set_d(arb_radref(acb_imagref(Lc)), 1e-15 * (1 + fabs(cy)));
        double tt0 = now();
        int r = evans_tm(E, W->P, Lc, u, W->b2lo, W->b2hi, W->Slo, W->Shi); W->ntm++;
        int npieces = 0;
        if (r) for (int m1 = 1; m1 <= 32 && !npieces; m1 *= 2) {
            int all = 1;
            for (int q = 0; q < m1 && all; q++) all = poly_excl_params(-1 + 2.0 * q / m1, -1 + 2.0 * (q + 1) / m1, -1, 1, -1, 1, 8);
            if (all) npieces = m1;
        }

        if (npieces) { ok = 1; *how = npieces; }
        if (getenv("EVANS_VERBOSE")) { fprintf(stderr, "seg [%.4f%+.4fi, %.4f%+.4fi] pieces=%d tm ok=%d zero_in=%d |E|~%.3e rad~%.3e (poly %.2e, lohner %.2e) %.2fs\n", ax, ay, bx, by, npieces, r,
            r ? acb_contains_zero(E) : -1, r ? acb_ub(E) : 0.0,
            r ? mag_get_d(arb_radref(acb_realref(E))) : 0.0, LAST_RPOLY, LAST_RLOH, now() - tt0); }
        acb_clear(Lc); acb_clear(u); acb_clear(L); acb_clear(E);
        return ok;
    }
    /* (B) direct ball */
    evans_chain(E, 1, W->P, L, W->sig, W->alp); W->nball++;
    if (acb_is_finite(E) && !acb_contains_zero(E)) { ok = 1; *how = 0; }
    if (!ok && W->K >= 0) {
        /* (T) Taylor model of order K in Lambda */
        int K = W->K;
        acb_ptr Dc = _acb_vec_init(K + 1), Db = _acb_vec_init(K + 2);
        acb_t c, eps, P; acb_init(c); acb_init(eps); acb_init(P);
        arb_set_d(acb_realref(c), cx); arb_set_d(acb_imagref(c), cy);
        evans_chain(Dc, K + 1, W->P, c, W->sig, W->alp);
        evans_chain(Db, K + 2, W->P, L, W->sig, W->alp); W->ntay++;
        /* eps ball: |Re| <= hx, |Im| <= hy ; r = sqrt(hx^2+hy^2) */
        acb_zero(eps); mag_set_d(arb_radref(acb_realref(eps)), hx * (1 + 1e-12)); mag_set_d(arb_radref(acb_imagref(eps)), hy * (1 + 1e-12));
        acb_zero(P);
        for (int k = K; k >= 0; k--) { acb_mul(P, P, eps, PREC); acb_add(P, P, Dc + k, PREC); }    /* Horner */
        arf_t u; arf_init(u); acb_get_abs_ubound_arf(u, Db + K + 1, 64);
        double R = arf_get_d(u, ARF_RND_UP) * pow(sqrt(hx * hx + hy * hy) * (1 + 1e-12), K + 1);
        arf_set_d(u, R * (1 + 1e-10)); acb_add_error_arf(P, u); arf_clear(u);
        if (acb_is_finite(P) && !acb_contains_zero(P)) { ok = 1; *how = 1; }
        _acb_vec_clear(Dc, K + 1); _acb_vec_clear(Db, K + 2); acb_clear(c); acb_clear(eps); acb_clear(P);
    }
    acb_clear(L); acb_clear(E);
    return ok;
}
/* endpoint value at the centre parameters (point), cached for the previous endpoint */
static void Dpoint(acb_t d, wctx *W, double x, double y)
{
    acb_t L; acb_init(L); arb_set_d(acb_realref(L), x); arb_set_d(acb_imagref(L), y);
    evans_chain(d, 1, W->P, L, W->sigc, W->alpc); acb_clear(L);
}
/* winding of D along the closed polygon (axis-parallel sides); returns 1 on success, *wind enclosure in turns */
static int winding(wctx *W, double *vx, double *vy, int nv, arb_t wind, double *minrat)
{
    arb_t tot, ar; acb_t Da, Db, q; arb_init(tot); arb_init(ar); acb_init(Da); acb_init(Db); acb_init(q);
    arb_zero(tot);
    typedef struct { double ax, ay, bx, by; int depth; } seg;
    int cap = 1 << 16; seg *st = malloc(sizeof(seg) * cap);
    int okall = 1;
    for (int s = 0; s < (W->open ? nv - 1 : nv) && okall; s++) {
        double ax = vx[s], ay = vy[s], bx = vx[(s + 1) % nv], by = vy[(s + 1) % nv];
        double len = hypot(bx - ax, by - ay); int n = (int) ceil(len / (2 * W->h0)); if (n < 1) n = 1;
        /* push pieces in reverse so that they are processed in order */
        int sp = 0;
        for (int k = n - 1; k >= 0; k--) { st[sp].ax = ax + (bx - ax) * k / n; st[sp].ay = ay + (by - ay) * k / n;
            st[sp].bx = (k + 1 == n) ? bx : ax + (bx - ax) * (k + 1) / n; st[sp].by = (k + 1 == n) ? by : ay + (by - ay) * (k + 1) / n;
            st[sp].depth = 0; sp++; }     /* exact closure: last piece ends exactly at the next vertex */
        while (sp > 0 && okall) {
            seg g = st[--sp]; int how;
            if (seg_ok(W, g.ax, g.ay, g.bx, g.by, &how)) {
                W->nseg++;
                int np = (W->use_tm && how >= 1) ? how : 1;
                for (int q2 = 0; q2 < np; q2++) {
                    double xa = g.ax + (g.bx - g.ax) * q2 / np, ya = g.ay + (g.by - g.ay) * q2 / np;
                    double xb = (q2 + 1 == np) ? g.bx : g.ax + (g.bx - g.ax) * (q2 + 1) / np, yb = (q2 + 1 == np) ? g.by : g.ay + (g.by - g.ay) * (q2 + 1) / np;
                    /* the piece [xa,xb] lies inside the t1-piece certified by poly_excl (both computed from the same
                       double expressions up to rounding <= 1e-15 relative; the TM centre Lc carries a 1e-15 radius) */
                    Dpoint(Da, W, xa, ya); Dpoint(Db, W, xb, yb);
                    acb_div(q, Db, Da, PREC); acb_arg(ar, q, PREC); arb_add(tot, tot, ar, PREC);
                }
                W->npieces += np;
            } else {
                if (g.depth >= W->maxdepth || sp + 2 >= cap) { okall = 0; W->nfail++; break; }
                double mx = 0.5 * (g.ax + g.bx), my = 0.5 * (g.ay + g.by);
                st[sp].ax = mx; st[sp].ay = my; st[sp].bx = g.bx; st[sp].by = g.by; st[sp].depth = g.depth + 1; sp++;
                st[sp].ax = g.ax; st[sp].ay = g.ay; st[sp].bx = mx; st[sp].by = my; st[sp].depth = g.depth + 1; sp++;
            }
        }
    }
    free(st);
    arb_const_pi(ar, PREC); arb_mul_2exp_si(ar, ar, 1); arb_div(wind, tot, ar, PREC);
    arb_clear(tot); arb_clear(ar); acb_clear(Da); acb_clear(Db); acb_clear(q);
    (void) minrat;
    return okall;
}

/* Lambda* table of Theorem 2.1 (G2 v2), (1;1;0.5) */
static double lam_star(double sig_hi_unused, double sig_lo)
{
    if (sig_lo < 5) return 15;
    if (sig_lo < 8) return 12;
    if (sig_lo < 12) return 10;
    if (sig_lo < 20) return 8;
    if (sig_lo < 40) return 6;
    if (sig_lo < 80) return 5;
    return 4;
}

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static void usage(void)
{
    fprintf(stderr,
        "usage:\n"
        "  evans point N J G alpha sigma ReLam ImLam [levels]   -> D_k, k<levels (enclosures)\n"
        "  evans ball  N J G alpha sigma ReLam ImLam r          -> D over the square ball of radius r\n"
        "  evans box   N J G alpha_lo alpha_hi sigma_lo sigma_hi [K h0 maxdepth Lstar]\n"
        "        -> rigorous zero count of D in B(Lam*) = [-delta',Lam*] x [-Lam*, S_hi+Lam*] for ALL (alpha,sigma) in the box\n"
        "           (one JSON line on stdout).  K = Taylor order in Lambda (-1: balls only), default 6.\n");
    exit(1);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage();
    trip_init();
    if (getenv("EVANS_PREC")) PREC = atol(getenv("EVANS_PREC"));
    if (getenv("EVANS_DEG")) DEG = atoi(getenv("EVANS_DEG"));
    if (getenv("EVANS_HSCALE")) HSCALE = atof(getenv("EVANS_HSCALE"));
    /* defaults by regime (measured, report v2): sigma < 15 -> plain Taylor models of degree 3 (cheaper, Lohner does not pay);
       sigma >= 15 -> Lohner remainder + anisotropic degrees (5 total, 3 in the parameters).  Env variables override. */
    double sig_hint = (argc > 7 && (!strcmp(argv[1], "box") || !strcmp(argv[1], "tmprobe"))) ? atof(argv[7]) : 100.0;
    int lowsig = sig_hint < 15.0;
    tm_setup2(getenv("EVANS_TK") ? atoi(getenv("EVANS_TK")) : (lowsig ? 3 : 5), getenv("EVANS_TKP") ? atoi(getenv("EVANS_TKP")) : (lowsig ? 99 : 3));
    LOHNER = lowsig ? 0 : 1;
    if (getenv("EVANS_LOHNER")) LOHNER = atoi(getenv("EVANS_LOHNER"));
    if (getenv("EVANS_MSHIFT")) MSHIFT = atof(getenv("EVANS_MSHIFT"));
    params_t P;
    if (!strcmp(argv[1], "point") || !strcmp(argv[1], "ball")) {
        if (argc < 9) usage();
        params_init(&P, argv[2], argv[3], argv[4]);
        arb_t al, sg; acb_t L; arb_init(al); arb_init(sg); acb_init(L);
        /* inputs read as IEEE doubles (exactly the binary numbers python-flint receives from Python floats) */
        arb_set_d(al, atof(argv[5])); arb_set_d(sg, atof(argv[6]));
        arb_set_d(acb_realref(L), atof(argv[7])); arb_set_d(acb_imagref(L), atof(argv[8]));
        int levels = 1;
        if (!strcmp(argv[1], "point")) { if (argc > 9) levels = atoi(argv[9]); }
        else { double r = atof(argv[9]); mag_set_d(arb_radref(acb_realref(L)), r); mag_set_d(arb_radref(acb_imagref(L)), r); }
        acb_ptr D = _acb_vec_init(levels);
        double t0 = now(); int ok = evans_chain(D, levels, &P, L, sg, al); double t1 = now();
        for (int k = 0; k < levels; k++) { printf("D_%d = ", k); acb_printn(D + k, 20, ARB_STR_NO_RADIUS);
            printf("   rad(re,im) = %.3e %.3e\n", mag_get_d(arb_radref(acb_realref(D + k))), mag_get_d(arb_radref(acb_imagref(D + k)))); }
        printf("ok=%d time=%.4f s\n", ok, t1 - t0);
        return 0;
    }
    if (!strcmp(argv[1], "tmprobe")) {   /* tmprobe N J G a1 a2 s1 s2 ReL ImL dRe dIm : one TM evaluation on the segment [L, L+d] */
        params_init(&P, argv[2], argv[3], argv[4]);
        double a1 = atof(argv[5]), a2 = atof(argv[6]), s1 = atof(argv[7]), s2 = atof(argv[8]);
        double x = atof(argv[9]), y = atof(argv[10]), dx = atof(argv[11]), dy = atof(argv[12]);
        double b2lo = a1 * a1 / pow(s2, 2.0 / 3) * (1 - 1e-12), b2hi = a2 * a2 / pow(s1, 2.0 / 3) * (1 + 1e-12);
        double Slo = 2 * cbrt(s1) * (1 - 1e-12), Shi = 2 * cbrt(s2) * (1 + 1e-12);
        acb_t Lc, u, E; acb_init(Lc); acb_init(u); acb_init(E);
        arb_set_d(acb_realref(Lc), x + dx / 2); arb_set_d(acb_imagref(Lc), y + dy / 2);
        arb_set_d(acb_realref(u), dx / 2); arb_set_d(acb_imagref(u), dy / 2);
        double t0 = now(); int r = evans_tm(E, &P, Lc, u, b2lo, b2hi, Slo, Shi); double t1 = now();
        printf("c0 mid=%.3e rad=%.3e ; ", r ? hypot(arf_get_d(arb_midref(acb_realref(LAST_DT)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(LAST_DT)), ARF_RND_NEAR)) : 0.0, r ? mag_get_d(arb_radref(acb_realref(LAST_DT))) : 0.0);
        printf("ok=%d |c0|=%.4e  R_poly=%.3e R_lohner=%.3e  ratio=%.3e  t=%.2fs\n", r, r ? acb_ub(LAST_DT) : 0.0, LAST_RPOLY, LAST_RLOH,
               r ? LAST_R / acb_ub(LAST_DT) : 0.0, t1 - t0);
        return 0;
    }
    if (!strcmp(argv[1], "box")) {
        if (argc < 9) usage();
        params_init(&P, argv[2], argv[3], argv[4]);
        double a1 = atof(argv[5]), a2 = atof(argv[6]), s1 = atof(argv[7]), s2 = atof(argv[8]);
        wctx W; W.P = &P; W.K = argc > 9 ? atoi(argv[9]) : 6; W.h0 = argc > 10 ? atof(argv[10]) : 0.25;
        W.maxdepth = argc > 11 ? atoi(argv[11]) : 14; W.nseg = W.nball = W.ntay = W.nfail = 0;
        double Ls = argc > 12 ? atof(argv[12]) : lam_star(s2, s1);
        arb_init(W.sig); arb_init(W.alp); arb_init(W.sigc); arb_init(W.alpc);
        /* boxes as balls (outward) */
        arb_set_str(W.sig, argv[7], PREC); { arb_t t; arb_init(t); arb_set_str(t, argv[8], PREC); arb_union(W.sig, W.sig, t, PREC);
            arb_set_str(W.alp, argv[5], PREC); arb_set_str(t, argv[6], PREC); arb_union(W.alp, W.alp, t, PREC); arb_clear(t); }
        arb_get_mid_arb(W.sigc, W.sig); arb_get_mid_arb(W.alpc, W.alp);
        W.use_tm = getenv("EVANS_TM") ? atoi(getenv("EVANS_TM")) : (a2 > a1 || s2 > s1); W.ntm = 0; W.npieces = 0;
        { /* (beta^2, S) rectangle covering the image of the (alpha, sigma) box, outward */
          arb_t t, q; arb_init(t); arb_init(q); arf_t u; arf_init(u);
          arb_set_d(t, s2); arb_root_ui(t, t, 3, PREC); arb_mul(t, t, t, PREC); arb_set_d(q, a1); arb_mul(q, q, q, PREC); arb_div(q, q, t, PREC);
          arb_get_lbound_arf(u, q, 64); W.b2lo = arf_get_d(u, ARF_RND_DOWN);
          arb_set_d(t, s1); arb_root_ui(t, t, 3, PREC); arb_mul(t, t, t, PREC); arb_set_d(q, a2); arb_mul(q, q, q, PREC); arb_div(q, q, t, PREC);
          arb_get_ubound_arf(u, q, 64); W.b2hi = arf_get_d(u, ARF_RND_UP);
          arb_set_d(t, s1); arb_root_ui(t, t, 3, PREC); arb_mul_2exp_si(t, t, 1); arb_get_lbound_arf(u, t, 64); W.Slo = arf_get_d(u, ARF_RND_DOWN);
          arb_set_d(t, s2); arb_root_ui(t, t, 3, PREC); arb_mul_2exp_si(t, t, 1); arb_get_ubound_arf(u, t, 64); W.Shi = arf_get_d(u, ARF_RND_UP);
          if (W.b2lo < 0) W.b2lo = 0; arb_clear(t); arb_clear(q); arf_clear(u); }
        /* contour: delta' and S_hi (upper bounds, outward) */
        arf_t u; arf_init(u);
        arb_get_ubound_arf(u, P.delta_p, 64); double dp = arf_get_d(u, ARF_RND_UP);
        arb_t S; arb_init(S); arb_set_d(S, s2); arb_root_ui(S, S, 3, PREC); arb_mul_2exp_si(S, S, 1);
        arb_get_ubound_arf(u, S, 64); double Shi = arf_get_d(u, ARF_RND_UP);
        double vx[4] = {-dp, Ls, Ls, -dp}, vy[4] = {-Ls, -Ls, Shi + Ls, Shi + Ls};
        /* mode "line" (continuation, report v5): only the segment Re Lambda = -dl, Im Lambda in [-Ls, S_hi + Ls], with
           dl = LOWER bound of delta' (so the line lies in Re Lambda >= -delta', where the a priori region confines zeros).
           Certifies D != 0 on that segment for all parameters of the box; NO zero count. */
        const char *mode = getenv("EVANS_MODE") ? getenv("EVANS_MODE") : "contour";
        int linemode = !strcmp(mode, "line");
        arb_get_lbound_arf(u, P.delta_p, 64); double dl = arf_get_d(u, ARF_RND_DOWN);
        if (linemode) { vx[0] = -dl; vy[0] = -Ls; vx[1] = -dl; vy[1] = Shi + Ls; W.open = 1; } else W.open = 0;
        arb_t wind; arb_init(wind);
        double t0 = now(); NEVAL_VEC = 0;
        int ok = winding(&W, vx, vy, linemode ? 2 : 4, wind, NULL);
        double t1 = now();
        if (linemode) {
            printf("{\"NJG\":[\"%s\",\"%s\",\"%s\"],\"alpha\":[%.17g,%.17g],\"sigma\":[%.17g,%.17g],\"Lstar\":%g,\"mode\":\"line\","
                   "\"re_line\":%.17g,\"S_hi\":%.17g,\"ok\":%d,\"line_ok\":%d,\"segments\":%ld,\"tm_evals\":%ld,\"TK\":%d,\"TKP\":%d,"
                   "\"lohner\":%d,\"pieces\":%ld,\"seconds\":%.2f}\n", argv[2], argv[3], argv[4], a1, a2, s1, s2, Ls, -dl, Shi, ok, ok,
                   W.nseg, W.ntm, TK, TKP, LOHNER, W.npieces, t1 - t0);
            return ok ? 0 : 2;
        }
        /* certified integer count */
        long cnt = 0; int isint = 0;
        if (ok) { arf_get_d(arb_midref(wind), ARF_RND_NEAR); double m = arf_get_d(arb_midref(wind), ARF_RND_NEAR);
            cnt = lround(m); arb_t z; arb_init(z); arb_set_si(z, cnt);
            arb_t lo; arb_init(lo); arb_sub_si(lo, wind, cnt, PREC);
            isint = arb_contains(wind, z) && mag_cmp_2exp_si(arb_radref(wind), -2) < 0;
            arb_clear(z); arb_clear(lo); }
        printf("{\"NJG\":[\"%s\",\"%s\",\"%s\"],\"alpha\":[%.17g,%.17g],\"sigma\":[%.17g,%.17g],\"Lstar\":%g,\"delta_p\":%.17g,"
               "\"S_hi\":%.17g,\"K\":%d,\"ok\":%d,\"count\":%ld,\"count_certified\":%d,\"wind_mid\":%.6f,\"wind_rad\":%.3g,"
               "\"segments\":%ld,\"ball_evals\":%ld,\"taylor_evals\":%ld,\"tm_evals\":%ld,\"TK\":%d,\"TKP\":%d,\"lohner\":%d,\"pieces\":%ld,\"vec_steps\":%ld,\"seconds\":%.2f,\"mode\":\"contour\"}\n",
               argv[2], argv[3], argv[4], a1, a2, s1, s2, Ls, dp, Shi, W.K, ok, cnt, isint,
               arf_get_d(arb_midref(wind), ARF_RND_NEAR), mag_get_d(arb_radref(wind)), W.nseg, W.nball, W.ntay, W.ntm, TK, TKP, LOHNER, W.npieces, NEVAL_VEC, t1 - t0);
        return ok && isint ? 0 : 2;
    }
    usage();
    return 1;
}
