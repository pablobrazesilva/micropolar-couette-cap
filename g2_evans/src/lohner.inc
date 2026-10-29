/* lohner.inc -- Lohner (QR) control of the Taylor-model REMAINDER.
 *
 * Representation:  Y_true(x; p) = P(x; p) + Q e,   for every parameter p in the box (Lam on the segment, beta^2, S),
 *   P = polynomial part of the Taylor models (exact coefficients up to tiny ball radii),
 *   Q = point 20x20 matrix (exact doubles), e = complex vector with |e_i| <= rho_i.
 * Step x_k -> x_{k+1}:  Y_true(x_{k+1}) = Phi(p) Y_true(x_k) with Phi(p) the step propagator of the linear system.
 *   Phi(p) P  is computed by Taylor-model arithmetic: P' + local errors tau (componentwise disks, truncation + tails);
 *   Phi(p) Q e = Q' (Q'^{-1} Phi(p) Q) e,  and  Q'^{-1} Phi(p) Q  is contained entrywise in the ball matrix
 *   Q'^{-1} (Phi_ball Q), where Phi_ball encloses Phi(p) for ALL p in the box (one Taylor step in ball arithmetic with
 *   the coefficient matrix replaced by its range over the box, plus the majorant tail).
 *   Hence Y_true(x_{k+1}) = P' + Q' e',  |e'_i| <= rho'_i := sum_j |Q'^{-1} Phi_ball Q|_ij rho_j + sum_j |Q'^{-1}|_ij tau_j.
 *   Q' is ANY exactly representable invertible matrix (its inverse is enclosed rigorously by acb_mat_inv): the choice
 *   (QR of mid(Phi_ball Q) with columns sorted by rho-weighted norm, Lohner 1987) affects only the sharpness.
 * At the end D = b.Y = b.P + (b^T Q) e, so the remainder disk of D has radius sum_j |b^T Q|_j rho_j.
 * Why it helps: the componentwise propagation r' = |B| r of the plain Taylor models turns the rotation of the oscillatory
 * modes into growth (factor ~ exp(int |sqrt q|) instead of exp(int Re sqrt q)); in the Lohner frame the rotation is
 * carried by Q and only the (small) local errors are measured in absolute value. */
#include <complex.h>
#include <flint/acb_mat.h>

static int LOHNER = 1;

/* range of a TM over [-1,1]^3 as an acb ball (crude: monomials -> [-1,1] or [0,1]) */
static void tm_range(acb_t out, acb_srcptr c, double r)
{
    arb_t rng; arb_init(rng); acb_t w; acb_init(w);
    acb_zero(out);
    for (int k = 0; k < NM; k++) {
        if (acb_is_zero(c + k)) continue;
        int odd = (MI[k][0] & 1) || (MI[k][1] & 1) || (MI[k][2] & 1), deg = MI[k][0] + MI[k][1] + MI[k][2];
        if (deg == 0) arb_one(rng); else if (odd) { arb_zero(rng); mag_one(arb_radref(rng)); }
        else { arb_set_d(rng, 0.5); mag_set_d(arb_radref(rng), 0.5); }
        acb_mul_arb(w, c + k, rng, PREC); acb_add(out, out, w, PREC);
    }
    if (r > 0) { arf_t e; arf_init(e); arf_set_d(e, r * (1 + 1e-12)); acb_add_error_arf(out, e); arf_clear(e); }
    arb_clear(rng); acb_clear(w);
}
static void tsmat_to_ball(smat *Bb, const tsmat *T)
{
    Bb->n = T->n;
    for (int e = 0; e < T->n; e++) { Bb->row[e] = T->row[e]; Bb->col[e] = T->col[e]; tm_range(Bb->val + e, T->val + e * NM, T->r[e]); }
}
/* Phi (20x20, column-major in a vec of 400) for one step of length h (point) or delta in [-rS, rS] (tmstep) */
static int ball_prop(acb_ptr Phi, const smat *Bb, const smat *C1b, const arb_t h, int tmstep, double rS)
{
    int N = NW * NW;
    double hd = tmstep ? rS : arf_get_d(arb_midref(h), ARF_RND_UP) * (1 + 1e-15) + mag_get_d(arb_radref(h));
    double M0 = smat_norm(Bb), M1 = smat_norm(C1b);
    double th = (M0 * hd + M1 * hd * hd) / (DEG + 1);
    if (th >= 0.5) return 0;
    acb_ptr Ycur = _acb_vec_init(N), Yprev = _acb_vec_init(N), Ynew = _acb_vec_init(N);
    for (int j = 0; j < NW; j++) acb_one(Ycur + j * NW + j);
    _acb_vec_set(Phi, Ycur, N);
    arb_t hp; arb_init(hp); acb_t t; acb_init(t);
    for (int m = 0; m < DEG; m++) {
        _acb_vec_zero(Ynew, N);
        for (int j = 0; j < NW; j++) {
            smat_mulvec_add(Ynew + j * NW, Bb, Ycur + j * NW);
            if (m > 0) smat_mulvec_add(Ynew + j * NW, C1b, Yprev + j * NW);
        }
        _acb_vec_scalar_div_ui(Ynew, Ynew, N, m + 1, PREC);
        if (!tmstep) arb_pow_ui(hp, h, m + 1, PREC);
        else { arb_zero(hp); mag_set_d(arb_radref(hp), pow(rS, m + 1) * (1 + 1e-12)); }     /* delta^{m+1}, |delta| <= rS */
        for (int i = 0; i < N; i++) { acb_mul_arb(t, Ynew + i, hp, PREC); acb_add(Phi + i, Phi + i, t, PREC); }
        acb_ptr sw = Yprev; Yprev = Ycur; Ycur = Ynew; Ynew = sw;
    }
    double cN = vec_maxabs(Ycur, N) * pow(hd, DEG), cN1 = vec_maxabs(Yprev, N) * pow(hd, DEG - 1);
    double tail = fmax(cN, cN1) * 2 * th / (1 - th) * (1 + 1e-10);
    arf_t err; arf_init(err); arf_set_d(err, tail);
    for (int i = 0; i < N; i++) acb_add_error_arf(Phi + i, err);
    arf_clear(err); arb_clear(hp); acb_clear(t);
    _acb_vec_clear(Ycur, N); _acb_vec_clear(Yprev, N); _acb_vec_clear(Ynew, N);
    return 1;
}
static double acb_ub(const acb_t z) { arf_t u; arf_init(u); acb_get_abs_ubound_arf(u, z, 64); double d = arf_get_d(u, ARF_RND_UP); arf_clear(u); return d; }
/* Q (acb_mat, point entries), rho -> Q', rho' ;  Phi column-major vec; tau local errors */
static int lohner_update(acb_mat_t Q, double *rho, acb_srcptr Phi, const double *tau)
{
    acb_mat_t P, M, Qn, Qi, A; acb_mat_init(P, NW, NW); acb_mat_init(M, NW, NW); acb_mat_init(Qn, NW, NW); acb_mat_init(Qi, NW, NW); acb_mat_init(A, NW, NW);
    for (int i = 0; i < NW; i++) for (int j = 0; j < NW; j++) acb_set(acb_mat_entry(P, i, j), Phi + j * NW + i);
    acb_mat_mul(M, P, Q, PREC);
    /* QR (modified Gram-Schmidt, twice) of mid(M) with columns sorted by rho-weighted norm */
    double complex Am[NW][NW]; double w[NW]; int ord[NW]; double rmax = 0;
    for (int j = 0; j < NW; j++) if (rho[j] > rmax) rmax = rho[j];
    for (int j = 0; j < NW; j++) {
        double s = 0;
        for (int i = 0; i < NW; i++) { acb_srcptr z = acb_mat_entry(M, i, j);
            Am[i][j] = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR) + I * arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
            s += creal(Am[i][j] * conj(Am[i][j])); }
        w[j] = sqrt(s) * (rmax > 0 ? fmax(rho[j], 1e-300 * rmax) : 1.0); ord[j] = j;
    }
    for (int a = 0; a < NW; a++) for (int b = a + 1; b < NW; b++) if (w[ord[b]] > w[ord[a]]) { int x = ord[a]; ord[a] = ord[b]; ord[b] = x; }
    double complex V[NW][NW];
    for (int k = 0; k < NW; k++) {
        for (int i = 0; i < NW; i++) V[i][k] = Am[i][ord[k]];
        for (int pass = 0; pass < 2; pass++)
            for (int q = 0; q < k; q++) { double complex d = 0; for (int i = 0; i < NW; i++) d += conj(V[i][q]) * V[i][k];
                for (int i = 0; i < NW; i++) V[i][k] -= d * V[i][q]; }
        double s = 0; for (int i = 0; i < NW; i++) s += creal(V[i][k] * conj(V[i][k])); s = sqrt(s);
        if (!(s > 1e-280) || !isfinite(s)) { for (int i = 0; i < NW; i++) V[i][k] = (i == k); }   /* degenerate: unit vector */
        else for (int i = 0; i < NW; i++) V[i][k] /= s;
    }
    for (int i = 0; i < NW; i++) for (int k = 0; k < NW; k++) acb_set_d_d(acb_mat_entry(Qn, i, k), creal(V[i][k]), cimag(V[i][k]));
    int ok = acb_mat_inv(Qi, Qn, PREC);
    if (!ok) {     /* fall back: Q' = identity (always valid) */
        acb_mat_one(Qn); acb_mat_one(Qi);
    }
    acb_mat_mul(A, Qi, M, PREC);
    double rn[NW];
    for (int i = 0; i < NW; i++) {
        double s = 0;
        for (int j = 0; j < NW; j++) { if (rho[j] > 0) s += acb_ub(acb_mat_entry(A, i, j)) * rho[j];
                                       if (tau[j] > 0) s += acb_ub(acb_mat_entry(Qi, i, j)) * tau[j]; }
        rn[i] = s * (1 + 1e-12);
    }
    for (int i = 0; i < NW; i++) rho[i] = rn[i];
    acb_mat_set(Q, Qn);
    acb_mat_clear(P); acb_mat_clear(M); acb_mat_clear(Qn); acb_mat_clear(Qi); acb_mat_clear(A);
    return 1;
}
