"""pares (alpha_k, sigma_k): F_A(alpha) >= sigma_k para todo alpha in [alpha_k, 7] (a certificar em fa_cert.py);
   alpha >= 7: F_A >= 0.796 alpha^2 (certificado, agente Arb v4).  alpha_max(sigma) para (alpha,sigma) em K'."""
import math
CA = 0.796; A_CERT = 7.0; A0P = 19.0
# sigma_k = floor(F_A(alpha_k)*(1-2e-4), 4 decimais) -- valores de thmA.F com folga
PAIRS = [(0.0, 3.7498), (0.25, 3.7787), (0.5, 3.8772), (0.75, 4.0596), (1.0, 4.3406), (1.25, 4.7320), (1.5, 5.2395),
         (1.75, 5.8629), (2.0, 6.5980), (2.25, 7.4399), (2.5, 8.3840), (2.75, 9.4271), (3.0, 10.5672), (3.25, 11.8035),
         (3.5, 13.1355), (3.75, 14.5633), (4.0, 16.0870), (4.25, 17.7071), (4.5, 19.4238), (4.75, 21.2375), (5.0, 23.1484),
         (5.25, 25.1569), (5.5, 27.2632), (5.75, 29.4677), (6.0, 31.7703), (6.25, 34.1713), (6.5, 36.6710), (6.75, 39.2695),
         (7.0, 39.0040)]   # ultimo: c_A*49 = 39.004 (alpha >= 7)
SIGMA_MIN = PAIRS[0][1]
def alpha_max(s_hi):
    """cota superior de alpha para (alpha, sigma) em K' com sigma < s_hi"""
    if s_hi > 39.004: return min(A0P, math.sqrt(s_hi/CA))
    for a, s in PAIRS:
        if s >= s_hi: return a
    return A_CERT
