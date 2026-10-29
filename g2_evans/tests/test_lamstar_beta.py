"""Testes do leitor do formato combinado lamstar_beta (arquivo de EXEMPLO tests/lamstar_beta_exemplo.txt)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
import sweep as S
S.load_lamstar(os.path.join(here, '..', 'lamstar.txt'))
S.load_lamstar_beta(os.path.join(here, 'lamstar_beta_exemplo.txt'))
c = lambda b, s: b*s**(1/3)          # alpha correspondente a beta em sigma
ok = True
def chk(name, cond):
    global ok; ok &= bool(cond); print('%-60s %s' % (name, 'OK' if cond else 'FALHA'))
# 1. caixa inteira em NONE (pior canto: alpha1 * sigma2^{-1/3} >= 0.8 e alpha2 sigma1^{-1/3} <= 2)
a1, a2 = c(0.85, 20), c(1.9, 10)
chk('caixa inteira em NONE -> sai do G2', S.lam_beta_none(a1, a2, 10, 20) and S.lam_box(a1, a2, 10, 20) == 0.0)
# 2. caixa tocando a fronteira beta=0.8 por baixo -> nao sai; Lambda* = 3 (retangulo b)
a1 = c(0.75, 20)
chk('caixa cruzando beta=0.8 nao sai', not S.lam_beta_none(a1, a2, 10, 20))
chk('  e Lambda* = 3 (parte nao-NONE so no retangulo b)', S.lam_box(a1, a2, 10, 20) == 3)
# 3. sobreposicao: beta em [0.45,0.55] -> min(6,3,2) = 2
chk('sobreposicao vale o menor: beta em [0.42,0.58] -> 2', S.lam_box(c(0.42, 10.001), c(0.58, 10), 10, 10.001) == 2)
chk('superconjunto em beta (caixa larga em sigma) e conservador -> 6', S.lam_box(c(0.42, 10), c(0.58, 20), 10, 20) == 6)
# 4. buraco: sigma em [20,30], beta ate 1.5 -> nao coberto -> fallback em sigma (lamstar.txt) e max com a parte coberta
L = S.lam_box(c(0.5, 20), c(1.5, 30), 20, 30)
chk('buraco -> fallback para a tabela em sigma (>= 5)', L is not None and L >= 5 and L == max(5, S.lam_star(20, 30)))
# 5. certificados de exemplo inexistentes -> verificador recusa todos
chk('verificador recusa retangulos sem certificado', len(S.check_lamb_certs()) == 5)
# 6. sem arquivo -> fallback total
S.load_lamstar_beta(os.path.join(here, 'nao_existe.txt'))
chk('arquivo ausente -> Lambda* da tabela em sigma', S.lam_box(1, 2, 10, 11) == S.lam_star(10, 11))
sys.exit(0 if ok else 1)
