"""Testes da ancoragem do modo continuacao (V2 da Prop. 5.13): union-find sobre adjacencia por FACE (comprimento > 0)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..'))
import sweep as S
ok = True
def chk(name, cond):
    global ok; ok &= bool(cond); print('%-72s %s' % (name, 'OK' if cond else 'FALHA'))
L, A, N = (0, 1, 10, 11, 'line'), (1, 2, 10, 11), (2, 3, 10, 11, 'none')
# (a) a reta so' alcanca a ancora atraves de uma caixa descartada pelo Teorema A (A nao e' no do grafo) -> rejeitada
chk('(a) ponte so pelo Teorema A -> reta REJEITADA', 0 not in S.anchor_components([L, N]))
# (b) casos validos
chk('(b1) reta com face comum a uma caixa NONE -> aceita', 0 in S.anchor_components([L, (1, 2, 10, 11, 'none')]))
chk('(b2) reta - reta - contorno (transitivo) -> ambas aceitas',
    S.anchor_components([L, (1, 2, 10, 11, 'line'), (2, 3, 10, 11, 'contour')]) == {0, 1})
chk('(b3) face comum PARCIAL (comprimento > 0) -> aceita', 0 in S.anchor_components([L, (1, 2, 10.5, 12, 'none')]))
chk('(b4) sobreposicao de area (reta dentro de regiao NONE) -> aceita', 0 in S.anchor_components([L, (0.5, 3, 9, 12, 'none')]))
chk('(b5) face horizontal (sigma comum) -> aceita', 0 in S.anchor_components([L, (0.5, 1.5, 11, 12, 'none')]))
# (c) contato so' por canto NAO conta
chk('(c) contato so por canto -> rejeitada', 0 not in S.anchor_components([L, (1, 2, 11, 12, 'none')]))
# (d) componente sem ancora: duas retas isoladas
chk('(d) componente so de retas -> ambas rejeitadas', S.anchor_components([L, (1, 2, 10, 11, 'line')]) == set())
sys.exit(0 if ok else 1)
