import os as _os  # [repo] paths relative to the repository root (replaces absolute paths of the original runs)
_R = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
_CERT = _os.path.join(_R, 'certificates'); _WORK = _os.environ.get('CAP_WORK', _os.path.join(_R, 'work'))
def _wout(*p):
    f = _os.path.join(_WORK, *p); _os.makedirs(_os.path.dirname(f), exist_ok=True); return f
# negative controls: the verifier must NOT certify a slightly false target (stops at first unresolved box)
import sys; sys.path.insert(0, _os.path.join(_R, 'g0_sector_lemma'))
from g0lib import *
def attempt(t1,t2,r1,r2,eta,th1,th2,mmin,budget=3000):
    st=[(t1,t2,r1,r2)]; n=0
    while st:
        a1,a2,b1,b2=st.pop(); n+=1
        if n>budget: return 'NOT CERTIFIED (budget exhausted after %d boxes)'%n
        d=box_affine(a1,a2,b1,b2,eta); ok,_=check_box(d,th1,th2,mmin)
        if ok: continue
        if a2-a1<1e-6 and b2-b1<1e-6: return 'NOT CERTIFIED (unresolvable box t=[%.9f,%.9f], rho=[%.2e,%.2e])'%(a1,a2,b1,b2)
        if (d and d['st']>=d['sr']) or b2-b1<1e-9: m=(a1+a2)/2; st+=[(a1,m,b1,b2),(m,a2,b1,b2)]
        else: m=(b1+b2)/2; st+=[(a1,a2,b1,m),(a1,a2,m,b2)]
    return 'CERTIFIED (!!)'
if __name__ == "__main__":
  for args,lab in [((-3.6,-3.5,0,0.05,'0','-0.0835','0.53803',None),'eta=0, th2=0.53803 < max 0.5380375'),
                 ((-1.9,-1.8,0,0.05,'0','-0.0835','0.53805',ub(arb('0.6338'))),'eta=0, mmin=0.6338 > min 0.633763'),
                 ((-3.65,-3.55,0,0.05,'0.75','-0.0935','0.8052',None),'eta=.75, th2=0.8052 < max 0.805268'),
                 ((2.55,2.65,0,0.05,'0.75','-0.0934','0.8054',None),'eta=.75, th1=-0.0934 > min -0.0934067')]:
      print(lab,'->',attempt(*args),flush=True)
