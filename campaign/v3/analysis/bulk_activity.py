"""Direct late flip fractions and bond-density checks from the largest bulk runs."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'campaign/v3/results/bulk_activity'
def mean(t,y,a,b):
 z=np.r_[a,t[(t>a)&(t<b)],b]
 return float(np.trapezoid(np.interp(z,t,y),z)/(b-a))
def main():
 rows=[];inputs=[]
 for f in sorted((ROOT/'data/v2/gpu/runs/flagship_bulk').glob('*L16384x16384*/out.txt')):
  d=np.loadtxt(f);assert d[-1,0]==5000000 and d[-1,2]>0
  q=float(f.parent.name.split('_q')[1].split('_s')[0])
  rho=mean(d[:,0],d[:,2],2.5e6,5e6);flip=mean(d[:,0],d[:,3],2.5e6,5e6)
  assert flip>=.25*rho
  rows.append(dict(run=f.parent.name,q=q,L=16384,T=5000000,rho=rho,flip=flip,rho_drift=mean(d[:,0],d[:,2],3.75e6,5e6)-mean(d[:,0],d[:,2],2.5e6,3.75e6)))
  inputs.append(dict(path=str(f.relative_to(ROOT)),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 assert len(rows)==12
 summaries=[]
 for q in sorted({r['q'] for r in rows}):
  g=[r for r in rows if r['q']==q];assert len(g)==4
  summaries.append(dict(q=q,n=4,**{k:float(np.mean([r[k] for r in g])) for k in ['rho','flip','rho_drift']}))
 OUT.mkdir(exist_ok=True,parents=True)
 for name,data in [('runs',rows),('summary',summaries)]:
  with (OUT/(name+'.csv')).open('w') as w:
   z=csv.DictWriter(w,fieldnames=data[0]);z.writeheader();z.writerows(data)
 (OUT/'provenance.json').write_text(json.dumps(dict(inputs=inputs,window=[2500000,5000000],methods='Time integration within each run; equal independent-run weights',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
 print(summaries)
if __name__=='__main__':main()
