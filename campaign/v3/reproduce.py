"""Recreate v3 reductions and all publication figures from the supplied project.

The v2 reduction tables are fixed inputs for retained analyses. This command
recomputes revised analyses from their archived raw inputs and then renders
every figure used by either document. It never writes to v2.
"""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[2]
ANALYSES=['bulk_activity','fronts_review','film_endpoint','molofsky']
FIGURES=['fig_snap','fig_flagship','fig_drop2','fig_drop','fig_wetting_extended','fig_phase',
         'fig_gpu_near20m','fig_theory_v2','fig_molofsky','fig_orient','fig_gpu_v2',
         'fig_spatial_v2','fig_bulk05','fig_lifetimes_v2','fig_wall05','fig_rough','fig_voter','fig_eta']
def digest(p):
 return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 cli=argparse.ArgumentParser(description=__doc__)
 cli.add_argument('--figures-only',action='store_true')
 cli.add_argument('--tests',action='store_true')
 a=cli.parse_args();report=ROOT/'reports/v3';report.mkdir(exist_ok=True,parents=True)
 before={str(f.relative_to(ROOT)):digest(f) for t in ['paper/v2','campaign/v2','reports/v2']
         for f in (ROOT/t).rglob('*') if f.is_file() and not {'__pycache__','.mplconfig'}&set(f.parts)}
 env=dict(os.environ,MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1')
 steps=[]
 commands=[]
 if not a.figures_only:
  commands += [[sys.executable,'-B','-m','campaign.v3.analysis.'+m] for m in ANALYSES]
 commands += [[sys.executable,'-B',str(ROOT/'paper/v3/figs'/(m+'.py'))] for m in FIGURES]
 if a.tests:commands += [[sys.executable,'-B','-m','unittest','discover','-s','campaign/v2/analysis','-t','.','-q']]
 with (report/'reproduce.log').open('w') as log:
  for cmd in commands:
   print('Running',cmd[-1],flush=True);start=time.monotonic()
   subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
   steps.append(dict(command=cmd,seconds=round(time.monotonic()-start,2)))
 changed=[p for p,h in before.items() if not (ROOT/p).is_file() or digest(ROOT/p)!=h]
 if changed:raise RuntimeError('v2 files changed: '+repr(changed))
 (report/'reproduce_verification.json').write_text(json.dumps(dict(steps=steps,
  v2_files_checked=len(before),v2_unchanged=True,
  figures=[dict(path=str(f.relative_to(ROOT)),sha256=digest(f)) for f in sorted((ROOT/'paper/v3/figs').glob('*.pdf'))]),indent=2)+'\n')
 print('Reproduction complete; v2 unchanged.',flush=True)
if __name__=='__main__':main()
