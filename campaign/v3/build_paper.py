"""Build both linked v3 documents, resolving their cross-document labels."""
from pathlib import Path
import argparse, os, shutil, subprocess
ROOT=Path(__file__).resolve().parents[2]
cli=argparse.ArgumentParser(description=__doc__)
cli.add_argument('--engine',default=os.environ.get('TECTONIC','tectonic'))
cli.add_argument('--only-cached',action='store_true')
a=cli.parse_args();engine=shutil.which(a.engine)
if not engine:raise SystemExit('Tectonic not found; set TECTONIC or pass --engine /path/to/tectonic')
paper=ROOT/'paper/v3';report=ROOT/'reports/v3';report.mkdir(exist_ok=True,parents=True)
with (report/'paper_build.log').open('w') as log:
 for name in ('main','supplement','main','supplement'):
  cmd=[engine,'--keep-logs','--keep-intermediates']
  if a.only_cached:cmd+=['--only-cached']
  subprocess.run([*cmd,name+'.tex'],cwd=paper,stdout=log,stderr=subprocess.STDOUT,check=True)
for name in ('main','supplement'):
 s=(paper/(name+'.log')).read_text()
 errors=[x for x in s.splitlines() if 'undefined' in x.lower() or 'multiply defined' in x.lower() or 'Overfull' in x]
 if errors:
  print(name+':');print('\n'.join(errors))
  raise SystemExit('Final build has unresolved references or overflowing boxes.')
 print(paper/(name+'.pdf'))
