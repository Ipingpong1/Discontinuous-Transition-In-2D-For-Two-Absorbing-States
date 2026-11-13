"""Package a v3 source ZIP, reproducibility ZIP and optional compressed data archive.

No simulations, uploads or job submissions are performed.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,re,subprocess,tarfile,zipfile,os,shutil
ROOT=Path(__file__).resolve().parents[2];PAPER=ROOT/'paper/v3';DEST=ROOT/'releases/v3'
IGNORE={'__pycache__','.mplconfig','.git','.claude'}
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
 return h.hexdigest()
def inputs():
 files=set()
 def add(tree,avoid=set()):
  for f in (ROOT/tree).rglob('*'):
   if f.is_file() and not f.is_symlink() and not (set(f.parts)&IGNORE) and f.name not in avoid:
    files.add(f)
 for tree in ['data/figdata','data/runs','data/snap','data/v2/cpu']:
  add(tree)
 add('data/v2/gpu/runs',{'checkpoint.bin'})
 for tree in ['campaign/molofsky_sim','campaign/strip_p099','campaign/round2/local']:
  for f in (ROOT/tree).rglob('*'):
   if f.is_file() and not f.is_symlink() and not(set(f.parts)&IGNORE) and f.suffix in {'.txt','.csv','.tsv','.json','.bin','.prof','.snap','.ckpt'} and len(f.relative_to(ROOT/tree).parts)>1:files.add(f)
 for tree in ['data/v2/gpu/repair_evidence','data/v2/gpu/repair_verifications']:
  add(tree)
 # No duplicate restart checkpoints for completed dense runs; final.bin remains.
 return sorted(files)
def records(files):
 return [dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in files]
def meta(scope,files):
 return dict(release='v3',date='2026-11-04',scope=scope,files=records(files))
def zipwrite(path,files,relroot,extra):
 tmp=path.with_suffix('.zip.tmp')
 with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files:z.write(p,str(p.relative_to(relroot)))
  for key,value in extra.items():z.writestr(key,value)
 with zipfile.ZipFile(tmp) as z:
  if z.testzip():raise RuntimeError('Archive integrity failure')
 os.replace(tmp,path)
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--data',action='store_true',help='Build immutable input-data archive; otherwise source+code ZIPs')
 ap.add_argument('--zstd',default=shutil.which('zstd'))
 a=ap.parse_args();DEST.mkdir(exist_ok=True,parents=True)
 data=inputs()
 if a.data:
  if not a.zstd:raise SystemExit('zstd is required to package the data.')
  manifest=DEST/'INPUTS_MANIFEST.json'
  m=meta('All locally available publication input trees. Dense final states retained; redundant checkpoint copies, bundle mirrors, recovery snapshots and remote archives omitted. Eligibility is defined by scientific manifests; incomplete legacy records remain explicitly identifiable. Some original wall/roughness snapshots are absent; their archived reductions are included.',data)
  manifest.write_text(json.dumps(m,indent=2)+'\n')
  listing=DEST/'input-files.txt';listing.write_text(''.join(str(p.relative_to(ROOT))+'\n' for p in data+[manifest]))
  output=DEST/'paper_v3_inputs.tar.zst';tmp=output.with_suffix('.zst.tmp')
  print('Packaging',len(data),'inputs,',sum(p.stat().st_size for p in data),'bytes',flush=True)
  proc=subprocess.Popen(['tar','-C',str(ROOT),'-cf','-','-T',str(listing)],stdout=subprocess.PIPE)
  with tmp.open('wb') as f:subprocess.run([a.zstd,'-T4','-3','--no-progress','-c'],stdin=proc.stdout,stdout=f,check=True)
  proc.stdout.close()
  if proc.wait():raise RuntimeError('tar failed')
  os.replace(tmp,output)
  subprocess.run([a.zstd,'-t',str(output)],check=True)
  print(output,output.stat().st_size,flush=True)
  return
 # Sources need only TeX, figures, bibliography, link-label cache and portable Makefile.
 mainfiles={PAPER/'Makefile',PAPER/'README.md',PAPER/'CHANGES.md',PAPER/'evidence_index.json',PAPER/'refs.bib'}
 mainfiles.update(PAPER.glob('*.tex'));mainfiles.update((PAPER/'sections').glob('*.tex'))
 for n in ('main.bbl','supplement.bbl','main.aux','supplement.aux'):mainfiles.add(PAPER/n)
 figure_names=set()
 for f in mainfiles:
  if f.suffix=='.tex':figure_names.update(re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',f.read_text()))
 mainfiles.update(PAPER/(n if Path(n).suffix else n+'.pdf') for n in figure_names)
 mainfiles=sorted(mainfiles)
 m=meta('Standalone typesetting sources and exact figure PDFs for main paper and supplement. No raw numerical data or simulator in this source-only ZIP.',mainfiles)
 for rec in m['files']:rec['path']=str(Path(rec['path']).relative_to('paper/v3'))
 zipwrite(DEST/'paper_v3_sources.zip',mainfiles,PAPER,{'SOURCE_MANIFEST.json':json.dumps(m,indent=2)+'\n','SHA256SUMS':''.join(x['sha256']+'  '+x['path']+'\n' for x in m['files'])})
 data_set=set(data);code=[]
 allowed={'.py','.c','.h','.cu','.cuh','.sh','.slurm','.md','.txt','.json','.csv','.tsv','.bib','.tex','.pdf','.png','.log','.yaml','.yml','.toml'}
 trees=['campaign','paper/v3','reports/v3','reports/v2']
 for tree in trees:
  for f in (ROOT/tree).rglob('*'):
   if not f.is_file() or f.is_symlink() or set(f.parts)&IGNORE or f in data_set:continue
   if 'lit' in f.parts or 'host_validation' in f.parts:continue
   if tree=='campaign' and f.suffix in {'.pdf','.png'}:continue
   if tree=='reports/v3' and (f.name.startswith('v2_before') or 'qa' in f.parts):continue
   if f.suffix in allowed or f.name in ('Makefile',):code.append(f)
 code += [ROOT/'dnsim.c',ROOT/'Makefile',ROOT/'README.md']
 # Original helper figure sources occasionally used by the retained reducers.
 code += [f for f in (ROOT/'paper/figs').glob('*.py')]
 code=sorted(set(code))
 manifest=meta('Code, fixed scientific reduction tables, revised manuscript/figures and verification records. Extract companion input archive into this same root for raw-input regeneration. No executable binaries or cluster submissions are needed for analysis.',code)
 readme=(ROOT/'campaign/v3/README.md').read_text()
 zipwrite(DEST/'paper_v3_reproducibility.zip',code,ROOT,{'REPRODUCE.md':readme,'REPRODUCIBILITY_MANIFEST.json':json.dumps(manifest,indent=2)+'\n','SHA256SUMS':''.join(x['sha256']+'  '+x['path']+'\n' for x in manifest['files'])})
 print('Created source and reproducibility archives.')
 for path in sorted(DEST.glob('paper_v3_*')):
  if path.is_file():print(path,path.stat().st_size)
 checks=[p for p in DEST.glob('paper_v3_*') if p.suffix in {'.zip','.zst'}]
 (DEST/'SHA256SUMS').write_text(''.join(sha(p)+'  '+p.name+'\n' for p in sorted(checks)))
if __name__=='__main__':main()
