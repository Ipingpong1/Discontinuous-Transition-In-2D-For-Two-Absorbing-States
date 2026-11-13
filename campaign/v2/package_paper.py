#!/usr/bin/env python3
"""Package the current manuscript and its exact figures; no simulations are run."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import re
import zipfile

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT/'paper/v2'
text = (PAPER/'main.tex').read_text()
files = [PAPER/name for name in ('main.tex','refs.bib','main.bbl','Makefile','README.md','evidence_index.json')]
for name in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', text):
    path = PAPER/name
    if not path.suffix:
        path=path.with_suffix('.pdf')
    if not path.resolve().is_relative_to(PAPER.resolve()) or not path.is_file():
        raise ValueError(f'Missing or external manuscript figure: {path}')
    files.append(path)
files=sorted(set(files))
manifest=[]
archive=PAPER/'paper_v2_sources.zip'
temporary=archive.with_suffix('.zip.tmp')
with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as output:
    for path in files:
        contents=path.read_bytes()
        relative=str(path.relative_to(PAPER))
        manifest.append(dict(path=relative,bytes=len(contents),sha256=hashlib.sha256(contents).hexdigest()))
        output.writestr(relative,contents)
    output.writestr('SOURCE_MANIFEST.json',json.dumps(dict(
        created_utc=datetime.now(timezone.utc).isoformat(),files=manifest,
        scope='Manuscript sources and exact figure PDFs. Analysis, raw data and review remain in the project directories described in README.md.'),indent=2)+'\n')
    output.writestr('SHA256SUMS',''.join(f"{f['sha256']}  {f['path']}\n" for f in manifest))
with zipfile.ZipFile(temporary) as check:
    if check.testzip() is not None: raise ValueError('ZIP integrity check failed')
    for item in manifest:
        if hashlib.sha256(check.read(item['path'])).hexdigest()!=item['sha256']:
            raise ValueError('Packaged content differs from source')
os.replace(temporary,archive)
(PAPER/'paper_v2_sources.sha256').write_text(hashlib.sha256(archive.read_bytes()).hexdigest()+'  '+archive.name+'\n')
print(archive)
print(f'{len(files)} verified files, {archive.stat().st_size} bytes')
