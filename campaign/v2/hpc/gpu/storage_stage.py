"""Durable per-task bundles; active simulation files stay on node-local storage."""
from pathlib import Path
import json
import math
import os
import shutil
import struct
import tarfile


def grid_time(path, task):
    p=Path(path)
    if not p.exists() or p.stat().st_size != 32+task['Lx']*task['Ly']//8:
        return None
    with p.open('rb') as f: h=f.read(32)
    if h[:8] != b'DNSIM1\0\0': return None
    x,y,t=struct.unpack('<QQQ',h[8:])
    return t if (x,y)==(task['Lx'],task['Ly']) and t<=task['T'] else None


def members(directory):
    # Checkpoint must be captured first. Text/profile prefixes were flushed
    # before its atomic rename; later appended records are trimmed on resume.
    d=Path(directory)
    ordered=[d/x for x in ('checkpoint.bin','out.txt','profile.bin','initial.bin','final.bin','done.json')]
    ordered += sorted(d.glob('run_*.json'))+sorted(d.glob('process_*.log'))
    return [p for p in ordered if p.is_file() and not p.is_symlink()]


def bundle(directory, destination):
    destination=Path(destination); destination.parent.mkdir(parents=True,exist_ok=True)
    tmp=destination.with_name(destination.name+f'.tmp.{os.getpid()}')
    with tarfile.open(tmp,'w:gz',compresslevel=1) as tar:
        for p in members(directory): tar.add(p,arcname=p.name,recursive=False)
    with tmp.open('rb') as stream: os.fsync(stream.fileno())
    os.replace(tmp,destination)
    fd=os.open(destination.parent,os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)


def unpack(source,directory):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    with tarfile.open(source,'r:gz') as tar:
        for item in tar:
            if not item.isfile() or Path(item.name).name != item.name or item.name in ('.','..'):
                raise ValueError(f'Unsafe bundle member: {item.name}')
            with tar.extractfile(item) as inp, (directory/item.name).open('wb') as out:
                shutil.copyfileobj(inp,out)


def copy_legacy(source,directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    for p in members(source): shutil.copy2(p,directory/p.name)


def expected_times(task, end):
    """Exact ppd=2 plus periodic schedule used by production_runner_v3."""
    times={0}
    for d in range(11):
        for i in range(2):
            t=math.floor(10.0**(d+i/2.0)+.5)  # C llround, nonnegative
            if t<=end: times.add(t)
    every=int(task['every'])
    if every>0: times.update(range(every,int(end)+1,every))
    if task['T']<=end: times.add(int(task['T']))
    return sorted(times)


def validate_prefix(directory, task, end, complete=False):
    """Validate requested parameters and all observations through a checkpoint.

    Prefix mode permits appended or partially written records after ``end``.
    Completion mode requires all files to end exactly there. No state is changed.
    """
    d=Path(directory); lines=[]; rows=[]; meta={}; footers=[]
    if not isinstance(end,int) or end<0 or end>task['T']:
        raise ValueError('invalid observation endpoint')
    raw=(d/'out.txt').read_text()
    for line in raw.splitlines(keepends=True):
        if line.startswith('#'):
            fields=line.split()
            try:
                if line.startswith('# Lx '):
                    for key in ('Lx','Ly','T','seed'):
                        value=int(fields[fields.index(key)+1])
                        if key in meta and meta[key]!=value: raise ValueError('conflicting header '+key)
                        meta[key]=value
                elif line.startswith('# rule '):
                    for key,tag in (('p','voter'),('alpha','alpha')):
                        meta[key]=float(fields[fields.index(tag)+1])
                    extra=fields[fields.index('jflip')+1].split(',')
                    if len(extra)!=1 or extra[0].split(':')[0]!='2':
                        raise ValueError('unexpected jflip model')
                    meta['q']=float(extra[0].split(':')[1])
                elif line.startswith('# profile_dir '):
                    meta['orientation']=[int(x) for x in fields[2].split(',')]
                elif line.startswith('# done t_end '):
                    footers.append((int(fields[3]),bool(int(fields[5]))))
            except (IndexError,ValueError) as exc:
                raise ValueError('invalid header/footer: '+line.strip()) from exc
            if not line.startswith('# done'): lines.append(line)
            continue
        if not line.strip(): raise ValueError('empty numerical record')
        # A live prefix may include a partial next record beyond its checkpoint.
        try: time_value=float(line.split()[0])
        except (ValueError,IndexError) as exc: raise ValueError('invalid record time') from exc
        if not complete and math.isfinite(time_value) and time_value>end: break
        try: row=[float(x) for x in line.split()]
        except ValueError as exc: raise ValueError('invalid numeric record') from exc
        if len(row)!=5 or not all(math.isfinite(x) for x in row) or not line.endswith('\n'):
            raise ValueError('nonfinite, partial, or malformed numeric record')
        if row[0]!=int(row[0]) or row[0]<0 or row[0]>end:
            raise ValueError('numeric time outside requested prefix')
        if rows and row[0]<=rows[-1][0]: raise ValueError('non-increasing numeric times')
        if any(x<0 or x>1 for x in row[1:]): raise ValueError('fraction outside [0,1]')
        rows.append(row);lines.append(line)
    for key in ('Lx','Ly','T','seed','p','q','alpha'):
        expected=float(task[key]) if key in ('p','q','alpha') else int(task[key])
        if meta.get(key)!=expected: raise ValueError('manifest/header mismatch '+key)
    if task['mode']=='band' and meta.get('orientation')!=list(task['orientation']):
        raise ValueError('manifest/header orientation mismatch')
    if task['mode']!='band' and 'orientation' in meta:
        raise ValueError('unexpected directional profile')
    times=[int(row[0]) for row in rows]
    if not times or times[-1]!=end or times!=expected_times(task,end):
        raise ValueError('missing, extra, or incomplete scheduled text records')
    profile_end=None
    if task['profile']:
        bins=task['Lx'] if task['mode']=='band' else task['Ly']
        stride=8+8*bins; path=d/'profile.bin'; totals=task['Lx']*task['Ly']
        count_format=struct.Struct('<'+'I'*(2*bins))
        with path.open('rb') as stream:
            for index,row in enumerate(rows):
                raw=stream.read(stride)
                if len(raw)!=stride: raise ValueError('short profile prefix')
                t=struct.unpack('<q',raw[:8])[0]
                if t!=int(row[0]): raise ValueError('profile/text timeline mismatch')
                counts=count_format.unpack(raw[8:])
                if any(x>task['Lx'] for x in counts): raise ValueError('profile bin count exceeds capacity')
                if (abs(sum(counts[:bins])/totals-row[1])>6e-10 or
                    abs(sum(counts[bins:])/totals-row[4])>6e-10):
                    raise ValueError('profile/text fraction mismatch')
            profile_end=stream.tell()
            if complete and stream.read(1): raise ValueError('profile has trailing records or bytes')
    return dict(lines=lines,rows=rows,profile_end=profile_end,footers=footers)


def validate_complete(directory,task,expected_end,absorbed):
    """Nonmutating scientific integrity check before publishing done.json."""
    d=Path(directory)
    if expected_end!=task['T'] and not absorbed: raise ValueError('nonterminal completion')
    info=validate_prefix(d,task,int(expected_end),complete=True)
    if info['footers']!=[(expected_end,bool(absorbed))]: raise ValueError('completion footer mismatch')
    if grid_time(d/'final.bin',task)!=expected_end: raise ValueError('invalid final grid header/length/time')
    last=info['rows'][-1]
    consensus=last[1] in (0.,1.) and last[2]==0 and last[4]==0
    if absorbed and not consensus: raise ValueError('absorbed marker without consensus observation')
    return True


def prepare_resume(directory,task):
    """Recover a validated prefix, otherwise restart from the same fixed seed."""
    d=Path(directory);t=grid_time(d/'checkpoint.bin',task)
    try:
        if t is None: raise ValueError('missing or invalid checkpoint')
        info=validate_prefix(d,task,t)
    except (ValueError,OSError,struct.error) as exc:
        for name in ('checkpoint.bin','out.txt','profile.bin','final.bin','done.json'):
            (d/name).unlink(missing_ok=True)
        if grid_time(d/'initial.bin',task)!=0: (d/'initial.bin').unlink(missing_ok=True)
        return dict(resume='fresh_same_seed',invalid_or_missing_prefix=True,reason=str(exc))
    (d/'out.txt').write_text(''.join(info['lines']))
    if info['profile_end'] is not None:
        with (d/'profile.bin').open('r+b') as f:f.truncate(info['profile_end'])
    last=info['rows'][-1]
    absorbed=last[1] in (0.,1.) and last[2]==0 and last[4]==0
    terminal=t==task['T'] or absorbed
    if terminal:
        shutil.copy2(d/'checkpoint.bin',d/'final.bin')
        with (d/'out.txt').open('a') as f:
            f.write(f'# recovered terminal checkpoint after prefix validation\n# done t_end {t} absorbed {int(absorbed)} wall_s 0 cell_updates_per_s 0\n')
        validate_complete(d,task,t,absorbed)
    (d/'done.json').unlink(missing_ok=True)
    return dict(resume='validated_checkpoint',checkpoint_time=t,recovered_terminal=terminal)


def completed_bundle(path,task):
    if not Path(path).is_file():return False
    try:
        with tarfile.open(path,'r:gz') as tar:
            data=json.load(tar.extractfile('done.json'))
        return data.get('id')==task['id'] and (data.get('t_end')==task['T'] or data.get('absorbed') is True)
    except (OSError,ValueError,KeyError,tarfile.TarError):return False
