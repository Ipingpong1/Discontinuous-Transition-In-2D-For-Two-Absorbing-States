from pathlib import Path
import json,subprocess,hashlib,os
root=Path('/home/mo2273/dn_v2_stage/cpu/legacy_replay');m=json.loads((root/'manifest.json').read_text())
assert hashlib.sha256(Path(m['binary']).read_bytes()).hexdigest()==m['binary_sha256']
for t in m['tasks']:
 subprocess.run([m['binary']]+t['args'],check=True,env=dict(os.environ,OMP_NUM_THREADS='1'))
