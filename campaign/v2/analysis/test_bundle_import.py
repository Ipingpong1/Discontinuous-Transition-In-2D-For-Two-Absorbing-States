"""Publication bundles must be complete and must not escape their task directory."""
import hashlib
import io
import json
from pathlib import Path
import struct
import tarfile
import tempfile
import unittest
from .import_gpu_bundles import import_one, safe_unpack


class BundleTests(unittest.TestCase):
    def fixture(self,root):
        task=dict(id='example',family='flagship_bulk',Lx=64,Ly=64,T=100,seed=17,
                  p='.5',alpha='.5',q='.04592',mode='bulk',orientation=None,profile=False,every=50)
        state=b'DNSIM1\0\0'+struct.pack('<QQQ',64,64,100)+bytes([85])*512
        done=dict(id='example',t_end=100,absorbed=False,final_sha256=hashlib.sha256(state).hexdigest())
        text=('# Lx 64 Ly 64 T 100 seed 17 threads 1\n'
              '# rule B3678/S34678 alpha 0.5 voter 0.5 jflip 2:.04592\n'
              '# init random:0.5 bg 0 load - t0 0\n'
              +''.join(f'{t} .5 .5 0 1\n' for t in [0,1,3,10,32,50,100])
              +'# done t_end 100 absorbed 0\n')
        files={'out.txt':text.encode(),'final.bin':state,'done.json':json.dumps(done).encode()}
        bundle=root/'bundles'/'flagship_bulk'/'example.tar.gz';bundle.parent.mkdir(parents=True)
        with tarfile.open(bundle,'w:gz') as archive:
            for name,data in files.items():
                member=tarfile.TarInfo(name);member.size=len(data);archive.addfile(member,io.BytesIO(data))
        return task,bundle
    def test_complete_bundle_import_and_idempotence(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);task,bundle=self.fixture(root);records={}
            self.assertEqual(import_one(bundle,task,root,records),'imported')
            self.assertTrue((root/'runs/flagship_bulk/example/out.txt').exists())
            self.assertEqual(import_one(bundle,task,root,records),'unchanged')
    def test_checkpoint_archive_is_never_publication_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);task,bundle=self.fixture(root)
            target=root/'checkpoints/flagship_bulk/example.tar.gz';target.parent.mkdir(parents=True);bundle.rename(target)
            with self.assertRaisesRegex(ValueError,'Checkpoint'):
                import_one(target,task,root,{})
    def test_traversal_and_links_rejected(self):
        for name,kind in [('../escaped',tarfile.REGTYPE),('./..',tarfile.REGTYPE),('out.txt',tarfile.SYMTYPE)]:
            with self.subTest(name=name,kind=kind), tempfile.TemporaryDirectory() as directory:
                root=Path(directory);bundle=root/'bad.tar.gz';stage=root/'stage';stage.mkdir()
                with tarfile.open(bundle,'w:gz') as archive:
                    member=tarfile.TarInfo(name);member.type=kind
                    if kind==tarfile.SYMTYPE:member.linkname='../escaped'
                    archive.addfile(member)
                with self.assertRaises(ValueError):safe_unpack(bundle,stage)
                self.assertFalse((root/'escaped').exists())
    def test_corrupt_old_marker_is_archived(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);task,bundle=self.fixture(root)
            destination=root/'runs/flagship_bulk/example';destination.mkdir(parents=True)
            (destination/'done.json').write_text('{truncated')
            self.assertEqual(import_one(bundle,task,root,{}),'imported')
            old=list((root/'previous_runs/flagship_bulk').glob('*/done.json'))
            self.assertEqual(len(old),1)
            self.assertEqual(old[0].read_text(),'{truncated')

    def test_repair_overlap_gate_precedes_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);task,bundle=self.fixture(root)
            evidence=root/'repair_evidence/example';evidence.mkdir(parents=True)
            safe_unpack(bundle,evidence)
            original=(evidence/'out.txt').read_text()
            marker=json.loads((evidence/'done.json').read_text())
            (evidence/'repair_request.json').write_text(json.dumps(dict(expected_final_sha256=marker['final_sha256'],evidence=[])))
            (evidence/'out.txt').write_text(original.replace('100 .5 .5','100 .4 .5'))
            with self.assertRaisesRegex(ValueError,'preserved text'):
                import_one(bundle,task,root,{})
            self.assertFalse((root/'runs/flagship_bulk/example').exists())
            (evidence/'out.txt').write_text(original.replace('50 .5 .5 0 1\n',''))
            records={};self.assertEqual(import_one(bundle,task,root,records),'imported')
            self.assertTrue(records['example']['repair_verified'])
            report=json.loads((root/'runs/flagship_bulk/example/repair_verification.json').read_text())
            self.assertEqual(report['filled_text_times'],[50])

    def test_existing_complete_digest_conflict_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);task,bundle=self.fixture(root)
            destination=root/'runs/flagship_bulk/example';destination.mkdir(parents=True)
            (destination/'done.json').write_text(json.dumps(dict(t_end=100,final_sha256='conflicting')))
            with self.assertRaisesRegex(ValueError,'Conflict'):
                import_one(bundle,task,root,{})
            self.assertEqual(json.loads((destination/'done.json').read_text())['final_sha256'],'conflicting')

if __name__=='__main__':unittest.main()
