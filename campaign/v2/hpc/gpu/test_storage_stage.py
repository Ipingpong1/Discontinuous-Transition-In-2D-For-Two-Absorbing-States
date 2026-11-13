#!/usr/bin/env python3
"""Storage integrity regressions; fixtures are temporary, never paper data."""
from pathlib import Path
import json
import shutil
import struct
import tempfile
import unittest
import storage_stage as stage

TASK=dict(id='fixture',Lx=64,Ly=64,T=10,seed=18446744073709551615,
          p='0.5',alpha='0.5',q='0.0459',every=5,profile=True,mode='wall')


def fixture(directory,end=10):
    directory.mkdir(parents=True,exist_ok=True)
    header=(f"# dnsim v1.5\n# Lx 64 Ly 64 T 10 seed {TASK['seed']} threads 2\n"
            "# rule B3678/S34678 alpha 0.5 voter 0.5 jflip 2:.0459\n")
    times=stage.expected_times(TASK,end)
    (directory/'out.txt').write_text(header+''.join(f'{t} 0.5 0.125 0.2 0.25\n' for t in times))
    with (directory/'profile.bin').open('wb') as stream:
        for t in times:
            stream.write(struct.pack('<q',t)+struct.pack('<128I',*([32]*64+[16]*64)))
    (directory/'checkpoint.bin').write_bytes(b'DNSIM1\0\0'+struct.pack('<QQQ',64,64,end)+bytes(512))
    return times


class StorageTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.root=Path(self.temporary.name)
        self.data=self.root/'data'
    def tearDown(self):self.temporary.cleanup()
    def test_good_nonterminal_trims_only_tail(self):
        fixture(self.data,5)
        with (self.data/'out.txt').open('a') as f:f.write('10 .5 .1')
        with (self.data/'profile.bin').open('ab') as f:f.write(b'partial')
        result=stage.prepare_resume(self.data,TASK)
        self.assertEqual(result['checkpoint_time'],5)
        self.assertFalse(result['recovered_terminal'])
        self.assertEqual((self.data/'profile.bin').stat().st_size,4*520)
        self.assertNotIn('10 .5',(self.data/'out.txt').read_text())
    def test_terminal_recovery_and_atomic_bundle_roundtrip(self):
        fixture(self.data)
        result=stage.prepare_resume(self.data,TASK)
        self.assertTrue(result['recovered_terminal'])
        self.assertTrue(stage.validate_complete(self.data,TASK,10,False))
        (self.data/'done.json').write_text(json.dumps(dict(id='fixture',t_end=10,absorbed=False)))
        archive=self.root/'bundles'/'fixture.tar.gz'
        stage.bundle(self.data,archive)
        self.assertTrue(stage.completed_bundle(archive,TASK))
        stage.unpack(archive,self.root/'restored')
        self.assertTrue(stage.validate_complete(self.root/'restored',TASK,10,False))
    def test_profile_gap_rejected(self):
        fixture(self.data)
        path=self.data/'profile.bin';raw=path.read_bytes()
        path.write_bytes(raw[:2*520]+raw[3*520:])
        result=stage.prepare_resume(self.data,TASK)
        self.assertEqual(result['resume'],'fresh_same_seed')
        self.assertIn('timeline',result['reason'])
    def test_matching_gap_in_both_files_rejected(self):
        fixture(self.data)
        path=self.data/'out.txt';path.write_text(path.read_text().replace('3 0.5 0.125 0.2 0.25\n',''))
        path=self.data/'profile.bin';raw=path.read_bytes();path.write_bytes(raw[:2*520]+raw[3*520:])
        result=stage.prepare_resume(self.data,TASK)
        self.assertEqual(result['resume'],'fresh_same_seed')
        self.assertIn('scheduled',result['reason'])
    def test_corrupt_header_and_profile_counts(self):
        for corruption in ('seed','counts'):
            with self.subTest(corruption=corruption):
                if self.data.exists():shutil.rmtree(self.data)
                fixture(self.data)
                if corruption=='seed':
                    path=self.data/'out.txt';path.write_text(path.read_text().replace(str(TASK['seed']),'42'))
                else:
                    path=self.data/'profile.bin';raw=bytearray(path.read_bytes());struct.pack_into('<I',raw,8,0);path.write_bytes(raw)
                self.assertEqual(stage.prepare_resume(self.data,TASK)['resume'],'fresh_same_seed')
    def test_completion_validation_is_nonmutating(self):
        fixture(self.data);stage.prepare_resume(self.data,TASK)
        path=self.data/'profile.bin';path.write_bytes(path.read_bytes()+b'x')
        before={p.name:p.read_bytes() for p in self.data.iterdir()}
        with self.assertRaises(ValueError):stage.validate_complete(self.data,TASK,10,False)
        self.assertEqual(before,{p.name:p.read_bytes() for p in self.data.iterdir()})

if __name__=='__main__':unittest.main()
