import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from .gpu_repair import verify_replay

class ReplayTests(unittest.TestCase):
    def fixture(self,root):
        stage=root/'new/task1';evidence=root/'evidence/task1'
        stage.mkdir(parents=True);evidence.mkdir(parents=True)
        header='# Lx 64 Ly 64 T 100 seed 7 threads 1\n# rule B3678/S34678 alpha 0.5 voter 0.5 jflip 2:.04592\n'
        times=[0,1,3,10,32,50,100]
        for directory,records in [(stage,times),(evidence,[t for t in times if t!=50])]:
            (directory/'out.txt').write_text(header+''.join(f'{t} .5 .5 0 1\n' for t in records)+'# done t_end 100 absorbed 0\n')
            (directory/'done.json').write_text(json.dumps(dict(id='task1',t_end=100,absorbed=False,final_sha256='same')))
        dtype=np.dtype([('t','<i8'),('alive','<u4',(64,)),('iface','<u4',(64,))])
        profile=np.zeros(len(times),dtype=dtype);profile['t']=times;profile['alive']=32;profile['iface']=64
        profile.tofile(stage/'profile.bin');profile.tofile(evidence/'profile.bin')
        (evidence/'repair_request.json').write_text(json.dumps(dict(expected_final_sha256='same',evidence=[])))
        return stage,evidence,dict(id='task1',profile=True,mode='band',Lx=64,Ly=64),dict(final_sha256='same')
    def test_missing_text_filled_and_all_profiles_match(self):
        with tempfile.TemporaryDirectory() as directory:
            stage,evidence,task,validation=self.fixture(Path(directory))
            report=verify_replay(stage,evidence,task,validation)
            self.assertEqual(report['filled_text_times'],[50])
            self.assertEqual(report['overlapping_text_records'],6)
            self.assertEqual(report['overlapping_profile_records'],7)
            self.assertEqual(report['filled_profile_times'],[])
            dtype=np.dtype([('t','<i8'),('alive','<u4',(64,)),('iface','<u4',(64,))])
            old=np.fromfile(evidence/'profile.bin',dtype=dtype)
            old[old['t']!=50].tofile(evidence/'profile.bin')
            report=verify_replay(stage,evidence,task,validation)
            self.assertEqual(report['filled_profile_times'],[50])
            self.assertEqual(report['overlapping_profile_records'],6)
    def test_profile_or_final_state_mismatch_blocks_replay(self):
        with tempfile.TemporaryDirectory() as directory:
            stage,evidence,task,validation=self.fixture(Path(directory))
            with self.assertRaisesRegex(ValueError,'final state'):
                verify_replay(stage,evidence,task,dict(final_sha256='other'))
            with (stage/'profile.bin').open('r+b') as stream:
                stream.seek(8);stream.write((31).to_bytes(4,'little'))
            with self.assertRaisesRegex(ValueError,'profile counts'):
                verify_replay(stage,evidence,task,validation)

if __name__=='__main__':unittest.main()
