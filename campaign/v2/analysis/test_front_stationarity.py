import unittest
import numpy as np
from .front_stationarity import paired_contrasts,bracket_zero,EARLY,LATE,GRID

class StationarityTests(unittest.TestCase):
    def fixture(self):
        rows=[]
        for q in GRID:
            for seed,noise in enumerate([-3e-5,-1e-5,1e-5,3e-5]):
                for window,shift in [(EARLY,0.),(LATE,5e-5)]:
                    rows.append(dict(p=.5,L=8192,Ly=8192,T=5_000_000,ab='1-0',threshold=.3,
                        smooth_normal=9.,estimator='envelope',q=q,run=f'q{q}s{seed}',
                        t0=window[0],t1=window[1],v=10*(q-.04594)+noise+shift))
        return rows
    def test_pairing_preserves_common_run_noise(self):
        drift,crossing=paired_contrasts(self.fixture(),boot=200)
        self.assertEqual(len(drift),3);self.assertEqual(len(crossing),1)
        for row in drift:
            self.assertAlmostEqual(row['drift_mean'],5e-5,places=12)
            self.assertAlmostEqual(row['drift_ci_low'],row['drift_ci_high'],places=12)
        self.assertAlmostEqual(crossing[0]['difference'],-5e-6,places=12)
    def test_no_pooling_of_original_cohort_or_missing_window(self):
        rows=self.fixture()
        original=[dict(r,T=2_000_000,v=9) for r in rows]
        a=paired_contrasts(rows+original,boot=20)
        self.assertAlmostEqual(a[0][0]['drift_mean'],5e-5,places=12)
        rows.pop()
        drift,crossing=paired_contrasts(rows,boot=20)
        self.assertEqual(len(drift),2);self.assertEqual(crossing,[])
    def test_zero_never_extrapolated(self):
        self.assertTrue(np.isnan(bracket_zero([1,2,3],[1,2,3])))
        self.assertTrue(np.isnan(bracket_zero([1,2,3],[1,-1,1])))
        self.assertTrue(np.isnan(bracket_zero([1,2,3],[0,-1,1])))
        self.assertTrue(np.isnan(bracket_zero([1,2,3],[1,0,1])))
        self.assertTrue(np.isnan(bracket_zero([1,2,3],[0,0,1])))
        self.assertEqual(bracket_zero([1,2,3],[-1,0,1]),2.)

if __name__=='__main__':unittest.main()
