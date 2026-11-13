import unittest
import numpy as np
from .axial10m_diagnostics import paired_diagnostics,EARLY,LATE,GRID

class Axial10MTests(unittest.TestCase):
    def rows(self):
        rows=[]
        for q in GRID:
            for seed in range(4):
                for window,root in [(EARLY,.04594),(LATE,.045942)]:
                    rows.append(dict(T=1e7,L=8192,ab='1-0',estimator='envelope',ensemble_complete=True,
                        threshold=.3,q=q,run=f'{q}_{seed}',t0=window[0],t1=window[1],
                        t_start=window[0],t_end=window[1],v=100*(q-root)+(-3e-5,-1e-5,1e-5,3e-5)[seed]))
        return rows
    def test_paired_shift_and_horizon_isolation(self):
        rows=self.rows();rows.append(dict(rows[0],T=5e6,v=100))
        drift,roots=paired_diagnostics(rows,boot=32)
        self.assertEqual(len(drift),3)
        self.assertTrue(all(abs(r['drift_mean']+.0002)<1e-14 for r in drift))
        for row in drift:
            self.assertAlmostEqual(row['drift_ci_low'],row['drift_mean'],places=14)
            self.assertAlmostEqual(row['drift_ci_high'],row['drift_mean'],places=14)
        root=next(r for r in roots if r['threshold']==.3)
        self.assertEqual(root['bootstrap_valid'],32)
        self.assertAlmostEqual(root['difference'],2e-6,places=14)
    def test_complete_seed_and_geometry_gates(self):
        for change in ('missing','geometry'):
            rows=self.rows()
            if change=='missing':rows.pop()
            else:rows[-1]['t_end']=9e6
            drift,roots=paired_diagnostics(rows,boot=32)
            self.assertEqual(len(drift),2)
            root=next(r for r in roots if r['threshold']==.3)
            self.assertEqual(root['status'],'incomplete_or_geometry_invalid_paired_grid')
            self.assertNotIn('difference',root)
    def test_unbracketed_grid_not_extrapolated(self):
        rows=self.rows()
        for row in rows:row['v']=abs(row['v'])+.001
        _,roots=paired_diagnostics(rows,boot=32)
        root=next(r for r in roots if r['threshold']==.3)
        self.assertEqual(root['status'],'unbracketed_or_nonunique_reversal')
        self.assertTrue(np.isnan(root['q_late']))

if __name__=='__main__':unittest.main()
