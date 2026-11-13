import unittest
import numpy as np
from .susceptibility import conditional_moments, moments_summary

class SusceptibilityTests(unittest.TestCase):
    def test_threshold_crossing_integral(self):
        # t=0,1,2; rho rises linearly through threshold at1. Recorded moments
        # m²=.0,.25,1 and m⁴=0,.0625,1 have exact trapezoidal integrals on[1,2].
        ts=np.array([[0,.5,0,0,0],[1,.75,.5,0,0],[2,1,1,0,0]])
        m,f=conditional_moments(ts,0,2,.5)
        np.testing.assert_allclose(m,[.625,.53125,.75]);self.assertAlmostEqual(f,.5)
    def test_no_retained_time(self):
        ts=np.array([[0,.5,0,0,0],[2,.5,0,0,0]])
        m,f=conditional_moments(ts,0,2,.1)
        self.assertTrue(np.isnan(m).all());self.assertEqual(f,0)
    def test_independent_runs_equal_weight(self):
        rows=[dict(Lx=10,Ly=10,T=1000,absorbed=False,active_fraction=1.,
                   all_m2=x,all_m4=x*x,all_abs_m=x**.5) for x in [.1,.9]]
        r=moments_summary(rows,'all',0,np.random.default_rng(1))
        self.assertAlmostEqual(r['chi'],50.)
        self.assertAlmostEqual(r['binder'],1-.41/(3*.5**2))
if __name__=='__main__':unittest.main()
