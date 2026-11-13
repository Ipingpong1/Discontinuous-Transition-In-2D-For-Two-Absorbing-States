import unittest
import numpy as np
from .droplets import critical_radius, radius_bootstrap

class DropletTests(unittest.TestCase):
    def test_logit_interpolation_and_exact_half(self):
        self.assertAlmostEqual(critical_radius([2,8],[.2,.8]),4.)
        self.assertAlmostEqual(critical_radius([2,8],[.5,.8]),2.)
    def test_unbracketed_or_ambiguous_radius_is_not_extrapolated(self):
        self.assertTrue(np.isnan(critical_radius([2,8],[.6,.8])))
        self.assertTrue(np.isnan(critical_radius([1,2,3,4],[.2,.8,.2,.8])))
    def test_replica_uncertainty_shrinks_with_larger_fixed_ensemble(self):
        a=radius_bootstrap([2,8],[.2,.8],[100,100],boot=1000)
        b=radius_bootstrap([2,8],[.2,.8],[10000,10000],boot=1000)
        self.assertEqual(a['bootstrap_valid'],1000)
        self.assertLess(b['ci_high']-b['ci_low'],.2*(a['ci_high']-a['ci_low']))

if __name__=='__main__':unittest.main()
