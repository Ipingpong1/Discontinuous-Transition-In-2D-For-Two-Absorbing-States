"""Exact empirical-bootstrap tests; exact enumeration is not exact CI coverage."""
import itertools
import unittest
import numpy as np
from .stats import mean_summary


class FourRunMeanBootstrapTests(unittest.TestCase):
    def test_known_binomial_distribution_is_rng_and_budget_invariant(self):
        # Resampling [0,0,0,1] gives K/4 with K~Binomial(4,1/4).
        # Counts at K=0..4 are 81,108,54,12,1 out of256, so the
        # central95% percentile endpoints are exactly0 and3/4.
        for boot,seed in [(1,3),(19,8),(1000,729)]:
            result=mean_summary([0,0,0,1],boot,np.random.default_rng(seed))
            self.assertEqual(result['mean'],.25)
            self.assertEqual(result['ci_low'],0.)
            self.assertEqual(result['ci_high'],.75)

    def test_arbitrary_and_paired_means_match_all_ordered_resamples(self):
        early=np.array([7.,3.,-2.,10.]);late=np.array([4.,3.5,0.,19.])
        values=late-early
        means=[sum(values[i] for i in indices)/4 for indices in itertools.product(range(4),repeat=4)]
        expected=np.quantile(means,[.025,.975])
        for permuted in (values,values[::-1]):
            result=mean_summary(permuted,7,np.random.default_rng(17))
            np.testing.assert_array_equal([result['ci_low'],result['ci_high']],expected)

    def test_disabled_intervals_other_sizes_and_rng_cursor_are_preserved(self):
        result=mean_summary([0,1,2,3],0)
        self.assertTrue(np.isnan(result['ci_low']))
        rng=np.random.default_rng(42);control=np.random.default_rng(42)
        mean_summary([0,1,2,3],23,rng)
        control.integers(0,4,(23,4))
        np.testing.assert_array_equal(rng.integers(0,2**31,12),control.integers(0,2**31,12))
        values=np.array([0.,2.,9.]);rng=np.random.default_rng(7)
        expected=np.quantile(values[rng.integers(0,3,(64,3))].mean(1),[.025,.975])
        result=mean_summary(values,64,np.random.default_rng(7))
        np.testing.assert_array_equal([result['ci_low'],result['ci_high']],expected)


if __name__=='__main__':unittest.main()
