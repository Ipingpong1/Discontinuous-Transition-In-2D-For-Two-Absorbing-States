import unittest
import numpy as np
from .spatial_correlations import correlator

class SpatialCorrelationTests(unittest.TestCase):
    def test_against_independent_site_rolls(self):
        rng=np.random.default_rng(7193)
        for Lx,Ly in [(64,130),(128,192),(192,256)]:
            states=rng.integers(0,2,(Ly,Lx),dtype=np.uint8)
            words=np.packbits(states,axis=1,bitorder='little').view('<u8')
            distances=sorted({0,1,2,31,min(Lx,Ly)//2}|({63,64,65} if Lx>=192 else set()))
            totals,records=correlator(words,distances)
            self.assertEqual(totals['alive'],int(states.sum()))
            spins=states.astype(float)*2-1
            for rec in records:
                r=rec['r']
                dx=np.count_nonzero(states!=np.roll(states,-r,axis=1))
                dy=np.count_nonzero(states!=np.roll(states,-r,axis=0))
                self.assertEqual(rec['different_x'],dx)
                self.assertEqual(rec['different_y'],dy)
                self.assertAlmostEqual(rec['C'],.5*np.mean(spins*(np.roll(spins,-r,axis=1)+np.roll(spins,-r,axis=0))))
    def test_ordered_configurations_are_retained(self):
        for value in (0,np.iinfo(np.uint64).max):
            totals,records=correlator(np.full((128,2),value,dtype=np.uint64),[0,1,63,64])
            self.assertTrue(all(r['C']==1 and r['C_centered']==0 for r in records))
            self.assertEqual(abs(totals['magnetization']),1)

if __name__=='__main__': unittest.main()
