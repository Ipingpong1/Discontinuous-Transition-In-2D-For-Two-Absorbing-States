"""Numerical regression tests for scientific estimators, not plot appearance."""
import unittest
import json
import tempfile
from pathlib import Path
import numpy as np
from .bulk import extend_absorbed, reduce_run, summarize
from .fronts import covers_requested_window, crossing_models, polynomial_crossing, velocity, widths, summarize as front_summary
from .io import profile_path, run_id, completion_status, complete_ensembles, read_timeseries
from .stats import kaplan_meier, km_median, restricted_mean, survival_summary, window_mean
from .wetting import fit_models
from .gpu_campaign import task_completion

class FrontTests(unittest.TestCase):
    def test_internal_hole_does_not_stop_envelope(self):
        d=np.zeros((2,64));d[:,:32]=1.;d[1,14:18]=0.
        w=widths(np.array([0,1000]),d,1,threshold=.3,smooth_bins=1,reference_time=0)
        self.assertEqual(w['envelope'][0],32)
        self.assertEqual(w['envelope'][1],32)
        self.assertLess(w['component'][1],32)
        self.assertAlmostEqual(w['hole_fraction'][1],4/32)
    def test_translation_keeps_width(self):
        d=np.zeros((2,64));d[0,:32]=1;d[1,3:35]=1
        w=widths(np.array([0,1000]),d,1,smooth_bins=1)
        np.testing.assert_array_equal(w['envelope'],[32,32])
    def test_normal_velocity(self):
        t=np.arange(20)*1000.;norm=np.sqrt(5);v=.001
        fit,_=velocity(t,1000+2*v*t*norm,4096,norm,t0=0)
        self.assertAlmostEqual(fit['v'],v,places=12)
    def test_curvature_reported(self):
        x=np.array([.99,.995,1.,1.005,1.01]);v=.3*(x-1)+4*(x-1)**2
        f=polynomial_crossing(x,v,np.full(5,1e-6),2)
        self.assertAlmostEqual(f['crossing'],1.,places=12)
        self.assertLess(f['chi2'],1e-12)
        self.assertGreater(polynomial_crossing(x,v,np.full(5,1e-6),1)['reduced_chi2'],100.)
    def test_no_extrapolated_crossing(self):
        self.assertIsNone(polynomial_crossing([1,2,3],[1,2,3],[.1,.1,.1],1))
    def test_bootstrap_bracket_can_switch_pairs(self):
        f=crossing_models([.9,1.,1.1],[[-1.1,-.9,-1.,-1.],[-.1,.1,-.05,.05],[.9,1.1,1.,1.]],boot=100)
        b=next(r for r in f if r['model']=='bracket')
        self.assertEqual(b['bootstrap_valid'],100)
        self.assertLess(b['ci_low'],1.)
        self.assertGreater(b['ci_high'],1.)
    def test_geometrically_truncated_fit_is_not_full_window(self):
        row=dict(v=.001,t_start=1e6,t_end=1.5e6,t0=1e6,t1=2e6,T=2e6)
        self.assertFalse(covers_requested_window(row))
        row['t_end']=2e6
        self.assertTrue(covers_requested_window(row))
    def test_independent_horizons_are_never_pooled_in_shared_window(self):
        rows=[]
        for horizon,zero in [(2e6,.04594),(5e6,.04595)]:
            for q in [.04592,.04594,.04596]:
                for noise in [-3e-5,-1e-5,1e-5,3e-5]:
                    rows.append(dict(p=.5,L=8192,Ly=8192,T=horizon,ab='1-0',threshold=.3,
                        smooth_normal=9.,t0=1e6,t1=2e6,estimator='envelope',q=q,v=10*(q-zero)+noise))
        points,fits=front_summary(rows,boot=20)
        self.assertEqual(len(points),6)
        self.assertTrue(all(r['n']==4 for r in points))
        roots={r['T']:r['crossing'] for r in fits if r['model']=='bracket'}
        self.assertAlmostEqual(roots[2e6],.04594,places=12)
        self.assertAlmostEqual(roots[5e6],.04595,places=12)

    def test_stable_cohort_bootstrap_ignores_unrelated_completion(self):
        rows=[]
        for q in [.04592,.04594,.04596]:
            for noise in [-3e-5,-1e-5,1e-5,3e-5]:
                rows.append(dict(p=.5,L=8192,Ly=8192,T=5e6,ab='1-0',threshold=.3,
                    smooth_normal=9.,t0=1e6,t1=2e6,estimator='envelope',q=q,v=10*(q-.045945)+noise))
        alone=front_summary(rows,boot=100,stable_groups=True)[1]
        combined=front_summary(rows+[dict(r,T=2e6) for r in rows],boot=100,stable_groups=True)[1]
        same=[r for r in combined if r['T']==5e6]
        for a,b in zip(alone,same):
            self.assertEqual(a['model'],b['model'])
            self.assertEqual(a['ci_low'],b['ci_low'])
            self.assertEqual(a['ci_high'],b['ci_high'])

    def test_gpu_nested_profile_names(self):
        self.assertEqual(run_id('/data/orientation/task1/out.txt'),'task1')
        self.assertEqual(str(profile_path('/data/orientation/task1/out.txt')),'/data/orientation/task1/profile.bin')

class CompletionTests(unittest.TestCase):
    def test_profile_lagging_complete_text_is_rejected(self):
        meta = dict(done=True,T=100,t_end=100,absorbed=False)
        ts = np.array([[0,.5,.4,0,.8],[50,.5,.3,0,.7],[100,.5,.3,0,.7]])
        profile = np.array([(0,),(50,)],dtype=[('t','<i8')])
        result = completion_status('legacy.txt',meta,ts,profile,require_marker=False)
        self.assertFalse(result['complete'])
        self.assertIn('profile_end_mismatch',result['reasons'])
    def test_partial_terminal_profile_is_rejected(self):
        meta = dict(done=True,T=100,t_end=100,absorbed=False)
        ts = np.array([[0,.5,.4,0,.8],[100,.5,.3,0,.7]])
        profile = np.array([(0,),(100,)],dtype=[('t','<i8')])
        self.assertFalse(completion_status('legacy.txt',meta,ts,profile,trailing=3,require_marker=False)['complete'])
    def test_gpu_marker_must_match_terminal_record(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'task1'/'out.txt';path.parent.mkdir()
            path.write_text('# Lx 64 Ly 64 T 100 seed 7 threads 1\n# rule B3678/S34678 alpha 0.5 voter 0.5 jflip 2:.04592\n0 .5 .4 0 .8\n100 .5 .3 0 .7\n# done t_end 100 absorbed 0\n')
            (path.parent/'done.json').write_text(json.dumps(dict(id='task1',t_end=90,absorbed=False)))
            meta,_=read_timeseries(path)
            self.assertFalse(meta['done'])
            self.assertIn('completion_marker_mismatch',meta['completion_errors'])
            (path.parent/'done.json').write_text(json.dumps(dict(id='task1',t_end=100,absorbed=False)))
            self.assertTrue(read_timeseries(path)[0]['done'])
    def test_absorption_before_origin_is_not_assigned_planned_horizon(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'early.txt'
            path.write_text('# Lx 64 Ly 64 T 100000 seed 7 threads 1\n# rule B3678/S34678 alpha 0.5 voter 0.5 jflip 2:.04592\n0 .5 .4 0 .8\n100 0 0 0 0\n# done t_end 100 absorbed 1\n')
            row=reduce_run(path)
            self.assertEqual(row['half_time'],1000)
            self.assertEqual(row['half_lower'],0)
            self.assertTrue(row['half_left_censored'])
            self.assertEqual(row['absorption_time'],100)
            self.assertEqual(row['rho_all'],0)
    def test_corrupt_complete_profile_counts_are_rejected(self):
        meta=dict(done=True,T=100,t_end=100,absorbed=False,Lx=64,Ly=64)
        ts=np.array([[0,.5,.4,0,1.],[100,.5,.4,0,1.]])
        profile=np.zeros(2,dtype=[('t','<i8'),('alive','<u4',(64,)),('iface','<u4',(64,))])
        profile['t']=[0,100];profile['alive']=32;profile['iface']=64
        self.assertTrue(completion_status('legacy.txt',meta,ts,profile,require_marker=False)['complete'])
        profile['iface'][1,0]=65
        status=completion_status('legacy.txt',meta,ts,profile,require_marker=False)
        self.assertFalse(status['complete'])
        self.assertIn('profile_count_exceeds_bin_capacity',status['reasons'])
        self.assertIn('profile_text_fraction_mismatch',status['reasons'])
    def test_nonfinite_completed_numeric_output_is_rejected(self):
        meta=dict(done=True,T=100,t_end=100,absorbed=False)
        ts=np.array([[0,.5,.4,0,.8],[100,.5,np.nan,0,.8]])
        self.assertFalse(completion_status('legacy.txt',meta,ts,require_marker=False)['complete'])

    def test_truncated_mean_retains_left_and_right_observation_bounds(self):
        with tempfile.TemporaryDirectory() as d:
            header='# Lx 64 Ly 64 T 100000 seed 7 threads 1\n# rule B3678/S34678 alpha 0.5 voter 0.5 jflip 2:.04592\n'
            a=Path(d)/'early.txt';a.write_text(header+'0 .5 .4 0 .8\n500 0 0 0 0\n# done t_end 500 absorbed 1\n')
            b=Path(d)/'censored.txt';b.write_text(header+'0 .5 .4 0 .8\n100000 .5 .3 0 .7\n# done t_end 100000 absorbed 0\n')
            result=summarize([reduce_run(a),reduce_run(b)],boot=0)[0][0]
            self.assertEqual(result['half_recording_lower_mean'],50000.)
            self.assertEqual(result['half_recording_upper_mean'],50500.)
            self.assertEqual(result['half_rmst'],50500.)
            self.assertEqual(result['half_left_censored'],1)

    def test_completed_gpu_schedule_gaps_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'task1/out.txt';path.parent.mkdir()
            task=dict(id='task1',family='flagship_bulk',Lx=64,Ly=64,T=100,seed=7,
                p='.5',q='.04592',alpha='.5',mode='bulk',orientation=None,profile=False,every=50)
            header='# Lx 64 Ly 64 T 100 seed 7 threads 1\n# rule B3678/S34678 alpha 0.5 voter 0.5 jflip 2:.04592\n# init random:0.5 bg 0 load - t0 0\n'
            (path.parent/'done.json').write_text(json.dumps(dict(id='task1',t_end=100,absorbed=False)))
            times=[0,1,3,10,25,32,50,100]  # Extra finer-cadence sample25 is allowed.
            for missing in (32,50):
                path.write_text(header+''.join(f'{t} .5 .3 0 .7\n' for t in times if t!=missing)+'# done t_end 100 absorbed 0\n')
                status,_,_=task_completion(task,path)
                self.assertFalse(status['complete'])
                self.assertEqual(status['missing_scheduled_times'],[missing])
            path.write_text(header+''.join(f'{t} .5 .3 0 .7\n' for t in times)+'# done t_end 100 absorbed 0\n')
            self.assertTrue(task_completion(task,path)[0]['complete'])

    def test_missing_manifest_seed_keeps_ensemble_ineligible(self):
        task=dict(p='.5',alpha='.5',q='.04592',Lx=64,Ly=64,T=100,mode='bulk',orientation=None)
        tasks={f'task{i}':dict(task,id=f'task{i}') for i in range(4)}
        status={f'task{i}':dict(complete=True) for i in range(3)}
        result=complete_ensembles(tasks,status)[0]
        self.assertEqual(result['n_planned'],4)
        self.assertEqual(result['n_complete'],3)
        self.assertFalse(result['ensemble_complete'])
        status['task3']=dict(complete=True)
        self.assertTrue(complete_ensembles(tasks,status)[0]['ensemble_complete'])


class StatisticsTests(unittest.TestCase):
    def test_time_weighting_irregular_grid(self):
        self.assertAlmostEqual(window_mean([0,1,10],[0,1,10],0,10),5)
        self.assertAlmostEqual(window_mean([0,1,10],[0,1,10],2,8),5)
        self.assertTrue(np.isnan(window_mean([0,1],[0,1],0,2)))
    def test_kaplan_meier_censoring(self):
        c=kaplan_meier([1,2,3],[True,False,True])
        np.testing.assert_allclose([r['survival'] for r in c],[1,2/3,2/3,0])
        self.assertEqual(km_median(c),3)
        self.assertAlmostEqual(restricted_mean(c,3),7/3)
        self.assertGreater(c[-1]['high'],0.)
    def test_tied_death_and_censor(self):
        c=kaplan_meier([1,1,2],[True,False,True])
        self.assertAlmostEqual(c[1]['survival'],2/3)
        self.assertEqual(c[1]['at_risk'],3)
    def test_all_censored_no_fake_precision(self):
        s=survival_summary([10,10,10],[False]*3,boot=20)
        self.assertTrue(np.isnan(s['median']))
        self.assertTrue(np.isnan(s['rmst_low']))
        self.assertLess(s['curve'][-1]['low'],1)
        self.assertEqual(s['rmst'],10)
    def test_absorption_extended_exactly(self):
        ts=np.array([[0,.5,.4,.1,.8],[10,1,0,0,0]])
        out=extend_absorbed(ts,100)
        self.assertEqual(out[-1,0],100)
        self.assertEqual(window_mean(out[:,0],out[:,2],50,100),0)
        np.testing.assert_array_equal(out[-1,1:],[1,0,0,0])

class WettingTests(unittest.TestCase):
    def test_known_power_recovered(self):
        x=np.array([.004,.007,.01,.02,.03,.05]);y=6*x**(-.52)
        f=fit_models(x,y,np.eye(len(x))*.1)
        self.assertAlmostEqual(f['power']['psi'],.52,places=7)
        self.assertAlmostEqual(f['power']['A'],6,places=6)
        self.assertLess(f['power']['chi2'],1e-10)
    def test_fixed_third_offset_recovered(self):
        x=np.array([.004,.007,.01,.02,.03,.05]);y=20*x**(-1/3)-30
        f=fit_models(x,y,np.eye(len(x))*.1+.02)
        self.assertAlmostEqual(f['third_offset']['A'],20,places=7)
        self.assertAlmostEqual(f['third_offset']['B'],-30,places=7)
        self.assertLess(f['third_offset']['chi2'],1e-10)

if __name__=='__main__':unittest.main()
