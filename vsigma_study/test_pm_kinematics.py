"""Scientific regression checks: python -m unittest vsigma_study.test_pm_kinematics."""
import unittest
from unittest.mock import patch
import json
import numpy as np
from .pm_kinematics import basis, transform_pm, select_bins, fit_bin, K


class GeometryTests(unittest.TestCase):
    def test_bulk_and_perspective_removed_with_covariance(self):
        theta = np.linspace(0, 2*np.pi, 1000, endpoint=False)
        n0, e0, north0 = basis(201., -47.)
        angle = .006
        n = np.cos(angle)*n0 + np.sin(angle)*(np.cos(theta)[:,None]*e0+np.sin(theta)[:,None]*north0)
        ra = np.degrees(np.arctan2(n[:,1], n[:,0])) % 360
        dec = np.degrees(np.arcsin(n[:,2]))
        n, e, north = basis(ra, dec)
        radial = (np.cos(angle)*n-n0)/np.sin(angle)
        tangent = np.cross(n, radial)
        velocity = 5*K*(-3*e0-6*north0)+232*n0 + .2*5*K*tangent
        pm = np.column_stack([np.sum(velocity*e, axis=1), np.sum(velocity*north, axis=1)])/(5*K)
        r, u, cov, bulk = transform_pm(ra,dec,pm[:,0],pm[:,1],np.full((1000,2),.1),
                                      np.zeros(1000),(201.,-47.),5.,232.)
        np.testing.assert_allclose(u[:,0], 0, atol=1e-9)
        np.testing.assert_allclose(u[:,1], .2, atol=1e-9)
        np.testing.assert_allclose(bulk, [-3,-6], atol=1e-9)
        np.testing.assert_allclose(cov, np.tile(np.eye(2)*.01,(1000,1,1)), atol=1e-10)

    def test_bins_cover_every_star_and_use_error_information(self):
        r = np.linspace(.1, 12, 5000)
        cov = np.tile(np.eye(2)*.02**2, (len(r),1,1))
        good, _, _ = select_bins(r,cov,3.)
        poor, _, diagnostics = select_bins(r,cov*1000,3.)
        self.assertLess(len(poor),len(good))
        np.testing.assert_array_equal(np.sort(np.concatenate(poor)),np.arange(len(r)))
        self.assertTrue(any(d['low_information'] for d in diagnostics))
        for mode in ('equal_number','equal_radius'):
            bins, _, _ = select_bins(r,cov,3.,mode,17)
            self.assertEqual(len(bins),17)
            self.assertEqual(sum(map(len,bins)),len(r))


class InferenceTests(unittest.TestCase):
    def simulate(self, dispersions, errors, n=800, seed=17):
        rng=np.random.default_rng(seed)
        e=rng.uniform(errors*.6,errors*1.4,size=(n,2))
        c=np.zeros((n,2,2)); c[:,0,0]=e[:,0]**2; c[:,1,1]=e[:,1]**2
        c[:,0,1]=c[:,1,0]=.35*e[:,0]*e[:,1]
        intrinsic=np.diag(np.array(dispersions)**2)
        noise=np.einsum('nij,nj->ni',np.linalg.cholesky(c+intrinsic),rng.normal(size=(n,2)))
        return noise+[0.,.2], c

    def test_recover_rotation_and_anisotropic_dispersion(self):
        u,c=self.simulate([.13,.08],.04)
        fit=fit_bin(u,c,max_steps=5000)
        self.assertEqual(fit['status'],'resolved',fit)
        for key, truth in [('mean_t',.2),('sigma_r',.13),('sigma_t',.08)]:
            lo,mid,hi=fit[key]
            self.assertLess(abs(mid-truth),3*(hi-lo),fit)

    def test_pure_rotation_does_not_create_intrinsic_dispersion(self):
        u,c=self.simulate([0.,0.],.1,seed=32)
        fit=fit_bin(u,c,max_steps=5000)
        self.assertEqual(fit['status'],'unresolved',fit)
        self.assertLess(fit['sigma_upper95'],.06)
        self.assertAlmostEqual(fit['mean_t'][1],.2,delta=.02)

    def test_truncated_sample_flagged(self):
        u,c=self.simulate([0.,0.],.15,n=2000)
        keep=np.all(np.abs(u-[0.,.2])<.04,axis=1)
        fit=fit_bin(u[keep],c[keep],max_steps=5000)
        self.assertEqual(fit['status'],'error_model_conflict',fit)

    def test_zero_rotation_is_not_detected_by_absolute_value(self):
        u,c=self.simulate([.12,.12],.04,seed=29)
        u[:,1] -= .2
        fit=fit_bin(u,c,max_steps=5000)
        self.assertTrue(fit['converged'],fit)
        self.assertGreater(min(fit['positive_probability'],1-fit['positive_probability']),.0015)


class IntegrationTests(unittest.TestCase):
    def test_legacy_missing_covariance_is_counted_and_json_is_finite(self):
        import pandas as pd
        from .pm_kinematics import analyze
        from .vsigma_pipeline import CLUSTERS
        rng = np.random.default_rng(12)
        df = pd.DataFrame(dict(ra=201.6968+rng.normal(0,.1,100), dec=-47.4796+rng.normal(0,.1,100),
                               pmra=rng.normal(-3,.1,100), pmdec=rng.normal(-6,.1,100),
                               pmra_error=.05, pmdec_error=.05, membership_prob=1.))
        with patch('pandas.read_csv',return_value=df), patch('vsigma_study.pm_kinematics.fit_bin',return_value=dict(status='not_converged',converged=False)):
            result = analyze('NGC5139','unused.csv',CLUSTERS['NGC5139'])
        self.assertEqual(result['diagnostics']['diagonal_covariance_stars'],100)
        self.assertEqual(result['n_members'],100)
        self.assertEqual(result['bins']['ratio'],[None])
        json.dumps(result,allow_nan=False)

    def test_rotation_api_defaults_and_validation(self):
        import server
        import threading
        import urllib.request
        from urllib.error import HTTPError
        with server.ReusableServer(('127.0.0.1',0),server.Handler) as httpd:
            worker=threading.Thread(target=httpd.serve_forever,daemon=True)
            worker.start()
            opener=urllib.request.build_opener(urllib.request.ProxyHandler({}))
            base=f'http://127.0.0.1:{httpd.server_address[1]}'
            try:
                with patch('vsigma_study.vsigma_pipeline.analyze_cluster',return_value={'model_version':'test'}) as fit:
                    result=json.load(opener.open(base+'/rotation?id=NGC5139'))
                    self.assertEqual(result['catalog_variant'],'dr3')
                    self.assertEqual(fit.call_args.kwargs['bin_mode'],'auto')
                    self.assertIsNone(fit.call_args.kwargs['nbins'])
                with self.assertRaises(HTTPError) as error:
                    opener.open(base+'/rotation?id=NGC5139&variant=invalid')
                self.assertEqual(error.exception.code,400)
            finally:
                httpd.shutdown()
                worker.join()


if __name__ == '__main__':
    unittest.main()
