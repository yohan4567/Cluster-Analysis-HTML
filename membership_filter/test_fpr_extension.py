import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from membership_filter import ClusterConfig, cluster_cfg
from fpr_extension import filter_fpr, build_catalog


class FprTests(unittest.TestCase):
    def setUp(self):
        self.cluster = ClusterConfig("NGC 5139", "Omega Cen", 201.6968, -47.4796,
                                     5.2, 2.37, 48.389, 48.389/60)
        self.cfg = cluster_cfg(self.cluster.name)
        self.raw = pd.DataFrame([dict(source_id="6083408411248100737", ref_epoch=2017.5,
            ra=201.6968, dec=-47.4796, ra_error=1., dec_error=1., parallax=.19,
            parallax_error=.1, pmra=-3.25, pmra_error=.2, pmdec=-6.746,
            pmdec_error=.2, pmra_pmdec_corr=.5, phot_g_mean_mag=19., ruwe=1.)])

    def test_missing_colour_is_neutral_and_epoch_is_preserved(self):
        result = filter_fpr(self.raw, self.cluster, self.cfg)
        self.assertEqual(len(result), 1)
        self.assertTrue(result.is_member.iloc[0])
        self.assertTrue(result.p_cmd.isna().all())
        self.assertTrue(result.bp_rp.isna().all())
        self.assertEqual(result.catalog_ref_epoch.iloc[0], 2017.5)
        self.assertEqual(result.ref_epoch.iloc[0], 2016.)
        self.assertNotEqual(result.ra.iloc[0], result.catalog_ra.iloc[0])
        self.assertEqual(result.source_id.iloc[0], "6083408411248100737")

    def test_invalid_quality_and_foreground_rejected(self):
        for col, value in [("pmra_error", np.nan), ("pmdec_error", 0),
                           ("phot_g_mean_mag", 21),
                           ("parallax", 10), ("pmra_pmdec_corr", 1.1)]:
            with self.subTest(col=col):
                raw = self.raw.copy()
                raw[col] = value
                self.assertEqual(len(filter_fpr(raw, self.cluster, self.cfg)), 0)

    def test_far_motion_not_member(self):
        raw = self.raw.copy()
        raw["pmra"] = 20
        result = filter_fpr(raw, self.cluster, self.cfg)
        self.assertFalse(result.is_member.iloc[0])
        self.assertEqual(result.membership_prob.iloc[0], 0)

    def test_covariance_changes_score(self):
        raw = self.raw.copy()
        raw["pmra"] += .5
        raw["pmdec"] += .5
        positive = filter_fpr(raw, self.cluster, self.cfg).p_pm.iloc[0]
        raw["pmra_pmdec_corr"] = -.5
        negative = filter_fpr(raw, self.cluster, self.cfg).p_pm.iloc[0]
        self.assertGreater(positive, negative)

    def test_output_preserves_baseline_and_does_not_duplicate_on_rerun(self):
        with tempfile.TemporaryDirectory() as temp:
            baseline = Path(temp)/"baseline.csv"
            original = self.raw.copy()
            original["source_id"] = "6083408411248100736"
            original["membership_prob"] = .82
            original["is_member"] = True
            original.to_csv(baseline, index=False)
            before = baseline.read_bytes()
            with patch("membership_filter.read_targets", return_value=[self.cluster]), \
                 patch("fpr_extension.download_fpr", return_value=self.raw):
                for _ in range(2):
                    result = build_catalog(baseline, Path(temp)/"out", "unused")
                    self.assertEqual(len(result), 2)
                    self.assertEqual(result.membership_prob.iloc[0], .82)
                    self.assertEqual(result.source_key.nunique(), 2)
            self.assertEqual(before, baseline.read_bytes())


if __name__ == "__main__":
    unittest.main()
