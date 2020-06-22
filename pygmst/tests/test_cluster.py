import unittest
import pkgutil
from pkg_resources import resource_filename
from os.path import exists
import json
from .. import pygmst


class TestClusterFunction(unittest.TestCase):
    def setUp(self):
        self.cluster_answers = json.loads(
            pkgutil.get_data("pygmst", "tests/cluster_answers.json")
        )
        self.assertTrue(exists(resource_filename("pygmst", "tests/initial.meta.list.feature")))

    def test_cluster(self):
        bin_num, cutoffs, seq_GC = pygmst.cluster(
            feature_f=resource_filename("pygmst", "tests/initial.meta.list.feature"), clusters=0, min_length=50000
        )
        self.assertEqual(
            first=bin_num,
            second=self.cluster_answers["bin_num"],
            msg=f"function 'cluster' failed: bin_num is incorrect. should be {self.cluster_answers['bin_num']} but was {bin_num}",
        )
        self.assertEqual(
            first=cutoffs,
            second=self.cluster_answers["cutoffs"],
            msg=f"function 'cluster' failed: bin_num is incorrect. should be {self.cluster_answers['cutoffs']} but was {cutoffs}",
        )
        self.assertEqual(
            first=seq_GC,
            second=self.cluster_answers["seq_GC"],
            msg=f"function 'cluster' failed: bin_num is incorrect.",
        )

if __name__ == "__main__":
    unittest.main()