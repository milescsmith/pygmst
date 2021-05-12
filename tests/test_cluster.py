import json
import logging
import pkgutil
import unittest
from os.path import exists

from pkg_resources import resource_filename

from pygmst.pygmst import cluster

logging.basicConfig(level=logging.CRITICAL)


class TestClusterFunction(unittest.TestCase):
    def setUp(self):
        self.cluster_answers = json.loads(
            pkgutil.get_data("tests", "cluster_answers.json")
        )
        print("json loaded")
        self.assertTrue(exists(resource_filename("tests", "initial.meta.list.feature")))

    def test_cluster(self):
        print("run cluster")
        bin_num, cutoffs, seq_GC = cluster(
            feature_f=resource_filename("tests", "initial.meta.list.feature"),
            clusters=0,
            min_length=50000,
        )

        self.assertEqual(
            bin_num,
            self.cluster_answers["bin_num"],
            msg=f"function 'cluster' failed: bin_num is incorrect. should be {self.cluster_answers['bin_num']} but was {bin_num}",
        )
        self.assertEqual(
            cutoffs,
            self.cluster_answers["cutoffs"],
            msg=f"function 'cluster' failed: cutoffs is incorrect. should be {self.cluster_answers['cutoffs']} but was {cutoffs}",
        )

        self.assertDictEqual(
            seq_GC,
            self.cluster_answers["seq_GC"],
            msg=f"function 'cluster' failed: seq_GC is incorrect.",
        )


if __name__ == "__main__":
    unittest.main()
