import logging
import tempfile
import unittest
from os.path import exists

from pkg_resources import resource_filename

from pygmst.pygmst import gmst

logging.basicConfig(level=logging.CRITICAL)


class TestGMSTFunction(unittest.TestCase):
    def setUp(self):
        self.fasta = resource_filename("tests", "test.fa")
        self.results = resource_filename("tests", "gmst_test")
        self.faa = resource_filename("tests", "gmst_test.faa")
        self.fnn = resource_filename("tests", "gmst_test.fnn")

        self.assertTrue(exists(self.fasta), msg="Cannot find test.fa")
        self.assertTrue(exists(self.results), msg="Cannot find gmst_test")
        self.assertTrue(exists(self.faa), msg="Cannot find gmst_test.faa")
        self.assertTrue(exists(self.fnn), msg="Cannot find gmst_test.fnn")

    def test_gmst(self):
        print("run gmst")
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = f"{tmpdir}/unittest"
            gmst(
                seqfile=self.fasta, output=prefix, faa=True, fnn=True, strand="direct"
            )  # test how we are actually going to use the module

            with open(f"{prefix}", "r") as res, open(
                f"{prefix}.faa", "r"
            ) as resfaa, open(f"{prefix}.fnn", "r") as resfnn:
                test_results = "".join(
                    res.readlines()[7:]
                )  # the first few lines include a timestamp and the tmpdir, so they will always differ.  skip those lines.
                test_faa = "".join(resfaa.readlines())
                test_fnn = "".join(resfnn.readlines())

            with open(self.results, "r") as exp, open(self.faa, "r") as expfaa, open(
                self.fnn, "r"
            ) as expfnn:
                expected_results = "".join(exp.readlines()[7:])
                expected_faa = "".join(expfaa.readlines())
                expected_fnn = "".join(expfnn.readlines())

        self.assertMultiLineEqual(
            test_results,
            expected_results,
            msg=f"function 'gmst' failed: results did not match expected",
        )

        self.assertMultiLineEqual(
            test_faa,
            expected_faa,
            msg=f" 'gmst' failed: the resulting test.faa file did not match expected",
        )

        self.assertMultiLineEqual(
            test_fnn,
            expected_fnn,
            msg=f"f 'gmst' failed: the resulting test.fnn file did not match expected",
        )


if __name__ == "__main__":
    unittest.main()
