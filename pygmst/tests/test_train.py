import unittest
import pkgutil
import tempfile
from numpy.testing import assert_approx_equal

from pkg_resources import resource_filename
from os.path import exists, splitext, basename, abspath
from subprocess import run
from pygmst import train

import logging

logging.basicConfig(filename="pygmst_test_train.log", filemode="w", level=logging.DEBUG)


class TestTrainFunction(unittest.TestCase):
    def setUp(self):
        # with open(
        #     resource_filename("pygmst", "tests/test_cluster.mod"), "r"
        # ) as test_cluster:
        #     self.expected_model = test_cluster.readlines()
        logging.basicConfig(
            filename="pygmst_test_train.log", filemode="w", level=logging.DEBUG
        )

        self.testfasta = resource_filename("pygmst", "tests/test.fa")
        self.testsequence = resource_filename("pygmst", "tests/test_sequence")
        self.probuild = resource_filename("pygmst", "genemark/probuild")
        self.gmhmmp = resource_filename("pygmst", "genemark/gmhmmp")
        self.par_1 = resource_filename("pygmst", "genemark/par_1.default")
        self.gibbs3 = resource_filename("pygmst", "genemark/Gibbs3")

        self.assertTrue(exists(self.testfasta), msg="Cannot find test.fa")
        self.assertTrue(exists(self.testsequence), msg="Cannot find test_sequence")
        self.assertTrue(exists(self.probuild), msg="Cannot find probuild")
        self.assertTrue(exists(self.gmhmmp), msg="Cannot find gmhmmp")
        self.assertTrue(exists(self.par_1), msg="Cannot find par_1.default")
        self.assertTrue(exists(self.gibbs3), msg="Cannot find Gibbs3")

    def test_cluster(self):

        with tempfile.TemporaryDirectory() as tmpdir:
            test_model, _ = train(
                input_seq=self.testfasta,
                seq=self.testsequence,
                motif=True,
                fixmotif=True,
                order=4,
                order_non=2,
                start_prefix="startseq.",
                gibbs_prefix="itr_",
                prestart=6,
                width=12,
                build_cmd=f"{self.probuild} --par",
                hmm_cmd=f"{self.gmhmmp} -s d",
                par=f"{self.par_1}",
                maxitr=10,
                identity=0.99,
                gibbs3=self.gibbs3,
                tmpdir=tmpdir,
            )

            test_lst = f"{abspath(test_model).split('.')[0]}.lst"
            logging.debug("test_lst")
            command = f"{self.probuild} --par {self.par_1} --compare --source {resource_filename('pygmst', 'tests/actual.lst')} --target {test_lst}"
            logging.debug(command)
            self.diff = float(
                str(run(command.split(), capture_output=True).stdout, "utf-8").strip(
                    "\n"
                )
            )
            logging.debug(self.diff)

            self.assertGreaterEqual(
                self.diff,
                0.99,
                msg=f"function 'train' failed: the resulting model is not similar to the expected one",
            )


if __name__ == "__main__":
    unittest.main()
