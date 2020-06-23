import unittest
import pkgutil
import tempfile
import json
from pygmst import combineModel

from pkg_resources import resource_filename
from os.path import exists, splitext, basename, abspath
from subprocess import run


class TestTrainFunction(unittest.TestCase):
    def setUp(self):
        self.model_1 = resource_filename("pygmst", "tests/test_model_1.mod")
        self.model_2 = resource_filename("pygmst", "tests/test_model_2.mod")
        self.expected_combinedModel = resource_filename(
            "pygmst", "tests/expected_combinedModel.mod"
        )

        self.assertTrue(exists(self.model_1), msg="Cannot find test_model_1.mod")
        self.assertTrue(exists(self.model_2), msg="Cannot find test_model_2.mod")
        self.assertTrue(
            exists(self.expected_combinedModel),
            msg="Cannot find expected_combinedModel.mod",
        )

    def test_cluster(self):

        with tempfile.TemporaryDirectory() as tmpdir:
            final_model = combineModel(
                mod=[self.model_1, self.model_2], cut_offs=[30, 70]
            )

            with open(final_model, "r") as fm, open(
                self.expected_combinedModel, "r"
            ) as em:
                final_model = "".join(fm.readlines())
                expected_model = "".join(em.readlines())

            self.assertMultiLineEqual(
                final_model,
                expected_model,
                msg=f"function 'combineModel' failed: the resulting concatination did not match the expected one",
            )


if __name__ == "__main__":
    unittest.main()
