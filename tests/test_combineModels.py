import unittest
import tempfile
from pygmst.pygmst import combineModels

from pkg_resources import resource_filename
from os.path import exists, basename, abspath


class TestCombineFunction(unittest.TestCase):
    def setUp(self):
        self.model_1 = resource_filename("tests", "test_model_1.mod")
        self.model_2 = resource_filename("tests", "test_model_2.mod")
        self.expected_combinedModels = resource_filename(
            "tests", "expected_combinedModels.mod"
        )

        self.assertTrue(exists(self.model_1), msg="Cannot find test_model_1.mod")
        self.assertTrue(exists(self.model_2), msg="Cannot find test_model_2.mod")
        self.assertTrue(
            exists(self.expected_combinedModels),
            msg="Cannot find expected_combinedModels.mod",
        )

    def test_cluster(self):

        # with tempfile.TemporaryDirectory() as tmpdir:
        final_model = combineModels(
            mod=[self.model_1, self.model_2], cut_offs=[30, 70], #tmpdir=tmpdir
        )

        with open(final_model, "r") as fm, open(
            self.expected_combinedModels, "r"
        ) as em:
            final_model = "".join(fm.readlines())
            expected_model = "".join(em.readlines())

        self.assertMultiLineEqual(
            final_model,
            expected_model,
            msg=f"function 'combineModels' failed: the resulting concatination did not match the expected one",
        )


if __name__ == "__main__":
    unittest.main()
