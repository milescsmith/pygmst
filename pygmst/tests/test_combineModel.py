# import unittest
# import pkgutil
# import tempfile
# import json


# from pkg_resources import resource_filename
# from os.path import exists, splitext, basename, abspath
# from subprocess import run
# from .. import pygmst


# class TestTrainFunction(unittest.TestCase):
#     def setUp(self):


#         self.testfasta = resource_filename("pygmst", "tests/test.fa")
#         self.assertTrue(exists(self.testfasta), msg="Cannot find test.fa")


#     def test_cluster(self):

#         with tempfile.TemporaryDirectory() as tmpdir:
#             final_model = combineModel(models, cutoffs)

            

#             self.assertGreaterEqual(
#                 self.diff,
#                 0.99,
#                 msg=f"function 'train' failed: the resulting model is not similar to the expected one",
#             )


# if __name__ == "__main__":
#     unittest.main()
