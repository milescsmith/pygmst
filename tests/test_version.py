import unittest
from pkg_resources import resource_filename
from os.path import exists
import logging
import tempfile

from pygmst.pygmst import main

logging.basicConfig(level=logging.CRITICAL)


class TestVersion(unittest.TestCase):
    
    def test_version(self):
        print(__version__)

if __name__ == "__main__":
    unittest.main()
