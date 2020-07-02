import unittest
from pkg_resources import resource_filename
from os.path import exists
import logging
import tempfile

try:
    from importlib import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as metadata

logging.basicConfig(level=logging.CRITICAL)


class TestVersion(unittest.TestCase):
    
    def test_version(self):
        print(metadata.version('pygmst'))

if __name__ == "__main__":
    unittest.main()
