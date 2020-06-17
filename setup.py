#!/usr/bin/env python
import sys

if sys.version_info < (3,):
    sys.exit("pygmst requires Python >= 3.7")
from pathlib import Path
from setuptools import setup, find_packages

try:
    from pygmst import __author__, __email__
except ImportError:  # Deps not yet installed
    __author__ = "Miles Smith"
    __email__ = "miles-smith@omrf.org"

setup(
    name="pygmst",
    version="0.1.0",
    description="Python reimplementation of GMST: Identification of protein coding regions in RNA transcripts",
    long_description=Path("README.rst").read_text("utf-8"),
    url="https://github.com/milescsmith/pygmst",
    author=__author__,
    author_email=__email__,
    license="MPL2",
    python_requires=">=3.7",
    install_requires=[
        l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    # extras_require=dict(doc=["sphinx", "sphinx_rtd_theme", "sphinx_autodoc_typehints"]),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "gmst = gmst.gmst:main",
        ]
    },
    packages=find_packages(),
    package_dir={"gmst": "gmst"},
    package_data={
        "": [
            "gmst/models/*.*",
            "gmst/utilities/*.*",
        ]
    },
)
