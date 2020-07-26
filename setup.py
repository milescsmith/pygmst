#!/usr/bin/env python
import os.path as osp
import sys
from pathlib import Path

from setuptools import setup, find_packages

if sys.version_info < (3,):
    sys.exit("pygmst requires Python >= 3.7")

with open(osp.join("src", "pygmst", "__about__.py")) as f:
    exec(f.read())


setup(
    name="pygmst",
    version=__version__,
    description="Python reimplementation of GMST: Identification of protein coding regions in RNA transcripts",
    long_description=Path("README.rst").read_text("utf-8"),
    url="https://github.com/milescsmith/pygmst",
    author=__author__,
    author_email=__email__,
    license="MPL2",
    python_requires=">=3.7",
    install_requires=[
        line.strip() for line in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    # extras_require=dict(doc=["sphinx", "sphinx_rtd_theme", "sphinx_autodoc_typehints"]),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "pygmst = pygmst.pygmst:main"
        ]
    },
    packages=find_packages(where="src"),
    package_dir={"pygmst": "src/pygmst"},
    package_data={"": ["pygmst/genemark/*.*",]},
)
