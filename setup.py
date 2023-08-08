from setuptools import setup, find_packages
import numpy as np

# Package metadata
NAME = 'wlcstat'
VERSION = '0.1.0'
AUTHOR = 'Andrew Spakowitz'
AUTHOR_EMAIL = 'ajspakow@stanford.edu'
DESCRIPTION = 'WLC Statistics Package'
URL = 'https://github.com/ajspakow/wlcstat.git'
KEYWORDS = 'thermodynamics, epigenetics, nucleosomes, euchromatin, polymer'
CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science Research',
    'Programming Language :: Python :: 3.9',
]

# Package dependencies
INSTALL_REQUIRES = [
    "matplotlib~=3.2.1",
    "numba~=0.44.1",
    "numpy~=1.18.4",
    "pandas~=1.0.3",
    "pathos~=0.3.1",
    "rotations~=0.0.2",
    "scipy~=1.4.1",
    "setuptools~=40.8.0",
    "sympy~=1.4"
]

# Read the long description from README.md
with open('README.md', 'r', encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

# Setup configuration
setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url=URL,
    keywords=KEYWORDS,
    classifiers=CLASSIFIERS,
    packages=find_packages(include=["wlcstat", "wlcstat/util"]),
    include_package_data=True,
    include_dirs=[np.get_include(), "wlcstat", "wlcstat/util"],
    install_requires=INSTALL_REQUIRES,
)
