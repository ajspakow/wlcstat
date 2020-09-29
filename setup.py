#!/usr/bin/env python

import versioneer
from setuptools import setup


def _readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='wlcstat',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Utilities for analyzing particle trajectories',
    long_description=_readme(),
    author='Bruno Beltran',
    author_email='brunobeltran0@gmail.com',
    packages=['wlcstat'],
    license='MIT',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Utilities'
    ],
    keywords='polymer statistics scientific',
    url='https://github.com/ajspakow/wlcstat',
    install_requires=['numpy', 'scipy', 'numba'],
)
