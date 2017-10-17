#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Setup script for SHTOOLS."""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

# the setuptools import dummy patches the distutil commands such that
# python setup.py develop works
import setuptools  # NOQA

from numpy.distutils.core import setup

# convert markdown README.md to restructured text .rst for pypi
# pandoc can be installed with
# conda install -c conda-forge pandoc pypandoc
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except(IOError, ImportError):
    print('no pandoc installed. Careful, pypi description will not be '
          'formatted correctly.')
    long_description = open('README.md').read()


VERSION = 1.0
print('obspy 3d visualization package version: {}'.format(VERSION))


CLASSIFIERS = [
    'Development Status :: 5 - Production/Stable',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',  # NOQA
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics'
]


KEYWORDS = ['Obspy', 'VTK', 'Visualization', 'Mayavi', '3D Ray Paths']


INSTALL_REQUIRES = [
    'numpy (>=1.0.0)',
    'scipy',
    'matplotlib (>=1.5)']


metadata = dict(
    name='vtkobs',
    version=VERSION,
    description='3d vtk based routines for obspy',
    long_description=long_description,
    url='https://github.com/MMesch/vtkobs',
    download_url='https://github.com/MMesch/cmap_builder/zipball/master',
    author='MMesch',
    author_email="MMesch@users.noreply.github.com",
    license='GPL v3',
    keywords=KEYWORDS,
    requires=INSTALL_REQUIRES,
    platforms='OS Independent',
    packages=['vtkobs'],
    package_dir={'vtkobs': 'vtkobs'},
    classifiers=CLASSIFIERS
)


setup(**metadata)
