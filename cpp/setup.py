from setuptools import setup, find_packages
from setuptools.extension import Extension

import os

# Third-party modules - we depend on numpy for everything
import numpy
# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

install_requires = [
    'numpy>=1.16',
    'matplotlib',
]

test_requires = [
    'nose2',
    'scipy',
]

_layout = Extension(
    name="_layout",
    headers=["./s_gd2/layout.hpp"],
    sources=["./s_gd2/layout.cpp", "./s_gd2/sparse.cpp", "./s_gd2/swig/layout_wrap.cxx"],
    extra_compile_args=["-std=c++11"],
    include_dirs=[numpy_include],
)

version_py = os.path.join(os.path.dirname(__file__), 's_gd2', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"', '').strip()

setup(
    name="s_gd2",
    version=version,
    author="Jonathan Zheng",
    author_email="jxz12@ic.ac.uk",
    url="https://www.github.com/jxz12/s_gd2",
    description="A package for performing stochastic gradient descent (arXiv:1710.04626) to layout graphs",
    install_requires=install_requires,
    extras_require={'test':test_requires},
    test_suite='nose2.collector.collector',
    packages=find_packages(),
    ext_modules=[_layout]
)
