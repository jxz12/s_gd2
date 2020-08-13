from setuptools import setup, find_packages
from setuptools.extension import Extension

import sys
import os

# Third-party modules - we depend on numpy for everything
class get_numpy_include(object):
    """Returns Numpy's include path with lazy import.
    """
    def __str__(self):
        import numpy
        # Obtain the numpy include directory.  This logic works across numpy versions.
        try:
            return numpy.get_include()
        except AttributeError:
            return numpy.get_numpy_include()

python_version = (int(sys.version[0]), int(sys.version[2]))
numpy_version = 'numpy>=1.16'
if python_version < (3, 5):
    numpy_version += ',<1.17'
elif python_version < (3, 6):
    numpy_version += ',<1.18'

setup_requires = [
    numpy_version,
    'Cython',
]

install_requires = [
    numpy_version,
]

test_requires = [
    'nose2',
    'scipy',
]

_layout = Extension(
    name="_layout",
    headers=["./s_gd2/layout.hpp", "./s_gd2/randomkit.h"],
    sources=["./s_gd2/layout.cpp", "./s_gd2/sparse.cpp", "./s_gd2/swig/layout_wrap.cxx", "./s_gd2/randomkit.c"],
    include_dirs=[get_numpy_include()],
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
    setup_requires=setup_requires,
    install_requires=install_requires,
    extras_require={'test':test_requires},
    test_suite='nose2.collector.collector',
    packages=find_packages(),
    ext_modules=[_layout]
)
