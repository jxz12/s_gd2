from setuptools import setup, find_packages
from setuptools.extension import Extension

# Third-party modules - we depend on numpy for everything
import numpy
# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

_layout = Extension(
    name="_layout",
    headers=["./s_gd2/layout.hpp"],
    sources=["./s_gd2/layout.cpp", "./s_gd2/sparse.cpp", "./s_gd2/swig/layout_wrap.cxx"],
    extra_compile_args=["-std=c++11"],
    include_dirs=[numpy_include],
)

setup(
    name="s_gd2",
    version="1.4",
    author="Jonathan Zheng",
    author_email="jxz12@ic.ac.uk",
    url="https://www.github.com/jxz12/s_gd2",
    description="A package for performing stochastic gradient descent (arXiv:1710.04626) to layout graphs",
    install_requires=['numpy', 'pycairo'],
    packages=find_packages(),
    ext_modules=[_layout]
)
