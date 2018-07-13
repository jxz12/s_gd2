from setuptools import setup
from setuptools.extension import Extension

# Third-party modules - we depend on numpy for everything
import numpy
# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


_s_gd2 = Extension(
    name="_s_gd2",
    sources=["layout_wrap.cxx", "layout.cpp"],
    include_dirs = [numpy_include]
)

setup(
    name="s_gd2",
    version="0",
    author="Jonathan Zheng",
    author_email="jxz12@ic.ac.uk",
    description="A package for performing stochastic gradient descent (arXiv:1710.04626) to layout graphs",
    py_modules=['s_gd2'],
    ext_modules=[_s_gd2]
)

