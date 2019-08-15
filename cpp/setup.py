from setuptools import setup
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
    sources=["layout.cpp", "sparse.cpp", "layout_wrap.cxx"],
    extra_compile_args=["-std=c++11"],
    include_dirs=[numpy_include]
)

setup(
    name="s_gd2",
    version="0.20",
    author="Jonathan Zheng",
    author_email="jxz12@ic.ac.uk",
    url="https://www.github.com/jxz12/s_gd2",
    description="A package for performing stochastic gradient descent (arXiv:1710.04626) to layout graphs",
    install_requires=['numpy'],
    py_modules=['s_gd2', 'layout'],
    ext_modules=[_layout]
)

