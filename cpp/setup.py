# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# sgd extension module
_s_gd2 = Extension("_s_gd2",
                   ["layout.i", "layout.cpp"],
                   include_dirs = [numpy_include],
                   )

# sgd setup
setup(  name        = "stochastic gradient descent for graph drawing",
        description = "sgd() takes an array of 2d positions, condensed distance matrices d and w, and a step size schedule eta, and minimizes stress.",
        author      = "Jonathan Zheng",
        version     = "1.0",
        ext_modules = [_s_gd2]
        )
