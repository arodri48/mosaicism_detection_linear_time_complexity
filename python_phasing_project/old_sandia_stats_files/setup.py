from setuptools import setup
from Cython.Build import cythonize
import numpy
setup(
    name='Sandia tools',
    ext_modules=cythonize("sandia_stats.pyx"),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)
