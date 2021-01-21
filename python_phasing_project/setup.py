from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Sandia tools',
    ext_modules=cythonize("sandia_stats.pyx"),
    zip_safe=False,
)
