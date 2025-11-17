from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension("Bakalarka",
              ["Bakalarka.pyx"],
              include_dirs=[np.get_include()],
              
              ),
              

]

setup(
    ext_modules=cythonize(extensions),
)