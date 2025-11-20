# from Cython.Build import cythonize
from setuptools import setup

setup(name='DeepHF',
      version='1.0',
      author='Jean-Philippe Fortin',
      author_email="fortin946@gmail.com",
      description=("Python package to run DeepHF algorithm"),
      packages=["deephf"],
      package_data={'deephf': ['models/*.hd5']},
      license="MIT",
      )
