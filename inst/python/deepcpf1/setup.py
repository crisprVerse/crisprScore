# from Cython.Build import cythonize
from setuptools import setup


setup(name='DeepCpf1',
      version='1.0',
      author='Jean-Philippe Fortin',
      author_email="fortin946@gmail.com",
      description=("Python package to run DeepCpf1 algorithm"),
      packages=["deepcpf1"],
      package_data={'deepcpf1': ['weights/*.h5']},
      license="MIT",
      )
