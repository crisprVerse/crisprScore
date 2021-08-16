# from Cython.Build import cythonize
from setuptools import setup


setup(name='DeepHF',
      version='1.0',
      author='Jean-Philippe Fortin',
      author_email="fortin946@gmail.com",
      description=("Python package to run DeepHF algorithm"),
      packages=["deephf"],
      package_data={'deephf': ['models/*.hd5']},
      install_requires=['scipy==1.1.0', 'numpy==1.14.0', 
      'h5py==2.9.0', 'tensorflow==1.8.0', 'keras==2.1.6',
      'scikit-learn==0.19.1','biopython==1.71','matplotlib',
      'DotMap','GPyOpt','pandas'],
      license="MIT",
      )
