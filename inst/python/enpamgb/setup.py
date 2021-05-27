#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['biopython>=1.74', 'joblib>=0.13.2', 'numpy>=1.16.4', 'pandas>=0.24.2',
                'scikit-learn==0.21.2', 'scipy>=1.3.0', 'tensorflow>=1.14.0',]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Peter DeWeirdt",
    author_email='petedeweirdt@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="On target modeling for CRISPR guides",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='sgrna_modeler',
    name='sgrna_modeler',
    packages=find_packages(include=['sgrna_modeler']),
    package_data={'sgrna_modeler': ['data/datasets/*.csv.zip', 'data/features/*.csv.zip',
                                    'data/saved_models/*.h5', 'data/saved_models/*.joblib']},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/peterdeweirdt/sgrna_modeler',
    version='0.1.0',
    zip_safe=False,
)
