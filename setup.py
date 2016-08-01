#!/usr/bin/env python

from distutils.core import setup

setup(
    name='srst2',
    version='0.2.0',
    author='Kathryn Holt',
    author_email='drkatholt@gmail.com',
    packages=['srst2'],
    scripts=['scripts/getmlst.py', 'scripts/slurm_srst2.py'],
    entry_points={
        'console_scripts': ['srst2 = srst2.srst2:main']
    },
    package_dir = {'srst2': 'scripts'},
    package_data={'srst2': ['data/resistance.*']},
    url='http://katholt.github.io/srst2/',
    license='LICENSE.txt',
    description='Short Read Sequence Typing for Bacterial Pathogens',
    long_description=('This program is designed to take Illumina'
                      'sequence data, a MLST database and/or a database'
                      'of gene sequences (e.g. resistance genes, virulence'
                      'genes, etc) and report the presence of STs and/or'
                      'reference genes.'),
    install_requires=[
        # Although we depend on scipy, which depends on numpy, we don't
        # specify the dependencies here because they don't play well with
        # any Python installing system, such as pip or easy_install.
        # So we assume the user has already installed the dependencies
        # themselves.
        #"numpy >= 1.7.1",
        #"scipy >= 0.12.0",
    ],
)
