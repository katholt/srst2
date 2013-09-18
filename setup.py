#!/usr/bin/env python

from distutils.core import setup

setup(
    name='srst2',
    version='0.1.0',
    author='Kathryn Holt',
    author_email='kholt@unimelb.edu.au',
    packages=['srst2'],
    scripts=['srst2/getmlst.py', 'srst2/scores_vs_expected.py', 'srst2/slurm_srst2.py'],
    entry_points={
        'console_scripts': ['srst2 = srst2.srst2:main']
    },
    package_data={'srst2': ['srst2/data/resistance.*']},
    url='https://github.com/katholt/srst',
    license='LICENSE.txt',
    description='Short Read Sequence Typing (for Bacterial Pathogens)',
    long_description=('This program is designed to take Illumina'
                      'sequence data, a MLST database and/or a database'
                      'of gene sequences (e.g. resistance genes, virulence'
                      'genes, etc) and report the presence of STs and/or'
                      'reference genes.'),
    install_requires=[
        "scipy >= 0.12.0",
    ],
)
