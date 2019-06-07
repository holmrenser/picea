#!/usr/bin/python
__author__ = 'rensholmer'
__created__ = '07/06/2019'

import sys
from setuptools import setup

def main():
    setup(name='picea',
        packages=['picea'],
        author='rens holmer',
        author_email='rens.holmer@wur.nl',
        description='A lightweight python library for working (phylogenetic) trees',
        version='0.0.1',
        url='https://github.com/holmrenser/picea'
        )

if __name__ == '__main__':
    main()