#!/usr/bin/python
from setuptools import setup
from sphinx.setup_command import BuildDoc

__author__ = 'Rens Holmer'
__created__ = '07/06/2019'

name = 'picea'
version = '0.0.4'


def main():
    setup(
        name=name,
        packages=[name],
        author=__author__,
        author_email='rens.holmer@wur.nl',
        description='A lightweight python library for working (phylogenetic) \
            trees',
        version=version,
        url='https://github.com/holmrenser/picea',
        cmdclass={'build_sphinx': BuildDoc},
        command_options={
            'build_sphinx': {
                'project': ('setup.py', name),
                'version': ('setup.py', version),
                'source_dir': ('setup.py', 'docs/source'),
                'build_dir': ('setup.py', 'docs/build')
            }
        },
        install_requires=[
            'sphinx',
            'numpy',
            'matplotlib'
        ]
    )


if __name__ == '__main__':
    main()
