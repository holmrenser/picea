#!/usr/bin/python
from setuptools import setup
from sphinx.setup_command import BuildDoc
import picea

name = 'picea'
version = picea.__version__
author = picea.__author__


def main():
    setup(
        name=name,
        packages=[name],
        author=author,
        author_email='rens.holmer@wur.nl',
        description='A lightweight python library for working with \
            trees and sequence collections',
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
