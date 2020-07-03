#!/usr/bin/python
from setuptools import setup
from sphinx.setup_command import BuildDoc
import picea

name = 'picea'
version = picea.__version__
author = picea.__author__

with open('./README.md') as filehandle:
    long_description = filehandle.read()


def main():
    setup(
        name=name,
        packages=[name],
        author=author,
        author_email='rens.holmer@wur.nl',
        description='A lightweight python library for working with \
            trees and sequence collections',
        # long_description=long_description,
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
        ],
        license='MIT',
        platforms=['Windows', 'Windows Cygwin', 'GNU/Linux', 'MacOS'],
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Natural Language :: English',
            'Operating System :: MacOS',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Utilities'
        ]
    )


if __name__ == '__main__':
    main()
