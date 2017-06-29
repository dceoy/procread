#!/usr/bin/env python

from setuptools import setup, find_packages
from procread import __version__


setup(
    name='procread',
    version=__version__,
    packages=find_packages(),
    install_requires=[
        'docopt',
        'luigi'
    ],
    entry_points={
        'console_scripts': ['procread=procread.main:main']
    },
)
