#!/usr/bin/env python

from setuptools import setup, find_packages
from procread import __version__


setup(
    name='procread',
    version=__version__,
    description='Read-to-variant pipeline tool for DNA-seq analyses',
    packages=find_packages(),
    author='Daichi Narushima',
    author_email='d.narsil@gmail.com',
    url='https://github.com/dceoy/procread',
    include_package_data=True,
    install_requires=[
        'biopython',
        'docopt',
        'pyyaml'
    ],
    entry_points={
        'console_scripts': ['procread=procread.main:main']
    },
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Environment :: Console',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers'
    ]
)
