#!/usr/bin/env python3

from setuptools import setup

setup(
    name='gpcgrlib',
    version='0.0.0-SNAPSHOT',

    description='common libraries for gpc-gr project',
    long_description='common libraries for gpc-gr project',
    url='',

    author='gpc-gr',
    author_email='gr@gpc.megabank.tohoku.ac.jp',
    license='Apache-2.0',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.6'
    ],
    keywords='',
    
    packages=[
        'gpcgrlib'
    ],
    entry_points="""
        [console_scripts]
        tqsub = gpcgrlib.bin.tqsub:main
    """,

    install_requires=[
        'jinja2'
    ]
)

