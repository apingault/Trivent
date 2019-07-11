#!/usr/bin/env python

from setuptools import setup
from setuptools import find_packages

install_deps = ['marlinpy_sdhcal']

setup(
    name='install_marlinpy_sdhcal',
    description='dummy file to install marlinpy_sdhcal',
    install_requires=install_deps,
    python_requires=">=2.7",
    dependency_links=['git+ssh://git@gitlab.cern.ch:7999/apingaul/marlinpy_sdhcal.git#egg=marlinpy_sdhcal'],
    platforms="Any")
