from setuptools import setup, find_packages

# information about the dingo library
version = "0.0.1"
license = ("LGPL3",)
packages = ["src/fluxpy"]
description = "A python library for metabolic networks sampling and analysis"
author = "Haris Zafeiropoulos"
author_email = "haris.zafeiropoulos@kuleuven.be"
name = "fluxpy"


setup(
    name='fluxpy',
    version='0.0.1',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    # packages=packages,  # find_packages()
    package_data={
        'fluxpy': ['data/*'],  # 'illustrations/*.json'
    },
    description='A set of tools to support metabolic modelling analysis related tasks.',
    author='Haris Zafeiropoulos',
    author_email='haris.zafeiropoulos@kuleuven.be',
    url='https://github.com/hariszaf/fluxpy',
    license='MIT',
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'plotly',
        'cobra',
        'mergem'
    ],
)
