from setuptools import setup, find_packages

setup(
    name='fluxpy',
    version='0.1',
    packages=find_packages(),
    package_data={
        'fluxpy': ['illustrations/*.json'],
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
