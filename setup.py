from setuptools import setup

setup(
    name='gsc',
    version='2.0',
    description='"gsc" or "graph state compass" is a Python package containing tools for mapping and depicting the '
    'local Clifford equivalence classes of quantum graph states.',
    url='https://github.com/sammorley-short/gsc',
    author='Sam Morley-Short',
    license='GNU General Public License v3.0',
    packages=['gsc'],
    install_requires=[
        'abp',
        'numpy',
        'sympy',
        'matplotlib==2.2.3',
        'networkx'
    ],
    zip_safe=False
)
