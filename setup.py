import os
from distutils.core import setup
import sys

if sys.version_info < (3,4,1):
    sys.exit('Sorry, Python < 3.4.1 is not supported')

setup(
    name='mod_sbml',
    packages=['mod_sbml'],
    package_data={'mod_sbml': [os.path.join('data', '*.obo'), os.path.join('data', '*.csv'),
                               os.path.join('data', '*.txt'), os.path.join('*', '*.py'),
                               os.path.join('*', '*', '*.py'),
                               os.path.join('..', 'README.md')]},
    long_description=open('README.md').read(),
    include_package_data=True,
    platform=['MacOS', 'Linux', 'Windows'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='0.2.3',
    description='Utilities for working with metabolic models in SBML format.',
    author='Anna Zhukova',
    author_email='zhutchok@gmail.com',
    url='https://github.com/annazhukova/mod_sbml',
    download_url='https://github.com/annazhukova/mod_sbml/archive/0.2.3.zip',
    keywords=['SBML', 'metabolic model', 'utility'],
    install_requires=['openpyxl', 'python-libsbml-experimental', 'pandas', 'matplotlib', 'pyparsing', 'natsort']
)
