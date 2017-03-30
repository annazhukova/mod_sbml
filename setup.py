import os
from distutils.core import setup

setup(
    name='mod_sbml',
    packages=['mod_sbml'],
    package_data={'mod_sbml': [os.path.join('data', '*.obo'), os.path.join('data', '*.csv'),
                               os.path.join('data', '*.txt')]},
    include_package_data=True,
    platform=['MacOS', 'Linux', 'Windows'],
    classifiers=[
      'Development Status :: 4 - Beta',
      'Environment :: Console',
      'Intended Audience :: Developers',
      'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Topic :: Software Development',
      'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    version='0.2',
    description='Utilities for working with metabolic models in SBML format.',
    author='Anna Zhukova',
    author_email='zhutchok@gmail.com',
    url='https://github.com/annazhukova/mod_sbml',
    download_url='https://github.com/annazhukova/mod_sbml/archive/v.0.2.zip',
    keywords=['SBML', 'metabolic model', 'utility'],
    install_requires=['openpyxl', 'python-libsbml-experimental', 'pandas', 'matplotlib', 'pyparsing', 'natsort']
)
