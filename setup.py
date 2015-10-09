__author__ = 'anna'

from setuptools import setup, find_packages
from sys import version

if version < '2.2.3':
    from distutils.dist import DistributionMetadata

    DistributionMetadata.classifiers = None
    DistributionMetadata.download_url = None


setup(name='mod_sbml',
      description='Utilities for working with SBML models.',
      long_description=open('README.md').read(),
      author='Anna Zhukova',
      author_email='zhutchok@gmail.com',
      version='0.1',
      packages=find_packages(),
      package_data={'mod_sbml': ['data/*.obo', 'data/*.csv', 'data/*.txt']},
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
      download_url='https://github.com/annazhukova/mod_sbml',
      url='https://github.com/annazhukova/mod_sbml',
      install_requires=['openpyxl', 'python-libsbml-experimental', 'pandas', 'matplotlib', 'colorsys']
      )

