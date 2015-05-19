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
      author_email='anna.zhukova@ibgc.cnrs.fr',
      version='0.1',
      packages=find_packages(),
      package_data={'': ['data/', 'data/*.obo', 'data/*.xlsx']},
      include_package_data=True,
      license='LICENSE',
      platform=['MacOS', 'Linux', 'Windows'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'License :: GNU GENERAL PUBLIC LICENSE',
          'Topic :: Systems Biology',
          'Topic :: Software Development',
      ], requires=['openpyxl', 'python-libsbml-experimental', 'numpy', 'matplotlib']
)