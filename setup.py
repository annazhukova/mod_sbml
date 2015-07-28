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
      package_data={'mod_sbml': ['data/*.obo', 'data/*.xlsx', 'data/*.txt']},
      include_package_data=True,
      license='LICENSE',
      platform=['MacOS', 'Linux', 'Windows'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'License :: CeCILL',
          'Topic :: Systems Biology',
          'Topic :: Software Development',
      ],
      download_url='https://github.com/annazhukova/mod_sbml',
      install_requires=['openpyxl', 'python-libsbml-experimental', 'numpy', 'matplotlib']
)
