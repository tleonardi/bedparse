from setuptools import setup

setup(name='bedParse',
      version='0.0.1a',
      description='A simple library and CLI tool to manipulate BED files',
      url='https://github.com/tleonardi/bedParse',
      author='Tommaso Leonardi',
      author_email='tom@tleo.io',
      license='MIT',
      packages=['bedParse'],
      install_requires=['argparse'],
      scripts=['bin/bedParse'],
      zip_safe=False)
