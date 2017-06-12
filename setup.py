from setuptools import setup

setup(name='bedParse',
      description='A simple library and CLI tool to manipulate BED files',
      version="0.0.2",
      url='https://github.com/tleonardi/bedParse',
      author='Tommaso Leonardi',
      author_email='tom@tleo.io',
      license='MIT',
      packages=['bedparse'],
      install_requires=['argparse'],
      entry_points={
          'console_scripts': [
              'bedparse = bedparse.bedparse:main'
          ]
      },
      zip_safe=False)
