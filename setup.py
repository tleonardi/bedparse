from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='bedparse',
      description='A simple library and CLI tool to manipulate BED files',
      long_description=long_description,
      long_description_content_type="text/markdown",
      version="0.2.2",
      url='https://github.com/tleonardi/bedparse',
      author='Tommaso Leonardi',
      author_email='tom@tleo.io',
      license='MIT',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6'
      ],
      packages=['bedparse'],
      install_requires=['argparse', 'setuptools'],
      python_requires='>=3.4',
      entry_points={
          'console_scripts': [
              'bedparse = bedparse.bedparse:main'
          ]
      },
      zip_safe=False)
