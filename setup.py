from setuptools import setup

setup(name='bedparse',
      description='A simple library and CLI tool to manipulate BED files',
      version="0.1.2",
      url='https://github.com/tleonardi/bedparse',
      author='Tommaso Leonardi',
      author_email='tom@tleo.io',
      license='MIT',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.2',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6'
      ],
      packages=['bedparse'],
      install_requires=['argparse', 'setuptools'],
      python_requires='>=3',
      entry_points={
          'console_scripts': [
              'bedparse = bedparse.bedparse:main'
          ]
      },
      zip_safe=False)
