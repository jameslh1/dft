from setuptools import setup, find_packages

print( 'DFT module \n\n setup.py script ')

setup(name='dft',
      version='1.0',
      description = 'DFT module',
      author = 'James Lloyd-Hughes',
      author_email = 'J.Lloyd-Hughes@warwick.ac.uk',
      install_requires=[
        'numpy',
        'matplotlib'
      ],
      packages=find_packages(),
      url='https://github.com/jameslh1/dft/',
      classifiers=[
        'Programming Language :: Python :: 2 or 3',
        'Operating System :: OS Independent',
      ],
      python_requires='>=3.6',      
      )
