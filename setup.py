import pathlib
from setuptools import setup, Extension, find_packages

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()


### Cython ###

"""
Ultraplex uses two cython functions from cutadapt 2.10.
It's considered best practice to distribute pre-cythonised files
rather than making each user use cython to cythonise the .pyx
files.

The setup therefore works in the following way:

1. First, we give the file of the cython scripts.

2. The installer checks whether there is a "'PKG-INFO file". If so, 
this means that the .pyx files are actually .c/.cpp files in disguise,
so the .c/.cpp files are renamed to .c/.cpp.

3. If 

"""

extensions = [
    Extension('_align', sources=['ultraplex/_align.c']),
    Extension('qualtrim', sources=['ultraplex/qualtrim.c']),
]

# This call to setup() does all the work
setup(
    name="ultraplex",
    version="1.0.4.2",
    description="fastq demultiplexer",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/ulelab/ultraplex.git",
    author="Oscar Wilkins",
    author_email="oscar.wilkins@crick.ac.uk",
    license="MIT",
    ext_modules=extensions,
   	#package_dir={'': 'ultraplex'},
    #packages=find_packages(''),
    install_requires=[
        'dnaio~=0.5.0',
        'xopen~=1.0.0',
        "dataclasses>=0.8; python_version<'3.9'",
    ],
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    include_package_data=True,
     entry_points={
        "console_scripts": [
            "ultraplex_test = ultraplex.__main__:main",
        ]
    },
)
