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
    Extension('ultraplex._align', sources=['ultraplex/_align.pyx']),
    Extension('ultraplex.qualtrim', sources=['ultraplex/qualtrim.pyx']),
]


def no_cythonize(extensions, **_ignore):
    """
    Change file extensions from .pyx to .c or .cpp.
    Copied from Cython documentation

    This is run only if a "PKG-INFO" file is present
    """
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in ('.pyx', '.py'):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources


# This call to setup() does all the work
setup(
    name="ultraplex",
    version="1.0.3",
    description="fastq demultiplexer",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/ulelab/ultraplex.git",
    author="Oscar Wilkins",
    author_email="oscar.wilkins@crick.ac.uk",
    license="MIT",
    ext_modules=extensions,
   	    package_dir={'': 'ultraplex'},
    packages=find_packages('ultraplex'),
    install_requires=[
        'dnaio~=0.5.0',
        'xopen~=1.0.0',
        "dataclasses>=0.8; python_version<'3.7'",
    ],
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    include_package_data=True,
    #install_requires=["cutadapt==2.10"],

     entry_points={
        "console_scripts": [
            "ultraplex = ultraplex.__main__:main",
        ]
    },
)
