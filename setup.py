import pathlib
from setuptools import setup, Extension, find_packages
import os


# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

## Compile Cython
if os.path.isfile("ultraplex/_align_new.c"):
    USE_CYTHON = False
else:
    USE_CYTHON = True

ext = ".pyx" if USE_CYTHON else ".c"

extensions = [
    Extension("_align_new", sources=["ultraplex/_align_new" + ext]),
    Extension("qualtrim_new", sources=["ultraplex/qualtrim_new" + ext]),
]

if USE_CYTHON:
    from Cython.Build import cythonize

    extensions = cythonize(extensions, compiler_directives={"language_level": "3"})

# This call to setup() does all the work
setup(
    name="ultraplex",
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    python_requires=">=3.7",
    description="fastq demultiplexer",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/ulelab/ultraplex.git",
    author="Oscar Wilkins",
    author_email="oscar.wilkins@crick.ac.uk",
    license="MIT",
    ext_modules=extensions,
    # package_dir={'': 'ultraplex'},
    # packages=find_packages(''),
    packages=find_packages(),
    install_requires=["dnaio>=0.5.0", "xopen~=1.0.0", "multiprocess~=0.70.15"],
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
    ],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "ultraplex = ultraplex.__main__:main",  # main() of ultraplex/__main__.py
        ]
    },
)
