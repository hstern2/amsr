import os
from glob import glob
from setuptools import setup

exec(open("amsr/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="amsr",
    version=__version__,
    description="Another molecular string representation",
    author="Harry Stern",
    url="https://hstern2.github.io/amsr/",
    license="MIT",
    packages=["amsr"],
    install_requires=["rdkit-pypi", "networkx", "anytree", "Levenshtein"],
    test_suite="tests",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
