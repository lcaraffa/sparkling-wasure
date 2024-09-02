from setuptools import setup, find_packages

import myproject

setup(
    name=myproject.__name__,
    version=myproject.__version__,
    url='https://github.com/maet3608/minimal-python-project',
    author='Author name',
    author_email='author@gmail.com',
    description='Template and example for a minimal Python project',
    packages=find_packages(),
    install_requires=[],  # e.g. ['numpy >= 1.11.1', 'matplotlib >= 1.5.1']
)

