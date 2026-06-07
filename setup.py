# setup.py
from setuptools import setup, find_packages

setup(
    name="brinkgrad",
    version="1.0.0",
    description="Adjoint-based topology optimisation of coupled Brinkman-convection-diffusion systems in FEniCSx",
    long_description=open("README.md", "r", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Nisong Monyimba, Vincent Pizziconi, Aurel Coza",
    author_email="nmonyimb@asu.edu",
    url="https://github.com/NisongMonyimba/brinkgrad",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.23.2",
        "scipy>=1.11.4",
        "matplotlib>=3.8.3",
    ],
    extras_require={
        "fenics": [
            "fenics-dolfinx>=0.7.3",
            "petsc4py>=3.19.2",
            "mpi4py>=3.1.5",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.9",
)
