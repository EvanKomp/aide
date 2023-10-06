from setuptools import setup, find_packages
from aide import __version__

setup(
    name="aide",
    version=__version__,
    description="AI-assisted Directed Evolution Unifying Framework (AIDE)",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourgithubusername/aide",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scikit-learn",
        "tensorflow",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    python_requires=">=3.7",
    # now the additional development dependancies
    extras_require={
        "dev": [
            "flake8",
            "jupyter",
            "pytest",
            "sphinx",
        ]
    },
)
