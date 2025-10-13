from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="deepsea-edna",
    version="0.1.0",
    author="Deep Divers Team",
    author_email="team@deepdivers.org",
    description="A comprehensive pipeline for processing environmental DNA sequences from deep-sea samples",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Deep-Divers/deepsea-edna",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "deepsea=deepsea_edna.cli.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "deepsea_edna": [
            "data/*",
            "module1_qc_asv/primers/*.json",
            "module1_qc_asv/*.yaml",
            "module2_features/*.yaml",
        ],
    },
)