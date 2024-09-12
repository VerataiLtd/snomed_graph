from setuptools import setup, find_packages

setup(
    name="snomed_graph",
    version="0.2.0",
    author="Will Hardman",
    author_email="will.hardman@veratai.co.uk",
    description="A simple Python library for working with SNOMED CT",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/VerataiLtd/snomed_graph",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Apache License 2.0",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "pandas",
        "tqdm",
        "networkx>=3.0",
    ],
)