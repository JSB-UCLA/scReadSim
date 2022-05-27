import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scReadSim",
    version="0.1.0",
    author="Guanao Yan",
    author_email="gayan@g.ucla.com",
    description="A single-cell multi-omics read simulator.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Dominic7227/scReadSim",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    include_package_data=True,
    packages=['scReadSim'],
    package_dir={'scReadSim': 'scReadSim'},
    package_data={'scReadSim': ['Rscript/*.R'],
    'scReadSim': ['data/*']}
)