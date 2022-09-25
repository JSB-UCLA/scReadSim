import setuptools
from pathlib import Path

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scReadSim",
    version="1.2.0",
    author="Guanao Yan",
    author_email="gayan@g.ucla.com",
    description="A single-cell RNA-seq and ATAC-seq read simulator.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JSB-UCLA/scReadSim",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
    "numpy",
    "pandas",
    "pysam",
    "tqdm",
    "rp2",
    "joblib", 
    "pathlib"
    ],
    # install_requires=[
    #     x.strip() for x in
    #     Path('requirements_test.txt').read_text('utf-8').splitlines()
    # ],
    # include_package_data=True,
    # # packages=['scReadSim'],
    packages=setuptools.find_packages(),
    # package_dir={'': 'scReadSim'},
    package_data={
    'scReadSim': ['data/*', 'Rscript/*']}
)