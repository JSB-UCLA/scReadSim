import setuptools
from pathlib import Path
import os
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

install_require_list = [
    "numpy",
    "pandas",
    "pysam",
    "tqdm",
    "joblib", 
    "pathlib",
    "Bio",
    "gffpandas"
    ]
# we need to exclude rpy2 when building the docs, and mock it for import in docs/source/conf.py
# using the autodoc_mock_imports parameter:
if not os.getenv('READTHEDOCS'):
    install_require_list.append('rpy2')

setuptools.setup(
    name="scReadSim",
    version="1.4.1",
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
    install_requires=install_require_list,
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