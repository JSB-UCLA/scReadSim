# scReadSim on 10x scRNA-seq from multiple samples

This tutorial demonstrates the application of scReadSim generating synthetic reads for **scRNA-seq from multiple replicates/samples**. For simplicity, we chose two of 10x datasets as the input multi-samples for this tutorial illustration. The datasets we chose include
- [10x Single Cell Multiome dataset (RNA modality) for Fresh Embryonic E18 Mouse Brain](https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0)
- [10x E18 Mouse Brain scRNA-seq](https://www.10xgenomics.com/resources/datasets/5-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-cells-3-1-standard-6-0-0)


This tutorial's main steps and corresponding estimated time usage are as follows (tested on a server with the 256x Intel Xeon Phi CPU 7210 at 1.30 GHz):

<!-- - [Step 1: Import packages and data files](#step-1-import-packages-and-data-files): < 1 min
- [Step 2: Generate features](#step-2-generate-features): < 1 min
- [Step 3: Generate real count matrices](#step-3-generate-real-count-matrices): ~ 5 mins
- [Step 4: Simulate synthetic count matrix](#step-4-simulate-synthetic-count-matrix): ~ 15 mins
- [Step 5: Output synthetic read](#step-5-output-synthetic-read): ~ 11 mins -->
- **Step 1: Import packages and data files**: < 1 min
- **Step 2: Generate features**: < 1 min
- **Step 3: Generate real count matrices**: ~ 5 mins
- **Step 4: Simulate synthetic count matrix**: ~ 15 mins
- **Step 5: Output synthetic read**: ~ 11 mins




## Required softwares for scReadSim
scReadSim requires users to pre-install the following softwares:
- [samtools](http://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [seqtk](https://github.com/lh3/seqtk)
- [fgbio](http://fulcrumgenomics.github.io/fgbio/)

Depending on users' choices, the following softwares are optional:
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)



## Pre-process input BAM file
**Note: This tutorial does not need this pre-process step since the processed BAM file is provided by the scReadSim package (see [Step 1: Import packages and data files](#step-1-import-packages-and-data-files)).**

Input BAM file for scReadSim needs pre-processing to add the cell barcode in front of the read name. For example, in 10x sequencing data, cell barcode `TGGACCGGTTCACCCA-1` is stored in the field `CB:Z:TGGACCGGTTCACCCA-1`. 

```{code-block} console
$ samtools view unprocess.bam | head -n 1
A00836:472:HTNW5DMXX:1:1372:16260:18129      83      chr1    4194410 60      50M     =       4193976 -484    TGCCTTGCTACAGCAGCTCAGGAAATGTCTTTGTGCCCACAGTCTGTGGT   :FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF      NM:i:0  MD:Z:50 AS:i:50 XS:i:0  CR:Z:TCCGGGACAGCTAACA   CY:Z:FFFFFFFFFFFFFFF:   CB:Z:TGGACCGGTTCACCCA-1 BC:Z:AAACTCAT        QT:Z::FFFFFFF   RG:Z:e18_mouse_brain_fresh_5k:MissingLibrary:1:HTNW5DMXX:1
```

The following code chunk adds the cell barcodes in front of the read names.

```{code-block} console
$ # extract the header file
$ mkdir tmp
$ samtools view unprocess.bam -H > tmp/unprocess.header.sam

$ # create a bam file with the barcode embedded into the read name
$ time(cat <( cat tmp/unprocess.header.sam ) \
 <( samtools view unprocess.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
 | samtools view -bS - > processed.bam) 
$ rm -dr tmp

$ samtools view processed.bam | head -n 1
TGGACCGGTTCACCCA-1:A00836:472:HTNW5DMXX:1:1372:16260:18129      83      chr1    4194410 60      50M     =       4193976 -484    TGCCTTGCTACAGCAGCTCAGGAAATGTCTTTGTGCCCACAGTCTGTGGT   :FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF      NM:i:0  MD:Z:50 AS:i:50 XS:i:0  CR:Z:TCCGGGACAGCTAACA   CY:Z:FFFFFFFFFFFFFFF:   CB:Z:TGGACCGGTTCACCCA-1 BC:Z:AAACTCAT        QT:Z::FFFFFFF   RG:Z:e18_mouse_brain_fresh_5k:MissingLibrary:1:HTNW5DMXX:1
```



## Download reference genome for test example
For user convienience, we prepared the indexed reference genome files (by bowtie2), which can be downloaded using the following bash commands:
- GENCODE reference genome FASTA file and index file(indexed by bowtie2): reference.genome.chr1.tar.gz
- GENCODE genome annotation gtf file: gencode.vM10.annotation.gtf

**Note**: users may need to edit the code by using their own path.


```{code-block} console
$ mkdir /home/users/example/refgenome_dir # may use users' own path
$ cd /home/users/example/refgenome_dir
$ wget http://compbio10data.stat.ucla.edu/repository/gayan/Projects/scReadSim/reference.genome.chr1.tar.gz # 292 MB
$ wget http://compbio10data.stat.ucla.edu/repository/gayan/Projects/scReadSim/gencode.vM10.annotation.gtf # 765 MB
$ tar -xf reference.genome.chr1.tar.gz
```



## Step 1: Import packages and data files

Import modules.

```{code-block} python3
import sys, os
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scATAC_GenerateBAM as scATAC_GenerateBAM
import pkg_resources
```

The real BAM files and other input files are listed and can be accessed by simply loading the code chunk below:
-  Sample 1 BAM file: 10X_ATAC_chr1_4194444_4599104_rep1.bam
-  Sample 2 BAM file: 10X_ATAC_chr1_4194444_4599104_rep2.bam
-  Sample 1 cell barcode file: barcodes_top1k_ATACrep1.tsv
-  Sample 2 cell barcode file: barcodes_top1k_ATACrep2.tsv
-  chromosome size file: mm10.chrom.sizes

```{code-block} python3
INPUT_genome_size_file = pkg_resources.resource_filename("scReadSim", 'data/mm10.chrom.sizes')

# List of input cell barcode files
INPUT_cells_barcode_file = [pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv'), 
                            pkg_resources.resource_filename("scReadSim", 'data/barcodes_RNArep2.tsv')]

# List of input BAM files
INPUT_bamfile = [pkg_resources.resource_filename("scReadSim", 'data/10X_RNA_chr1_3073253_4526737.bam'),
                 pkg_resources.resource_filename("scReadSim", 'data/10X_RNA_chr1_3073253_4526737_rep2.bam')]
```



## Step 2: Generate features
To pre-process real scRNA-seq data for training, scReadSim requires a BAM file (containing scRNA-seq reads in cells) and a gene annotation file (in GTF format). Based on the gene coordinates in the annotation file, scReadSim segregates the reference genome into two sets of features: genes and inter-genes.


### Specify scReadSim working directory for outputs
**Note**: users may need to edit the code by using their own path.


```{code-block} python3
outdirectory = "/home/users/example/outputs" # may use user's own path
os.mkdir(outdirectory)
```

### Specify pre-installed software paths
**Note**: users may need to edit the code by using their own path.

```{code-block} python3
samtools_directory="/home/users/Tools/samtools/bin" 
bedtools_directory="/home/users/Tools/bedtools/bedtools2/bin"
seqtk_directory="/home/users/Tools/seqtk/bin"
fgbio_jarfile="/home/users/Tools/fgbio/target/scala-2.13/fgbio-2.0.1-e884860-SNAPSHOT.jar"
```

### Prepare Features

To prepare features for the following analysis, scReadSim utilizes function `Utility.scRNA_CreateFeatureSets_MultiSample` with following arguments

- `INPUT_bamfile`: List of input BAM files (use absolute paths to the BAM files).
- `samtools_directory`: Directory of software samtools.
- `bedtools_directory`: Directory of software bedtools.
- `outdirectory`: Specify the working directory of scReadSim for generating intermediate and final output files.
- `genome_size_file`: Directory of Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicates the size. 
- `genome_annotation`: Genome annotation file for the reference genome that the input BAM aligned on or the synthetic BAM should align on.


Function `Utility.scRNA_CreateFeatureSets_MultiSample` creates multiple sub-directories under `outdirectory` for each sample/replicate. Within each sub-directory, the following two bed files will be generated for following analysis:

- gene bed file: *scReadSim.Gene.bed*
- inter-gene bed file: *scReadSim.InterGene.bed*

```{code-block} python3
Utility.scRNA_CreateFeatureSets_MultiSample(INPUT_bamfile=INPUT_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_size_file=INPUT_genome_size_file, genome_annotation=INPUT_genome_annotation) # Less than one minute 
```



## Step 3: Generate real count matrices
Based on the feature sets output in **Step 2**, scReasSim constructs the count matrices for peaks and non-peaks through function `Utility.scATAC_bam2countmat_paral_MultiSample`. This function needs user to specify

- `cells_barcode_file`: List of cell barcode files corresponding to the input BAM files.
- `INPUT_bamfile`: List of input BAM files (use absolute paths to the BAM files).
- `outdirectory`: Specify the working directory of scReadSim for generating intermediate and final output files.
- `n_cores`: (Optional, default: '1') Specify the number of cores for parallel computing when generating count matrix.
- `UMI_modeling`: (Optional, default: True) Specify whether scReadSim should model UMI count of the input BAM file.
- `UMI_tag`: (Optional, default: 'UB:Z') If UMI_modeling is set to True, specify the UMI tag of input BAM file, default value 'UB:Z' is the UMI tag for 10x scRNA-seq.

For eac sample/replicate in its corresponding sub-directory, the following files will be generated for following analysis:

- Gene-by-cell count matrix: *RepX.gene.countmatrix.txt*
- Inter-gene-by-cell count matrix: *RepX.intergene.countmatrix.txt*

```{code-block} python3
Utility.scRNA_bam2countmat_paral_MultiSample(cells_barcode_file=INPUT_cells_barcode_file, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, n_cores=8)
```



## Step 4: Simulate synthetic count matrix
In this tutorial, scReadSim implements [scDesign2](https://github.com/JSB-UCLA/scDesign2) to generate synthetic count matrix based on the constructed count matrix from the input BAM files. Use function `GenerateSyntheticCount.scRNA_GenerateSyntheticCount_MultiSample` to generate synthetic count matrix with following paramters

- `INPUT_bamfile`: List of input BAM files (use absolute paths to the BAM files).
- `outdirectory`: Specify the working directory of scReadSim for generating intermediate and final output files.
- `n_cell_new`: (Optional, default: 'None') Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
- `total_count_new`: (Optional, default: 'None') Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
- `n_cores`: (Optional, default: '1') Specify the number of cores for parallel computing when generating count matrix.

For eac sample/replicate in its corresponding sub-directory, the following files will be generated for following analysis:

- Synthetic peak-by-cell count matrix: *RepX.gene.countmatrix.scDesign2Simulated.txt*
- Synthetic peak-by-cell count matrix: *RepX.intergene.countmatrix.scDesign2Simulated.txt*
- Synthetic cell cluster/type labels: *RepX.gene.countmatrix.scDesign2Simulated.CellTypeLabel.txt*
- Real cells' Louvain clustering labels: *RepX.gene.countmatrix.LouvainClusterResults.txt*



```{code-block} python3
GenerateSyntheticCount.scRNA_GenerateSyntheticCount_MultiSample(INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, n_cores=8)
```



## Step 5: Output synthetic read

Use function `scRNA_GenerateBAM.scATAC_GenerateSyntheticRead_MultiSample` to generate synthetic reads based on the synthetic count matrices generated from last step. The function takes following input arguments:

- `INPUT_bamfile`: List of input BAM files (use absolute paths to the BAM files).
- `outdirectory`: Specify the working directory of scReadSim for generating intermediate and final output files.
- `bedtools_directory`: Directory of software bedtools.
- `seqtk_directory`: Directory of software seqtk.
- `fgbio_jarfile`: Path to software fgbio jar script.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `read_len`: (Optional, default: '90') Specify the length of synthetic reads. Default value is 90 bp.


### Build reference genome dictionary (optional)

Note that before using function `scRNA_GenerateBAM.scRNA_GenerateSyntheticRead_MultiSample`, please create the reference dictionary for the reference genome with function `CreateSequenceDictionary` using software Picard and make sure that the output *.dict* files are within the same directory to *`referenceGenome_name`.fa*. **For this tutorial, no dictionary building is needed since we have built for chr1.fa in reference.genome.chr1.tar.gz**. 

```{code-block} console
$ cd /home/users/example/refgenome_dir # may use users' own path
$ java -jar /home/users/picard/build/libs/picard.jar CreateSequenceDictionary \
$       -R chr1.fa \
$       -O chr1.fa.dict
```

### Generate synthetic read
For eac sample/replicate in its corresponding sub-directory, the following files will be generated by function `scRNA_GenerateBAM.scRNA_GenerateSyntheticRead_MultiSample`:

- Synthetic read 1 FASTQ file: *RepX.syntheticBAM.combined.read1.bed2fa.sorted.fq*
- Synthetic read 2 FASTQ file: *RepX.syntheticBAM.combined.read2.bed2fa.sorted.fq*
- Synthetic read 1 FASTQ file with substitutional error: *RepX.syntheticBAM.combined.ErrorIncluded.read1.bed2fa.sorted.fq*
- Synthetic read 2 FASTQ file with substitutional error: *RepX.syntheticBAM.combined.ErrorIncluded.read2.bed2fa.sorted.fq*
- Synthetic cell barcode file: *synthetic_cell_barcode_repX.txt*


**Note**: users may need to edit the code by using their own path.

```{code-block} python3
# Specify reference genome file path, may use users' own path
referenceGenome_name = "chr1"
referenceGenome_dir = "/home/users/example/refgenome_dir" 
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)

scRNA_GenerateBAM.scRNA_GenerateSyntheticRead_MultiSample(INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, fgbio_jarfile=fgbio_jarfile, referenceGenome_file=referenceGenome_file)
```


