# scReadSim for 10x scRNA-seq with condition effect


This tutorial's main steps and corresponding estimated time usage are as follows (tested on a server with the 256x Intel Xeon Phi CPU 7210 at 1.30 GHz):

<!-- - [Step 1: Import packages and data files](#step-1-import-packages-and-data-files): < 1 min
- [Step 2: Generate features](#step-2-generate-features): < 1 min
- [Step 3: Generate real count matrices](#step-3-generate-real-count-matrices): ~ 3 mins
- [Step 4: Simulate synthetic count matrix](#step-4-simulate-synthetic-count-matrix): ~ 8 mins
- [Step 5: Output synthetic read](#step-5-output-synthetic-read): ~ 3 mins -->
- **Step 1: Import packages and data files**: < 1 min
- **Step 2: Generate features**: < 1 min
- **Step 3: Generate real count matrices**: ~ 3 mins
- **Step 4: Simulate synthetic count matrix**: ~ 10 mins
- **Step 5: Output synthetic read**: ~ 3 mins

By default, this tutorial uses Python (Python >= 3.8). However, we also include code chunks using bash commands to preprocess necessary files. To avoid users' confusion, bash commands start with a symbol **$**. We also indicate when a following code chunk is using bash commands. 


## Required softwares for scReadSim
scReadSim requires users to pre-install the following softwares:
- [samtools >= 1.12](http://www.htslib.org/)
- [bedtools >= 2.29.1](https://bedtools.readthedocs.io/en/latest/)
- [seqtk >= 1.3](https://github.com/lh3/seqtk)
- [fgbio >= 2.0.1](https://github.com/fulcrumgenomics/fgbio)

Depending on users' choices, the following softwares are optional:
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)


## Pre-process BAM file before scReadSim
**Note**: This tutorial does not need this pre-process step since the processed BAM file is provided by the scReadSim package (see below **Step 1: Import packages and data files**).

Input BAM file for scReadSim needs pre-processing to add the cell barcode in front of the read name. For example, in 10x sequencing data, cell barcode `TGGACCGGTTCACCCA-1` is stored in the field `CB:Z:TGGACCGGTTCACCCA-1`. 

The following code chunk (**bash commands**) outputs a read record from the original BAM file.

```{code-block} console
$ samtools view unprocess.bam | head -n 1
A00836:472:HTNW5DMXX:1:1372:16260:18129      83      chr1    4194410 60      50M     =       4193976 -484    TGCCTTGCTACAGCAGCTCAGGAAATGTCTTTGTGCCCACAGTCTGTGGT   :FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF      NM:i:0  MD:Z:50 AS:i:50 XS:i:0  CR:Z:TCCGGGACAGCTAACA   CY:Z:FFFFFFFFFFFFFFF:   CB:Z:TGGACCGGTTCACCCA-1 BC:Z:AAACTCAT        QT:Z::FFFFFFF   RG:Z:e18_mouse_brain_fresh_5k:MissingLibrary:1:HTNW5DMXX:1
```

The following code chunk (**bash commands**) adds the cell barcodes in front of the read names.

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
The example deploys scReadSim on the [10x single cell RNA-seq](https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0) dataset. For user convienience, we prepared the indexed reference genome files (by bowtie2), which can be downloaded using the following bash commands:
- GENCODE reference genome FASTA file and index file(indexed by bowtie2): *reference.genome.chr1.tar.gz*
- GENCODE genome annotation gtf file: *gencode.vM10.annotation.gtf*

**Note**: users may need to edit the code by using their own path. The following code chunk is using **bash commands**.


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
import scReadSim.GenerateSyntheticCount_ConditionEffect as GenerateSyntheticCount_ConditionEffect
import scReadSim.scRNA_GenerateBAM as scRNA_GenerateBAM
import pkg_resources
```

The real BAM file and other input files are listed and can be accessed by simply loading the code chunk below:
-  BAM file: *10X_RNA_chr1_3073253_4526737.bam*
-  cell barcode file: *barcodes.tsv*
-  cell condition file: *10X_RNA_conditionlabel.txt*
-  chromosome size file: *mm10.chrom.sizes*

**Note**:  To generate synthetic cells with condition effects, scReadSim requires users to input a real-cell condition label file. Although the 10x scRNA-seq dataset used in this tutorial does not contain multiple cell conditions, we generated randomly a real-cell condition label file for illustration purpose. The condition label file is already embedded in scReadSim package data and can be access through the following code chunk.

```{code-block} python3
INPUT_cells_barcode_file = pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv') 
INPUT_cells_condition_file = pkg_resources.resource_filename("scReadSim", 'data/10X_RNA_conditionlabel.txt')

filename = "10X_RNA_chr1_3073253_4526737"
INPUT_bamfile = pkg_resources.resource_filename("scReadSim", 'data/%s.bam' % filename)
INPUT_genome_size_file = pkg_resources.resource_filename("scReadSim", 'data/mm10.chrom.sizes')
```


## Step 2: Generate features
To pre-process real scRNA-seq data for training, scReadSim requires a BAM file (containing scRNA-seq reads in cells) and a gene annotation file (in GTF format). Based on the gene coordinates in the annotation file, scReadSim segregates the reference genome into two sets of features: genes and inter-genes.


### Specify output directory
**Note**: users may need to edit the code by using their own path.

```{code-block} python3
outdirectory = "/home/users/example/outputs" # may use user's own path
os.mkdir(outdirectory)
```


### Specify pre-installed software paths
**Note**: users may need to edit the code by using their own path.


```{code-block} python3
samtools_directory="/home/users/Tools/samtools" 
macs3_directory="/home/users/Tools/MACS3/bin"
bedtools_directory="/home/users/Tools/bedtools/bedtools2/bin"
seqtk_directory="/home/users/Tools/seqtk"
fgbio_jarfile="/home/users/Tools/fgbio/target/scala-2.13/fgbio-2.0.1-e884860-SNAPSHOT.jar"
```

### Prepare features
Given the input BAM file and gene annotation file, scReadSim prepares the bed files for features using function `scRNA_CreateFeatureSets`. This function needs user to specify

- `INPUT_bamfile`: Input BAM file for anlaysis.
- `samtools_directory`: Path to software *samtools*.
- `bedtools_directory`: Path to software *bedtools*.
- `outdirectory`: Specify the output directory of the features files.
- `genome_annotation`: Genome annotation file for the reference genome that the input BAM aligned on or the synthetic BAM should align on.
- `genome_size_file`: Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicating the size.

This function will generate the following two bed files into directory `outdirectory` for following analysis:

- gene bed file: *scReadSim.Gene.bed*
- inter-gene bed file: *scReadSim.InterGene.bed*

**Note**: users may need to edit the code by using their own path.

```{code-block} python3
# Specify the absolute path to gene annotation file
INPUT_genome_annotation = "/home/users/example/refgenome_dir/gencode.vM10.annotation.gtf" # may use user's own path

# Generate features
Utility.scRNA_CreateFeatureSets(INPUT_bamfile=INPUT_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_annotation=INPUT_genome_annotation, genome_size_file=INPUT_genome_size_file)
```


## Step 3: Generate real count matrices

Based on the feature sets output in **Step 2**, scReasSim constructs the UMI count matrices for genes and inter-genes through function `Utility.scRNA_bam2countmat_paral`. This function needs user to specify

- `cells_barcode_file`: Cell barcode file corresponding to the input BAM file.
- `bed_file`: Features bed file to generate the count matrix.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Specify the output directory of the count matrix file.
- `count_mat_filename`: Specify the base name of output read (or UMI) count matrix.
- `UMI_modeling`: (Optional, default: True) Specify whether scReadSim should model UMI count of the input BAM file.
- `UMI_tag`: (Optional, default: 'UB:Z') If UMI_modeling is set to True, specify the UMI tag of input BAM file, default value 'UB:Z' is the UMI tag for 10x scRNA-seq.
- `n_cores`: (Optional, default: '1') Specify the number of cores for parallel computing when generating count matrix.

For the user specified `count_mat_filename`, scReadSim will generate a count matrix named *`count_mat_filename`.txt* to directory `outdirectory`.

```{code-block} python3
# Specify the path to bed files generated by Utility.scRNA_CreateFeatureSets
gene_bedfile = outdirectory + "/" + "scReadSim.Gene.bed"
intergene_bedfile = outdirectory + "/" + "scReadSim.InterGene.bed"
# Specify the output count matrices' prenames
UMI_gene_count_mat_filename = "%s.gene.countmatrix" % filename
UMI_intergene_count_mat_filename = "%s.intergene.countmatrix" % filename

# Construct count matrix for genes
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=gene_bedfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=UMI_gene_count_mat_filename, UMI_modeling=True, UMI_tag = "UB:Z", n_cores=8)
# Construct count matrix for inter-genes
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=intergene_bedfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=UMI_intergene_count_mat_filename, UMI_modeling=True, UMI_tag = "UB:Z", n_cores=8)
```


## Step 4: Synthetic count matrix simulation
In this tutorial, scReadSim implements [scDesign3](https://github.com/SONGDONGYUAN1994/scDesign3) to generate synthetic count matrix with condition effect based on the constructed count matrix from the input BAM file. Use function `GenerateSyntheticCount_ConditionEffect.ConditionEffect_GenerateSyntheticCount` to generate synthetic count matrix with following paramters

- `count_mat_filename`: Base name of the count matrix output by function `Utility.scRNA_bam2countmat_paral`.
- `cellcondition_file`: Specify the one-column text file containing the predefined cell conditions. Make sure that the order of cell conditions correspond to the cell barcode file (and the columns of real count matrix).
- `directory`: Path to the count matrix.
- `outdirectory`: Output directory of coordinate files.
- `n_cell_new`: (Optional, default: None) Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
- `total_count_new`: (Optional, default: None) Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
- `celllabel_file`: (Optional, default: None) Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.
- `n_cores`: (Optional, default: '1') Specify the number of cores for parallel computing when generating count matrix.

Given the input count matrix *`count_mat_filename`.txt*, scReadSim generates the syntheitic count matrix file to `outdirectory` for following analysis:

- Synthetic count matrix: *`count_mat_filename`.ConditionEffect.scDesign3Simulated.txt*
- Synthetic cell cluster/type labels: *ConditionEffect.scDesign3Simulated.CellTypeLabel.txt*
- Synthetic cell condition labels: *ConditionEffect.scDesign3Simulated.ConditionLabel.txt*

Additionaly, if no `celllabel_file` is specified, scReadSim automatically performs Louvain clustering from Seurat and outputs clustering labels to `outdirectory`:
- Real cells' Louvain clustering labels: *`count_mat_filename`.LouvainClusterResults.txt*


```{code-block} python3
# Generate synthetic count matrix for gene-by-cell count matrix
GenerateSyntheticCount_ConditionEffect.ConditionEffect_GenerateSyntheticCount(count_mat_filename=UMI_gene_count_mat_filename, cellcondition_file=INPUT_cells_condition_file, directory=outdirectory, outdirectory=outdirectory)

# Specify cluster labels obtained from peak-by-cell matrix
celllabel_file = outdirectory + "/" + "10X_RNA_chr1_3073253_4526737.gene.countmatrix.LouvainClusterResults.txt"
# Generate synthetic count matrix for non-gene-by-cell count matrix
GenerateSyntheticCount_ConditionEffect.ConditionEffect_GenerateSyntheticCount(count_mat_filename=UMI_intergene_count_mat_filename, cellcondition_file=INPUT_cells_condition_file, directory=outdirectory, outdirectory=outdirectory, celllabel_file=celllabel_file)
```

## Step 5: Output synthetic read

### Generate synthetic reads in BED format
Based on the synthetic count matrix, scReadSim generates synthetic reads by randomly sampling from the real BAM file input by users. First use function `scRNA_GenerateBAMCoord` to create the synthetic reads and output in BED file storing the coordinates information. Function `scRNA_GenerateBAMCoord` takes following input arguments:
- `bed_file`: Features' bed file to generate the synthetic reads (Generated by function `Utility.scRNA_CreateFeatureSets`).
- `UMI_count_mat_file`: The path to the **synthetic UMI count matrix** generated by `GenerateSyntheticCount.scRNA_GenerateSyntheticCount`.
- `synthetic_cell_label_file`: Synthetic cell label file generated by `scRNA_GenerateSyntheticCount`.
- `read_bedfile_prename`: Specify the base name of output bed file.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Specify the output directory for synthetic reads bed file.
- `OUTPUT_cells_barcode_file`: Specify the file name storing the synthetic cell barcodes.
- `jitter_size`: (Optional, default: '5') Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
- `read_len`: (Optional, default: '50') Specify the length of synthetic reads. Default value is 50 bp.

This function will output a bed file *`read_bedfile_prename`.read.bed* storing the coordinates information of synthetic reads and its cell barcode file `OUTPUT_cells_barcode_file` in directory `outdirectory`.

After generation of synthetic reads for genes and inter-genes, combine the two bed files using function `scRNA_GenerateBAM.scRNA_CombineBED`, which takes following input arguments:
- `outdirectory`: Directory of `gene_read_bedfile_prename`.txt and `intergene_read_bedfile_prename`.txt.
- `gene_read_bedfile_prename`: File prename of foreground synthetic reads bed file.
- `intergene_read_bedfile_prename`: File prename of background synthetic reads bed file.
- `BED_filename_combined_pre`: Specify the combined syntehtic reads bed file prename. The combined bed file will be output to `outdirectory`.

```{code-block} python3
# Specify the names of synthetic count matrices (generated by GenerateSyntheticCount.scRNA_GenerateSyntheticCount)
synthetic_countmat_gene_file = UMI_gene_count_mat_filename + ".scDesign2Simulated.txt"
synthetic_countmat_intergene_file = UMI_intergene_count_mat_filename + ".scDesign2Simulated.txt"
# Specify the base name of bed files containing synthetic reads
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
gene_read_bedfile_prename = "%s.syntheticBAM.gene" % filename
intergene_read_bedfile_prename = "%s.syntheticBAM.intergene" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined" % filename
synthetic_cell_label_file = UMI_gene_count_mat_filename + ".scDesign2Simulated.CellTypeLabel.txt"

# Create synthetic read coordinates for genes
scRNA_GenerateBAM.scRNA_GenerateBAMCoord(
        bed_file=gene_bedfile, UMI_count_mat_file=outdirectory + "/" + synthetic_countmat_gene_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=gene_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=90)
# Create synthetic read coordinates for intergenes
scRNA_GenerateBAM.scRNA_GenerateBAMCoord(
        bed_file=intergene_bedfile, UMI_count_mat_file=outdirectory + "/" + synthetic_countmat_intergene_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=intergene_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=90)

# Combine bed files
scRNA_GenerateBAM.scRNA_CombineBED(outdirectory=outdirectory, gene_read_bedfile_prename=gene_read_bedfile_prename, intergene_read_bedfile_prename=intergene_read_bedfile_prename, BED_filename_combined_pre=BED_filename_combined_pre)
```

### Convert BED files to FASTQ files
Use function `scRNA_GenerateBAM.scRNA_BED2FASTQ` to convert BED file to FASTQ file. This function takes the following arguments:
- `bedtools_directory`: Path to software bedtools.
- `seqtk_directory`: Path to software seqtk.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Output directory of the synthteic bed file and its corresponding cell barcodes file.
- `BED_filename_combined`: Base name of the combined bed file output by function `scRNA_CombineBED`.
- `synthetic_fastq_prename`: Specify the base name of the output FASTQ files.

This function will output paired-end reads in FASTQ files named as *`BED_filename_combined`.read1.bed2fa.sorted.fq*, *`BED_filename_combined`.read2.bed2fa.sorted,fq* to directory `outdirectory`.

**Note**: users may need to edit the code by using their own path.

```{code-block} python3
referenceGenome_name = "chr1"
referenceGenome_dir = "/home/users/example/refgenome_dir"  # may use users' own path
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)
synthetic_fastq_prename = BED_filename_combined_pre

# Convert combined bed file into FASTQ files
scRNA_GenerateBAM.scRNA_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=BED_filename_combined_pre, synthetic_fastq_prename=synthetic_fastq_prename)
```


### Introduce Error to synthetic data
Use function `scRNA_ErrorBase` to introduce random error to synthetic reads. 

**Build reference genome dictionary (optional)**

Note that before using function `scRNA_ErrorBase`, please create the reference dictionary for the reference genome with function `CreateSequenceDictionary` using software Picard and make sure that the output *.dict* files are within the same directory to *`referenceGenome_name`.fa*. 

**Note**: For this tutorial, no dictionary building is needed since we have built for *chr1.fa* in *reference.genome.chr1.tar.gz*. The following code chunk is using **bash commands**.

```{code-block} console
$ cd /home/users/example/refgenome_dir # may use users' own path
$ java -jar /home/users/picard/build/libs/picard.jar CreateSequenceDictionary \
$       -R chr1.fa \
$       -O chr1.fa.dict
```

**Introduce errors to synthetic reads**

Function `scRNA_ErrorBase` takes the following arguments:
- `fgbio_jarfile`: Path to software fgbio jar script.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Specify the output directory of the synthteic FASTQ file with random errors.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scRNA_BED2FASTQ`.

This function will output synthetic reads with random errors in FASTQ files named as *`synthetic_fastq_prename`.ErrorIncluded.read1.bed2fa.fq*, *`synthetic_fastq_prename`.ErrorIncluded.read2.bed2fa.fq* to directory `outdirectory`.


```{code-block} python3
# Generate reads with errors in FASTQs
scRNA_GenerateBAM.scRNA_ErrorBase(fgbio_jarfile=fgbio_jarfile, INPUT_bamfile=INPUT_bamfile, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, synthetic_fastq_prename=synthetic_fastq_prename)
```

### Convert FASTQ files to BAM file (optional)
The current version of scReadSim implicitly uses [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align the synthetic reads onto the reference genome. Use function `scRNA_GenerateBAM.AlignSyntheticBam_Single` to align FASTQ files onto reference genome. It takes the following arguments:
- `bowtie2_directory`: Path to software bowtie2.
- `samtools_directory`: Path to software samtools.
- `outdirectory`: Specify the output directory of the synthteic BAM file.
- `referenceGenome_name`: Base name of the reference genome FASTA file. For example, you should input "chr1" for file "chr1.fa".
- `referenceGenome_dir`: Path to the reference genome FASTA file.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scRNA_BED2FASTQ`.
- `output_BAM_pre`: Specify the base name of the output BAM file.

**Index reference genome (optional)** 

Before using function `AlignSyntheticBam_Pair`, the reference gemome FASTA file should be indexed by bowtie2 through following chunk and make sure the output index files are within the same directory to *`referenceGenome_name`.fa*. 

**Note**: For this tutorial, no indexing is needed since we have indexed *chr1.fa* in *reference.genome.chr1.tar.gz*. The following code chunk is using **bash commands**. 



```{code-block} console
$ cd /home/users/example/refgenome_dir # change to directory where your reference genome file is
$ bowtie2-build chr1.fa chr1
```

**Align synthetic reads** 

Now align the synthetic reads on to the reference genome with bowtie2.


```{code-block} python3
# Specify bowtie2 path
bowtie2_directory="/home/users/Tools/bowtie2/bin"
# Specify output BAM name
output_BAM_pre = "%s.syntheticBAM.CBincluded" % filename

# Synthetic reads alignment
scRNA_GenerateBAM.AlignSyntheticBam_Single(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename, output_BAM_pre=output_BAM_pre)

# Synthetic reads (with sequencing errors) alignment
scRNA_GenerateBAM.AlignSyntheticBam_Single(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename + ".ErrorIncluded" , output_BAM_pre=output_BAM_pre+ ".ErrorIncluded")
```