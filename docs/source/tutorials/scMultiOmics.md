# scReadSim on 10x MultiOme ATAC + Gene Expression Dataset 

This tutorial demonstrates the application of scReadSim generating synthetic reads for **single-cell multiomics** with [10x Single Cell Multiome ATAC + Gene Expression Dataset](https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0). The current version of scReadSim multiomic module requires the input of two modalities, scATAC-seq and scRNA-seq, to have matchinging single cells.

This tutorial's main steps and corresponding estimated time usage are as follows (tested on a server with the 256x Intel Xeon Phi CPU 7210 at 1.30 GHz):

<!-- - [Step 1: Import packages and data files](#step-1-import-packages-and-data-files): < 1 min
- [Step 2: Generate features](#step-2-generate-features): < 1 min
- [Step 3: Generate real count matrices](#step-3-generate-real-count-matrices): ~ 3 mins
- [Step 4: Simulate synthetic count matrix](#step-4-simulate-synthetic-count-matrix): ~ 6 mins
- [Step 5: Output synthetic read](#step-5-output-synthetic-read): ~ 8 mins -->
- **Step 1: Import packages and data files**: < 1 min
- **Step 2: Generate features**: < 1 min
- **Step 3: Generate real count matrices**: ~ 3 mins
- **Step 4: Simulate synthetic count matrix**: ~ 6 mins
- **Step 5: Output synthetic read**: ~ 8 mins

By default, this tutorial uses Python (Python >= 3.8). However, we also include code chunks using bash commands to preprocess necessary files. To avoid users' confusion, bash commands start with a symbol **$**. We also indicate when a following code chunk is using bash commands. 


## Required softwares for scReadSim
scReadSim requires users to pre-install the following softwares:
- [MACS3](https://github.com/macs3-project/MACS)
- [samtools](http://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [seqtk](https://github.com/lh3/seqtk)
- [fgbio](http://fulcrumgenomics.github.io/fgbio/)

Depending on users' choices, the following softwares are optional:
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)


## Pre-process input BAM file
**Note**: This tutorial does not need this pre-process step since the processed BAM file is provided by the scReadSim package (see **Step 1: Import packages and data files**).

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
The example deploys scReadSim on the [10x Single Cell Multiome ATAC + Gene Expression Dataset](https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0). For user convienience, we prepared the indexed reference genome files (by bowtie2), which can be downloaded using the following bash commands:
- GENCODE reference genome FASTA file and index file(indexed by bowtie2): reference.genome.chr1.tar.gz
- GENCODE genome annotation gtf file: gencode.vM10.annotation.gtf

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
import scReadSim.scRNA_GenerateBAM as scRNA_GenerateBAM
import scReadSim.scATAC_GenerateBAM as scATAC_GenerateBAM
import scReadSim.GenerateSyntheticCount_MultiOmics as GenerateSyntheticCount_MultiOmics
import pkg_resources
```

The real BAM file and other input files are listed and can be accessed by simply loading the code chunk below:
-  BAM file (RNA modality): 10X_RNA_chr1_3073253_4526737.bam
-  BAM file (ATAC modality): 10X_ATAC_chr1_4194444_4399104.bam
-  cell barcode file: barcodes.tsv
-  chromosome size file: mm10.chrom.sizes

```{code-block} python3
# Load cell barcode and chrom size files
INPUT_cells_barcode_file = pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv') 
INPUT_genome_size_file = pkg_resources.resource_filename("scReadSim", 'data/mm10.chrom.sizes')

# Load RNA bam file
filename_RNA = "10X_RNA_chr1_3073253_4526737"
INPUT_RNA_bamfile = pkg_resources.resource_filename("scReadSim", 'data/%s.bam' % filename_RNA)

# Load ATAC bam file
filename_ATAC = "10X_ATAC_chr1_4194444_4399104"
INPUT_ATAC_bamfile = pkg_resources.resource_filename("scReadSim", 'data/%s.bam' % filename_ATAC)
```


## Step 2: Generate features

### Specify output directory
**Note**: users may need to edit the code by using their own path.


```{code-block} python3
outdirectory = "/home/users/example/outputs" # may use user's own path
os.mkdir(outdirectory)
```

### Specify pre-installed software paths
**Note**: users may need to edit the code by using their own path.

```{code-block} python3
samtools_directory="/home/users/Tools/samtools/bin" 
macs3_directory="/home/users/Tools/MACS3/bin"
bedtools_directory="/home/users/Tools/bedtools/bedtools2/bin"
seqtk_directory="/home/users/Tools/seqtk/bin"
fgbio_jarfile="/home/users/Tools/fgbio/target/scala-2.13/fgbio-2.0.1-e884860-SNAPSHOT.jar"
```

### Prepare Features for RNA modality

Given the input BAM file (RNA modality) and gene annotation file, scReadSim prepares the bed files for features using function `scRNA_CreateFeatureSets`. This function needs user to specify

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
Utility.scRNA_CreateFeatureSets(INPUT_bamfile=INPUT_RNA_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_annotation=INPUT_genome_annotation, genome_size_file=INPUT_genome_size_file)
```

### Prepare Features for ATAC modality
To pre-process real scATAC-seq data for training, scReadSim segments the reference genome into trustworthy peaks, trustworthy non-peaks and gray ares. First scReadSim prepares the trustworthy peaks and non-peaks for the input BAM file. Then scReadSim defines gray areas as the genomic regions complementary to the trustworthy peaks and non-peaks. Three bed files recording peaks, non-peaks and gray areas will be prepared by scReadSim for following analysis.

To prepare features for the following analysis, scReadSim utilizes function `Utility.scATAC_CreateFeatureSets` with following arguments
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `samtools_directory`: Directory of software samtools.
- `bedtools_directory`: Directory of software bedtools.
- `outdirectory`: Output directory of the prepared features.
- `genome_size_file`: Directory of Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicates the size.  
- `peak_mode`: (Optional, default: "macs3") Specify mode for trustworthy peak and non-peak generation, must be one of the following: "macs3", "user", and "superset". 
- `macs3_directory`: (Optional, default: None) Path to software MACS3. Must be specified if `INPUT_peakfile` and `INPUT_nonpeakfile` are None.
- `INPUT_peakfile`: (Optional, default: None) Directory of user-specified input peak file.
- `INPUT_nonpeakfile`: (Optional, default: None) Directory of user-specified input non-peak file.
- `superset_peakfile`: (Optional, default: None) Directory of a superset of potential chromatin open regions, including sources such as ENCODE cCRE (Candidate Cis-Regulatory Elements) collection. Must be specified under peak_mode "superset".
- `OUTPUT_peakfile`: (Optional, default: None) Directory of user-specified output peak file. Synthetic scATAC-seq reads will be generated taking `OUTPUT_peakfile` as ground truth peaks. Note that `OUTPUT_peakfile` does not name the generated feature files by function `scATAC_CreateFeatureSets`.

Three modes are supported by scReadSim to prepare features: "macs3" (default), "user" and "superset". 

**Note: This tutorial provides an example with the default `peak_mode` "macs3". Thus the following two code chunks with `peak_mode` set to "user" or "superset" do not need to be implemented.**

Under default mode "macs3" (by setting argument `peak_mode` as the default values "macs3"), scReadSim uses [MACS3](https://github.com/macs3-project/MACS) with the stringent criteria to call trustworthy peaks (q-value `0.01`) and non-peaks (q-value `0.1`) from the input BAM file. This function will generate the following three bed files into directory `outdirectory` for following analysis:

- peak bed file: *scReadSim.MACS3.peak.bed*
- non-peak bed file: *scReadSim.MACS3.nonpeak.bed*
- gray area bed file: *scReadSim.grayareas.bed*

```{code-block} python3
Utility.scATAC_CreateFeatureSets(peak_mode="macs3", INPUT_bamfile=INPUT_ATAC_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_size_file=INPUT_genome_size_file, macs3_directory=macs3_directory, INPUT_peakfile=None, INPUT_nonpeakfile=None)
```



## Step 3: Generate real count matrices

### RNA modality
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
UMI_gene_count_mat_filename = "%s.gene.countmatrix" % filename_RNA
UMI_intergene_count_mat_filename = "%s.intergene.countmatrix" % filename_RNA

## Construct count matrix for genes
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=gene_bedfile, INPUT_bamfile=INPUT_RNA_bamfile, outdirectory=outdirectory, count_mat_filename=UMI_gene_count_mat_filename, UMI_modeling=True, UMI_tag = "UB:Z", n_cores=8)
## Construct count matrix for inter-genes
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=intergene_bedfile, INPUT_bamfile=INPUT_RNA_bamfile, outdirectory=outdirectory, count_mat_filename=UMI_intergene_count_mat_filename, UMI_modeling=True, UMI_tag = "UB:Z", n_cores=8)
```


### ATAC modality

Based on the feature sets output in **Step 2**, scReasSim constructs the count matrices for both foreground feautures and background features through function `Utility.scATAC_bam2countmat_paral`. This function needs user to specify

- `cells_barcode_file`: Cell barcode file corresponding to the input BAM file.
- `bed_file`: Features' bed file to generate the count matrix (Generated by function `Utility.scATAC_CreateFeatureSets`).
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Specify the output directory of the count matrix file.
- `count_mat_filename`: Specify the base name of output count matrix.
- `n_cores`: (Optional, default: '1') Specify the number of cores for parallel computing when generating count matrix.

For the user specified `count_mat_filename`, scReadSim will generate a count matrix named  *`count_mat_filename`.txt* to directory `outdirectory`.

```{code-block} python3
# Specify the path to bed files generated by Utility.scATAC_CreateFeatureSets
peak_bedfile = outdirectory + "/" + "scReadSim.MACS3.peak.bed"
nonpeak_bedfile = outdirectory + "/" + "scReadSim.MACS3.nonpeak.bed"
# Specify the output count matrices' prenames
count_mat_peak_filename = "%s.peak.countmatrix" % filename_ATAC
count_mat_nonpeak_filename = "%s.nonpeak.countmatrix" % filename_ATAC

# Construct count matrix for peaks
Utility.scATAC_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=peak_bedfile, INPUT_bamfile=INPUT_ATAC_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_peak_filename, n_cores=1)
# Construct count matrix for non-peaks
Utility.scATAC_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=nonpeak_bedfile, INPUT_bamfile=INPUT_ATAC_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_nonpeak_filename, n_cores=1)
```


## Step 4: Simulate synthetic count matrix

In this tutorial, scReadSim implements [scDesign3](https://songdongyuan1994.github.io/scDesign3/docs/index.html) to simulate single-cell multiomics count matrix based on the real count matrices obtained from the input BAM files. Use function `GenerateSyntheticCount_MultiOmics.scMultiOmics_GenerateSyntheticCount` to generate synthetic count matrices with following paramters

- `RNA_count_mat_filename`: Base name of the count matrix output by function `scRNA_bam2countmat_paral()`.
- `ATAC_count_mat_filename`: Base name of the count matrix output by function `scATAC_bam2countmat_paral()`.
- `directory`: Path to the count matrix.
- `outdirectory`: Specify the output directory of the synthetic count matrix file.
- `doub_classification_label_file`: (Optional, default: 'None') Specify the absolute path to the doublet classification result `doublet_classification.Rdata` generated by function `DoubletDetection.detectDoublet`.
- `n_cell_new`: (Optional, default: 'None') Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
- `celllabel_file`: (Optional, default: 'None') Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign3.
- `n_cores`: (Optional, default: 1): Number of cores for parallel computing.

Given the input RNA count matrix *`RNA_count_mat_filename`.txt* and ATAC count matrix *`ATAC_count_mat_filename`.txt*, scReadSim generates the syntheitic count matrices file to `outdirectory` for following analysis:

- Synthetic count matrix (RNA modality): *`RNA_count_mat_filename`.scMultiOmics.scDesign3Simulated.ATAC.txt*
- Synthetic count matrix (ATAC modality): *`ATAC_count_mat_filename`.scMultiOmics.scDesign3Simulated.ATAC.txt*
- Synthetic cell cluster/type labels: *scMultiOmics.scDesign3Simulated.CellTypeLabel.txt*

Additionaly, if no `celllabel_file` is specified, scReadSim automatically performs Louvain clustering (from Seurat) based on the RNA modality count matrix and outputs clustering labels to `outdirectory`:
- Real cells' Louvain clustering labels: *`RNA_count_mat_filename`.LouvainClusterResults.txt*


```{code-block} python3
# Generate multiomic count matrices for genes and peaks
GenerateSyntheticCount_MultiOmics.scMultiOmics_GenerateSyntheticCount(RNA_count_mat_filename=UMI_gene_count_mat_filename, ATAC_count_mat_filename=count_mat_peak_filename, directory=outdirectory, outdirectory=outdirectory, n_cores=10)

# Specify clustering labels obtained from gene-by-cell matrix
celllabel_file = outdirectory + "/" + "10X_RNA_chr1_3073253_4526737.gene.countmatrix.LouvainClusterResults.txt"
# Generate multiomic count matrices for intergenes and non-peaks
GenerateSyntheticCount_MultiOmics.scMultiOmics_GenerateSyntheticCount(RNA_count_mat_filename=UMI_intergene_count_mat_filename, ATAC_count_mat_filename=count_mat_nonpeak_filename, directory=outdirectory, outdirectory=outdirectory, celllabel_file=celllabel_file ,n_cores=10)
```


## Step 5: Output synthetic read

### Generate synthetic reads in BED format (RNA modality)
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
synthetic_cell_label_file = "scMultiOmics.scDesign3Simulated.CellTypeLabel.txt"
synthetic_countmat_gene_file = UMI_gene_count_mat_filename + ".scMultiOmics.scDesign3Simulated.RNA.txt"
synthetic_countmat_intergene_file = UMI_intergene_count_mat_filename + ".scMultiOmics.scDesign3Simulated.RNA.txt"

# Specify the base name of bed files containing synthetic reads
gene_read_bedfile_prename = "%s.scMultiOmics.scReadSim.RNA.gene" % filename_RNA
intergene_read_bedfile_prename = "%s.scMultiOmics.scReadSim.RNA.intergene" % filename_RNA
RNA_read_bedfile_prename = "%s.scMultiOmics.scReadSim.RNA.combined" % filename_RNA

# Specify the file name of synthetic cell barcodes
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"

# Create synthetic read coordinates for genes
scRNA_GenerateBAM.scRNA_GenerateBAMCoord(
        bed_file=gene_bedfile, UMI_count_mat_file=outdirectory + "/" + synthetic_countmat_gene_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=gene_read_bedfile_prename, INPUT_bamfile=INPUT_RNA_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=90)

# Create synthetic read coordinates for intergenes
scRNA_GenerateBAM.scRNA_GenerateBAMCoord(
        bed_file=intergene_bedfile, UMI_count_mat_file=outdirectory + "/" + synthetic_countmat_intergene_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=intergene_read_bedfile_prename, INPUT_bamfile=INPUT_RNA_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=90)

# Combine synthetic read bed files
scRNA_GenerateBAM.scRNA_CombineBED(outdirectory=outdirectory, gene_read_bedfile_prename=gene_read_bedfile_prename, intergene_read_bedfile_prename=intergene_read_bedfile_prename, BED_filename_combined_pre=RNA_read_bedfile_prename)
```

### Generate synthetic reads in BED format (ATAC modality)
Based on the synthetic count matrix, scReadSim generates synthetic reads by randomly sampling from the real BAM file input by users. First use function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord` to create the synthetic reads and output in BED file storing the coordinates information. Function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord` takes following input arguments:

- `bed_file`: Features' bed file to generate the synthetic reads (Generated by function `Utility.scATAC_CreateFeatureSets`).
- `count_mat_file`: The path to the **synthetic count matrix** generated by `GenerateSyntheticCount.scATAC_GenerateSyntheticCount`.
- `synthetic_cell_label_file`: Synthetic cell label file generated by `scATAC_GenerateSyntheticCount`.
- `read_bedfile_prename`: Specify the base name of output bed file.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Specify the output directory for synthetic reads bed file.
- `OUTPUT_cells_barcode_file`: Specify the file name storing the synthetic cell barcodes.
- `jitter_size`: (Optional, default: '5') Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
- `read_len`: (Optional, default: '50') Specify the length of synthetic reads. Default value is 50 bp.
- `random_noise_mode`: (Optional, default: 'False') Specify whether to use a uniform distribution of reads.
- `GrayAreaModeling`: (Optional, default: 'False') Specify whether to generate synthetic reads for Gray Areas when generaing reads for non-peaks. Do not specify 'True' when generating reads for peaks.

This function will output two bed files *`read_bedfile_prename`.read1.bed* and *`read_bedfile_prename`.read2.bed* storing the coordinates information of synthetic reads and its cell barcode file `OUTPUT_cells_barcode_file` in directory `outdirectory`.

After generation of synthetic reads for both peaks and non-peaks, combine the their bed files using function `scATAC_GenerateBAM.scATAC_CombineBED`, which takes following input arguments:
- `outdirectory`: Directory of `peak_read_bedfile_prename`.txt and `nonpeak_read_bedfile_prename`.txt.
- `peak_read_bedfile_prename`: Base name of the bed file containig synthetic reads for peaks (generated by function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord`).
- `nonpeak_read_bedfile_prename`: Base name of the bed file containig synthetic reads for non-peaks (generated by function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord`).
- `BED_filename_combined_pre`: Specify the base name for the combined syntehtic reads bed file. The combined bed file will be output to `outdirectory`.


```{code-block} python3
# Specify the names of synthetic count matrices (generated by GenerateSyntheticCount.scATAC_GenerateSyntheticCount)
synthetic_countmat_peak_file = count_mat_peak_filename + ".scMultiOmics.scDesign3Simulated.ATAC.txt"
synthetic_countmat_nonpeak_file = count_mat_nonpeak_filename + ".scMultiOmics.scDesign3Simulated.ATAC.txt"

# Specify the base name of bed files containing synthetic reads
peak_read_bedfile_prename = "%s.scMultiOmics.scReadSim.ATAC.peak" % filename_ATAC
nonpeak_read_bedfile_prename = "%s.scMultiOmics.scReadSim.ATAC.nonpeak" % filename_ATAC
ATAC_read_bedfile_prename = "%s.scMultiOmics.scReadSim.ATAC.combined" % filename_ATAC

# Create synthetic read bed file for peaks
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(bed_file=peak_bedfile, count_mat_file=outdirectory + "/" + synthetic_countmat_peak_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=peak_read_bedfile_prename, INPUT_bamfile=INPUT_ATAC_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50)

# Ceate synthetic read bed file for non-peaks
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(bed_file=nonpeak_bedfile, count_mat_file=outdirectory + "/" + synthetic_countmat_nonpeak_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=nonpeak_read_bedfile_prename, INPUT_bamfile=INPUT_ATAC_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50,  GrayAreaModeling=True)

# Combine bed files
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory=outdirectory, peak_read_bedfile_prename=peak_read_bedfile_prename, nonpeak_read_bedfile_prename=nonpeak_read_bedfile_prename, BED_filename_combined_pre=ATAC_read_bedfile_prename)
```

### Convert BED files to FASTQ files
Use functions `scRNA_BED2FASTQ` and `scATAC_BED2FASTQ` to convert BED file to FASTQ file. These functions take the following arguments:
- `bedtools_directory`: Path to software bedtools.
- `seqtk_directory`: Path to software seqtk.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Output directory of the synthteic bed file and its corresponding cell barcodes file.
- `BED_filename_combined`: Base name of the combined bed file output by function `scRNA_CombineBED` or `scATAC_CombineBED`.
- `synthetic_fastq_prename`: Specify the base name of the output FASTQ files.

This function will output paired-end reads in FASTQ files named as *`synthetic_fastq_prename`.read1.bed2fa.sorted.fq*, *`synthetic_fastq_prename`.read2.bed2fa.sorted.fq* to directory `outdirectory`.

**Note**: users may need to edit the code by using their own path.


```{code-block} python3
referenceGenome_name = "chr1"
referenceGenome_dir = "/home/users/example/refgenome_dir" # may use users' own path
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)

# RNA modality
RNA_synthetic_fastq_prename = RNA_read_bedfile_prename
# Convert combined bed file into FASTQ files
scRNA_GenerateBAM.scRNA_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=RNA_read_bedfile_prename, synthetic_fastq_prename=RNA_synthetic_fastq_prename)

# ATAC modality
ATAC_synthetic_fastq_prename = ATAC_read_bedfile_prename
# Convert combined bed file into FASTQ files
scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=ATAC_read_bedfile_prename, synthetic_fastq_prename=ATAC_synthetic_fastq_prename)
```


### Introduce Error to synthetic data 
Use function `scRNA_GenerateBAM` and `scATAC_ErrorBase` to introduce random error to synthetic reads. 

**Build reference genome dictionary (optional)**

Note that before using function `scRNA_GenerateBAM` and `scATAC_ErrorBase`, please create the reference dictionary for the reference genome with function `CreateSequenceDictionary` using software Picard and make sure that the output *.dict* files are within the same directory to *`referenceGenome_name`.fa*. 

**Note**: For this tutorial, no dictionary building is needed since we have built for *chr1.fa* in *reference.genome.chr1.tar.gz*. The following code chunk is using **bash commands**.

```{code-block} console
$ cd /home/users/example/refgenome_dir # may use users' own path
$ java -jar /home/users/picard/build/libs/picard.jar CreateSequenceDictionary \
$       -R chr1.fa \
$       -O chr1.fa.dict
```

**Introduce errors to synthetic reads**

Functions `scRNA_GenerateBAM` and `scATAC_ErrorBase` take the following arguments:
- `fgbio_jarfile`: Path to software fgbio jar script.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Specify the output directory of the synthteic FASTQ file with random errors.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scRNA_BED2FASTQ` or `scATAC_BED2FASTQ`.

This function will output synthetic reads with random errors in FASTQ files named as *`synthetic_fastq_prename`.ErrorIncluded.read1.bed2fa.fq*, *`synthetic_fastq_prename`.ErrorIncluded.read2.bed2fa.fq* to directory `outdirectory`.


```{code-block} python3
# RNA modality
scRNA_GenerateBAM.scRNA_ErrorBase(fgbio_jarfile=fgbio_jarfile, INPUT_bamfile=INPUT_RNA_bamfile, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, synthetic_fastq_prename=RNA_synthetic_fastq_prename)

# ATAC modality
scATAC_GenerateBAM.scATAC_ErrorBase(fgbio_jarfile=fgbio_jarfile, INPUT_bamfile=INPUT_ATAC_bamfile, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, synthetic_fastq_prename=ATAC_synthetic_fastq_prename)
```



