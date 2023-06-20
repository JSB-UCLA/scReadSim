# scReadSim on 10x scATAC-seq with user-designed open chromatin regions

This tutorial's main steps and corresponding estimated time usage are as follows (tested on a server with the 256x Intel Xeon Phi CPU 7210 at 1.30 GHz):

<!-- - [Step 1: Import packages and data files](#step-1-import-packages-and-data-files): < 1 min
- [Step 2: Generate features](#step-2-generate-features): < 1 min
- [Step 3: Generate real count matrices](#step-3-generate-real-count-matrices): < 1 min
- [Step 4: Simulate synthetic count matrix](#step-4-simulate-synthetic-count-matrix): ~ 6 mins
- [Step 5: Output synthetic read](#step-5-output-synthetic-read): ~ 2 mins -->
- **Step 1: Import packages and data files**: < 1 min
- **Step 2: Generate features**: < 1 min
- **Step 3: Generate real count matrices**: < 1 min
- **Step 4: Simulate synthetic count matrix**: ~ 6 mins
- **Step 5: Output synthetic read**: ~ 2 mins

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
**Note: This tutorial does not need this pre-process step since the processed BAM file is provided by the scReadSim package (see Step 1: Import packages and data files).**

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
The example deploys scReadSim on the [10x single cell ATAC-seq](https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0) dataset. For user convienience, we prepared the indexed reference genome files (by bowtie2), which can be downloaded using the following bash commands:
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

The following data files can be accessed by simply loading the code chunk below:
-  BAM file: 10X_ATAC_chr1_4194444_4399104.bam
-  cell barcode file: barcodes.tsv
-  chromosome size file: mm10.chrom.sizes

```{code-block} python3
INPUT_cells_barcode_file = pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv') 
filename = "10X_ATAC_chr1_4194444_4399104"
INPUT_bamfile = pkg_resources.resource_filename("scReadSim", 'data/%s.bam' % filename)
INPUT_genome_size_file = pkg_resources.resource_filename("scReadSim", 'data/mm10.chrom.sizes')
```



## Step 2: Generate features
scReadSim allows users to specify open chromatin regions (referred to as "output peaks") and then generate synthetic scATAC-seq reads accordingly. When users take this option, scReadSim requires users to input the BAM file, trustworthy peaks and non-peaks (i.e., input peaks and input non-peaks. Alternatively, if users do not specify input peaks and non-peaks, scReadSim by default uses [MACS3](https://github.com/macs3-project/MACS) with stringent criteria to call trustworthy peaks (q-value `0.01`) and non-peaks (q-value `0.1`) from the input BAM file) of the BAM file, and a list of output peaks. Given the specified output peaks, scReadSim takes the inter-output-peak as the output non-peaks. In summary, scReadSim defines two sets of peaks and non-peaks: the "input peak and input non-peak" set based on the user-specified (or scReadSim-generated) input peaks and input non-peak and the "output peak and output non-peak" set based on the user-specified output peaks. The following chunks show how to prepare the features for scReadSim with user-spepcified open chromatin regions when anlayzing scATAC-seq data.

### Specify output directory
Specify the absolute path of output directory. Create output directory if it does not exist.

```{code-block} python3
outdirectory = "/home/users/example/outputs" # use absolute path
os.mkdir(outdirectory)
```

### Prepare Features

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

Under default mode "macs3" (by setting argument `peak_mode` as the default values "macs3"), scReadSim uses [MACS3](https://github.com/macs3-project/MACS) with the stringent criteria to call trustworthy peaks (q-value `0.01`) and non-peaks (q-value `0.1`) from the input BAM file. This function will generate the following three bed files into directory `outdirectory` for following analysis:

- peak bed file: *scReadSim.MACS3.peak.bed*
- non-peak bed file: *scReadSim.MACS3.nonpeak.bed*
- gray area bed file: *scReadSim.grayareas.bed*

```{code-block} python3
# Specify the path to OUTPUT_peakfile
OUTPUT_peakfile = pkg_resources.resource_filename("scReadSim", 'data/10X_ATAC_chr1_4194444_4399104.output.peak.bed') 

# Mode: macs3
Utility.scATAC_CreateFeatureSets(peak_mode="macs3", INPUT_bamfile=INPUT_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_size_file=INPUT_genome_size_file, macs3_directory=macs3_directory, INPUT_peakfile=None, INPUT_nonpeakfile=None, OUTPUT_peakfile=OUTPUT_peakfile)
```





## Step 3: Generate real count matrices
To consturct count matrix for user-specified output peaks and non-peaks, scReadSim first uses function `Utility.FeatureMapping` to define the mappings between output peak and input peak, and output non-peak and input non-peak. 


Function `Utility.FeatureMapping` takes following input arguments:

- `INPUT_bamfile`: Input BAM file for anlaysis.
- `input_peaks`: BED file of user specified (or generated by scReadSim+MACS3) input peaks.
- `input_nonpeaks`: BED file of user specified (or generated by scReadSim+MACS3) input non-peaks.
- `output_peaks`: BED file of user specified output peaks.
- `output_nonpeaks`: BED file of user specified output non-peaks.
- `outdirectory`: Output directory of the features assingment file.
- `assignment_peak_file`: Specify the name of peak mapping file.
- `assignment_nonpeak_file`: Specify the name of non-peak mapping file.
- `n_top`: Specify the number of input peaks (or non-peaks) with the most similar length as the candidate mapped input peaks (or non-peaks) for each the output peak (or non-peak). From the candidate input peaks (or non-peaks), scReadSim further selects the one with largest read density for peak mapping (smallest read density for non-peak mapping).

Two mapping files `assignment_file` and `assignment_nonpeak_file` will be output into directory `outdirectory`.

```{code-block} python3
# Specify the path to bed files generated by Utility.scATAC_CreateFeatureSets
input_peaks = outdirectory + "/" + "scReadSim.MACS3.peak.bed"
input_nonpeaks = outdirectory + "/" + "scReadSim.MACS3.nonpeak.bed"
output_peaks = outdirectory + "/" + "scReadSim.output.peak.bed"
output_nonpeaks = outdirectory + "/" + "scReadSim.output.nonpeak.bed"
# Specify the names of peak mapping files
assignment_peak_file = filename + ".assigned.peaks.txt"
assignment_nonpeak_file = filename + ".assigned.nonpeaks.txt"

# Generate mappings for peaks and nonpeaks
Utility.FeatureMapping(INPUT_bamfile=INPUT_bamfile, input_peaks=input_peaks, input_nonpeaks=input_nonpeaks, output_peaks=output_peaks, output_nonpeaks=output_nonpeaks, outdirectory=outdirectory, assignment_peak_file=assignment_peak_file, assignment_nonpeak_file=assignment_nonpeak_file, n_top=50)
```

Based on the bed files output in **Step 2** and the mapping files output by `Utility.FeatureMapping`, scReasSim constructs the count matrices for output peaks and non-peaks through function `Utility.scATAC_bam2countmat_OutputPeak`. This function needs user to specify

- `cells_barcode_file`: Cell barcode file corresponding to the input BAM file.
- `assignment_file`: Features assignment file output by function `Utility.FeatureMapping`.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Specify the output directory of the count matrix file.
- `count_mat_filename`: Specify the base name of output count matrix.

For the user specified `count_mat_filename`, scReadSim will generate a count matrix named *`count_mat_filename`.txt* to directory `outdirectory`.

```{code-block} python3
# Specify the output count matrices' base names
count_mat_peak_filename = "%s.output.peak.countmatrix" % filename
count_mat_nonpeak_filename = "%s.output.nonpeak.countmatrix" % filename

# Construct count matrix for peaks
Utility.scATAC_bam2countmat_OutputPeak(cells_barcode_file=INPUT_cells_barcode_file, assignment_file=outdirectory + "/" + assignment_peak_file, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_peak_filename)
# Construct count matrix for non-peaks
Utility.scATAC_bam2countmat_OutputPeak(cells_barcode_file=INPUT_cells_barcode_file, assignment_file=outdirectory + "/" + assignment_nonpeak_file, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_nonpeak_filename)
```


## Step 4: Simulate synthetic count matrix


In this tutorial, scReadSim implements [scDesign2](https://github.com/JSB-UCLA/scDesign2) to generate synthetic count matrix based on the constructed count matrix from the input BAM file. Use function `GenerateSyntheticCount.scATAC_GenerateSyntheticCount` to generate synthetic count matrix with following paramters

- `count_mat_filename`: Base name of the count matrix output by function `Utility.scATAC_bam2countmat_paral`.
- `directory`: Path to the count matrix.
- `outdirectory`: Specify the output directory of the synthetic count matrix file.
- `doub_classification_label_file`: (Optional, default: 'None') Specify the absolute path to the doublet classification result `doublet_classification.Rdata` generated by function `DoubletDetection.detectDoublet`.
- `n_cell_new`: (Optional, default: 'None') Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
- `total_count_new`: (Optional, default: 'None') Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
- `celllabel_file`: (Optional, default: 'None') Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.
- `n_cores`: (Optional, default: '1') Specify the number of cores for parallel computing when generating count matrix.


Given the input count matrix *`count_mat_filename`.txt*, scReadSim generates the syntheitic count matrix file to `outdirectory` for following analysis:

- Synthetic count matrix: *`count_mat_filename`.scDesign2Simulated.txt*
- Synthetic cell cluster/type labels: *`count_mat_filename`.scDesign2Simulated.CellTypeLabel.txt*

Additionaly, if no `celllabel_file` is specified, scReadSim automatically performs Louvain clustering from Seurat and outputs clustering labels to `outdirectory`:
- Real cells' Louvain clustering labels: *`count_mat_filename`.LouvainClusterResults.txt*


```{code-block} python3
# Generate synthetic count matrix for peak-by-cell count matrix
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_peak_filename, directory=outdirectory, outdirectory=outdirectory)

# Specify cluster labels obtained from peak-by-cell matrix
celllabel_file = outdirectory + "/" + "10X_ATAC_chr1_4194444_4399104.output.peak.countmatrix.LouvainClusterResults.txt"
# Generate synthetic count matrix for nonpeak-by-cell count matrix
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_nonpeak_filename, directory=outdirectory, outdirectory=outdirectory, celllabel_file=celllabel_file)
```


## Step 5: Output synthetic read

### Generate synthetic reads in BED format
Based on the synthetic count matrix, scReadSim generates synthetic reads by randomly sampling from the real BAM file input by users. First use function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord_OutputPeak` to create the synthetic reads and output in BED file storing the coordinates information. Function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord_OutputPeak` takes following input arguments:

- `target_peak_assignment_file`: Mapping file between input peaks and output peaks, output by 'FeatureMapping'.
- `count_mat_file`: The path to the **synthetic count matrix** generated by `GenerateSyntheticCount.scATAC_GenerateSyntheticCount`.
- `synthetic_cell_label_file`: Synthetic cell label file generated by `scATAC_GenerateSyntheticCount`.
- `read_bedfile_prename`: Specify the base name of output bed file.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Specify the output directory for synthetic reads bed file.
- `OUTPUT_cells_barcode_file`: Specify the file name storing the synthetic cell barcodes.
- `jitter_size`: (Optional, default: '5') Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
- `read_len`: (Optional, default: '50') Specify the length of synthetic reads. Default value is 50 bp.
- `random_noise_mode`: (Optional, default: 'False') Specify whether to use a uniform distribution of reads.


This function will output two bed files *`read_bedfile_prename`.read1.bed* and *`read_bedfile_prename`.read2.bed* storing the coordinates information of synthetic reads and its cell barcode file `OUTPUT_cells_barcode_file` in directory `outdirectory`.

After generation of synthetic reads for both peaks and non-peaks, combine the their bed files using function `scATAC_GenerateBAM.scATAC_CombineBED`, which takes following input arguments:
- `outdirectory`: Directory of `peak_read_bedfile_prename`.txt and `nonpeak_read_bedfile_prename`.txt.
- `peak_read_bedfile_prename`: Base name of the bed file containig synthetic reads for peaks (generated by function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord`).
- `nonpeak_read_bedfile_prename`: Base name of the bed file containig synthetic reads for non-peaks (generated by function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord`).
- `BED_filename_combined_pre`: Specify the base name for the combined syntehtic reads bed file. The combined bed file will be output to `outdirectory`.


```{code-block} python3
# Specify the names of synthetic count matrices (generated by GenerateSyntheticCount.scATAC_GenerateSyntheticCount)
synthetic_countmat_peak_file = count_mat_peak_filename + ".scDesign2Simulated.txt"
synthetic_countmat_nonpeak_file = count_mat_nonpeak_filename + ".scDesign2Simulated.txt"
# Specify the base name of bed files containing synthetic reads
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
peak_read_bedfile_prename = "%s.syntheticBAM.peak" % filename
nonpeak_read_bedfile_prename = "%s.syntheticBAM.nonpeak" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined" % filename
synthetic_cell_label_file = count_mat_peak_filename + ".scDesign2Simulated.CellTypeLabel.txt"

# Create synthetic read bed file for peaks
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_OutputPeak(target_peak_assignment_file=outdirectory + "/"+ assignment_peak_file, count_mat_file=outdirectory + "/" + synthetic_countmat_peak_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=peak_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50)

# Create synthetic read bed file for non-peaks
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_OutputPeak(target_peak_assignment_file=outdirectory + "/"+ assignment_nonpeak_file, count_mat_file=outdirectory + "/" + synthetic_countmat_nonpeak_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=nonpeak_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50)

# Combine bed files
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory=outdirectory, peak_read_bedfile_prename=peak_read_bedfile_prename, nonpeak_read_bedfile_prename=nonpeak_read_bedfile_prename, BED_filename_combined_pre=BED_filename_combined_pre, GrayAreaModeling=False)
```

### Convert BED files to FASTQ files
Use function `scATAC_BED2FASTQ` to convert BED file to FASTQ file. This function takes the following arguments:
- `bedtools_directory`: Path to software bedtools.
- `seqtk_directory`: Path to software seqtk.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Output directory of the synthteic bed file and its corresponding cell barcodes file.
- `BED_filename_combined`: Base name of the combined bed file output by function `scATAC_CombineBED`.
- `synthetic_fastq_prename`: Specify the base name of the output FASTQ files.

This function will output paired-end reads in FASTQ files named as *`synthetic_fastq_prename`.read1.bed2fa.sorted.fq*, *`synthetic_fastq_prename`.read2.bed2fa.sorted.fq* to directory `outdirectory`.

**Note**: users may need to edit the code by using their own path.


```{code-block} python3
referenceGenome_name = "chr1"
referenceGenome_dir = "/home/users/example/refgenome_dir" # may use users' own path
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)
synthetic_fastq_prename = BED_filename_combined_pre

# Convert combined bed file into FASTQ files
scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=BED_filename_combined_pre, synthetic_fastq_prename=synthetic_fastq_prename)
```



### Introduce Error to synthetic data 
Use function `scATAC_ErrorBase` to introduce random error to synthetic reads. 

**Build reference genome dictionary (optional)**

Note that before using function `scATAC_ErrorBase`, please create the reference dictionary for the reference genome with function `CreateSequenceDictionary` using software Picard and make sure that the output *.dict* files are within the same directory to *`referenceGenome_name`.fa*. **For this tutorial, no dictionary building is needed since we have built for chr1.fa in reference.genome.chr1.tar.gz**. 

```{code-block} console
$ cd /home/users/example/refgenome_dir # may use users' own path
$ java -jar /home/users/picard/build/libs/picard.jar CreateSequenceDictionary \
$       -R chr1.fa \
$       -O chr1.fa.dict
```

**Introduce errors to synthetic reads**

Function `scATAC_ErrorBase` takes the following arguments:
- `fgbio_jarfile`: Path to software fgbio jar script.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Specify the output directory of the synthteic FASTQ file with random errors.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scATAC_BED2FASTQ`.

This function will output synthetic reads with random errors in FASTQ files named as *`synthetic_fastq_prename`.ErrorIncluded.read1.bed2fa.fq*, *`synthetic_fastq_prename`.ErrorIncluded.read2.bed2fa.fq* to directory `outdirectory`.


```{code-block} python3
# Generate reads with errors in FASTQs
scATAC_GenerateBAM.scATAC_ErrorBase(fgbio_jarfile=fgbio_jarfile, INPUT_bamfile=INPUT_bamfile, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, synthetic_fastq_prename=synthetic_fastq_prename)
```

### Convert FASTQ files to BAM file (optional)
The current version of scReadSim implicitly uses [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align the synthetic reads onto the reference genome. Use function `AlignSyntheticBam_Pair` to align FASTQ files onto reference genome. It takes the following arguments:
- `bowtie2_directory`: Path to software bowtie2.
- `samtools_directory`: Path to software samtools.
- `outdirectory`: Specify the output directory of the synthteic BAM file.
- `referenceGenome_name`: Base name of the reference genome FASTA file. For example, you should input "chr1" for file "chr1.fa".
- `referenceGenome_dir`: Path to the reference genome FASTA file.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scATAC_BED2FASTQ`.
- `output_BAM_pre`: Specify the base name of the output BAM file.

**Index reference genome (optional)** 

Before using function `AlignSyntheticBam_Pair`, the reference gemome FASTA file should be indexed by bowtie2 through following chunk and make sure the output index files are within the same directory to *`referenceGenome_name`.fa*. **For this tutorial, no indexing is needed since we have indexed chr1.fa in reference.genome.chr1.tar.gz**. 


```{code-block} console
$ cd /home/users/example/refgenome_dir # may use users' own path
$ bowtie2-build chr1.fa chr1
```

**Align synthetic reads** 

Now align the synthetic reads on to the reference genome with bowtie2.

```{code-block} python3
output_BAM_pre = "%s.syntheticBAM.CBincluded" % filename

# Synthetic reads alignment
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename, output_BAM_pre=output_BAM_pre)

# Synthetic reads (with sequencing errors) alignment
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename + ".ErrorIncluded" , output_BAM_pre=output_BAM_pre+ ".ErrorIncluded")
```
