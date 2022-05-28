# scReadSim on 10x scATAC-seq 

Import modules.

```{code-block} python3
import sys, os
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scATAC_GenerateBAM as scATAC_GenerateBAM
import pkg_resources
```


## Step 1: Download test sample
The example deploys scReadSim on the [10x single cell ATAC-seq](https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0) dataset. The demo BAM file and its corresponding cell barcode file could be accessed through the following chunk. This BAM file uses mm10 as reference genome, the required chromosome size file is also embedded within the package.  

```{code-block} python3
INPUT_cells_barcode_file = pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv') 
filename = "10X_ATAC_chr1_4194444_4399104"
INPUT_bamfile = pkg_resources.resource_filename("scReadSim", 'data/%s.bam' % filename)
INPUT_genome_size_file = pkg_resources.resource_filename("scReadSim", 'data/mm10.chrom.sizes')
```

Other required files for this example inlcuding the reference genome FASTA file and annotation gtf file are downloadable through the following chunk.  

```{code-block} bash
wget http://compbio10data.stat.ucla.edu/repository/gayan/Projects/scReadSim/reference.genmoe.tar.gz # 292 MB
wget http://compbio10data.stat.ucla.edu/repository/gayan/Projects/scReadSim/gencode.vM10.annotation.gtf # 765 MB
```


## Step 2: Feature space construction
For scATAC-seq, scReadSim uses the chromatin open regions (peaks) identified by [MACS3](https://github.com/macs3-project/MACS) as the feature space. Specifically, scReasSim takes peak regions as foreground features and the copmlementary regions along the reference genome as the background features. 

### Specify input parameters 
Create output directory if it does not exist.

```{code-block} python3
outdirectory = "example/outputs"
os.mkdir(outdirectory)

MACS3_peakname_pre = filename + ".MACS3"
ref_peakfile = "%s_peaks.bed" % MACS3_peakname_pre
ref_comple_peakfile = "%s_peaks.COMPLE.bed" % MACS3_peakname_pre
```

### Peak calling
To identify chromatin open regions for scATAC-seq, scReadSim utilizes MACS3 through function `CallPeak` with following arguments
- `macs3_directory`: Path to software MACS3.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Output directory of peak calling.
- `MACS3_peakname_pre`: Base name of peak calling results for MACS3.

The peak calling results by MACS3 would be output into directory `outdirectory`.

```{code-block} python3
Utility.CallPeak(macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre)
```

### Generate feature sets 
Use function `scATAC_CreateFeatureSets` to generate features. This function needs user to specify

- `INPUT_bamfile`: Input BAM file for anlaysis.
- `samtools_directory`: Path to software *samtools*.
- `bedtools_directory`: Path to software *bedtools*.
- `outdirectory`: Specify the output directory of the features files.
- `genome_size_file`: Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicating the size.
- `ref_peakfile`: Specify the name of output foreground feature bed file.
- `ref_comple_peakfile`: Specify the name of output background feature bed file.    
- `MACS3_peakname_pre`: Base name of peak calling results for MACS3.

```{code-block} python3
######################## Generate Feature Set ######################## 
Utility.scATAC_CreateFeatureSets(INPUT_bamfile=INPUT_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_size_file=INPUT_genome_size_file, ref_peakfile=ref_peakfile, ref_comple_peakfile=ref_comple_peakfile, MACS3_peakname_pre=MACS3_peakname_pre)
```


## Step 3: Count matrix construction
Based on the feature sets output in **Step 2**, scReasSim constructs the count matrices for both foreground feautures and background features through function `bam2countmat`. This function needs user to specify

- `cells_barcode_file`: Cell barcode file corresponding to the input BAM file.
- `bed_file`: Features bed file to generate the count matrix.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Specify the output directory of the count matrix file.
- `count_mat_filename`: Specify the base name of output count matrix.

For the user specified `count_mat_filename`, scReadSim will generate a count matrix named `count_mat_filename`.txt to directory `outdirectory`.

```{code-block} python3
count_mat_filename = "%s.countmatrix" % filename
count_mat_comple_filename = "%s.COMPLE.countmatrix" % filename

# Construct count matrix for foregroud features
Utility.bam2countmat(cells_barcode_file=INPUT_cells_barcode_file, bed_file=outdirectory + "/" + ref_peakfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_filename)
# Construct count matrix for background features
Utility.bam2countmat(cells_barcode_file=INPUT_cells_barcode_file, bed_file=outdirectory + "/" + ref_comple_peakfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_comple_filename)
```


## Step 4: Synthetic count matrix simulation
The current version of scReadSim implement scDesign2 (reference) to generate synthetic count matrix based on the constructed count matrix from the input BAM file. Use function `GenerateSyntheticCount.scATAC_GenerateSyntheticCount` to generate synthetic count matrix with following paramters

- `count_mat_filename`: Base name of the count matrix output by function bam2countmat().
- `directory`: Path to the count matrix.
- `outdirectory`: Output directory of coordinate files.
- `cluster_prestep`: Set `cluster_prestep=True` to perform a Louvain clustering before implementing scDesign2.
- `n_cell_new`: (Optional, default: 'None') Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
- `total_count_new`: (Optional, default: 'None') Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
- `celllabel_file`: (Optional, default: 'None') Specify the file containing the predefined cell labels. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.

Given the input count matrix `count_mat_filename`.txt, scReadSim generates two files to `outdirectory` for following analysis:

- **`count_mat_filename`.scDesign2Simulated.txt**: Synthetic count matrix.
- **`count_mat_filename`.scDesign2Simulated.nReadRegionmargional.txt**: The per-feature summation of counts for synthetic count matrix.

```{code-block} python3
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_filename, directory=outdirectory, outdirectory=outdirectory)
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_comple_filename, directory=outdirectory, outdirectory=outdirectory)
```


## Step 5: Synthetic BAM file generation

### Generate synthetic reads in BED format
Based on the synthetic count matrix, scReadSim generates synthetic reads by randomly sampling from the real BAM file input by users. First use function `scATAC_GenerateBAMCoord` to create the synthetic reads and output in BED file storing the coordinates information. Function `scATAC_GenerateBAMCoord` takes following input arguments:
- `count_mat_filename`: The base name of output count matrix in bam2countmat.
- `samtools_directory`: Path to software samtools.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `ref_peakfile`: Features bed file.
- `directory_cellnumber`: Directory of the marginal synthetic count vector file output in scATAC_GenerateSyntheticCount step.
- `outdirectory`: Specify the output directory for synthetic reads bed file.
- `BED_filename`: Specify the base name of output bed file.
- `OUTPUT_cells_barcode_file`: Specify the file name storing the synthetic cell barcodes.
- `read_len`: (Optional, default: '50') Specify the length of synthetic reads. Default value is 50 bp.
- `jitter_size`: (Optional, default: '5') Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.

This function will output a bed file `BED_filename`.bed storing the coordinates information of synthetic reads and its cell barcode file `OUTPUT_cells_barcode_file` in directory `outdirectory`.

```{code-block} python3
directory_cellnumber = outdirectory
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
BED_filename_pre = "%s.syntheticBAM.CBincluded" % filename
BED_COMPLE_filename_pre = "%s.syntheticBAM.COMPLE.CBincluded" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined.CBincluded" % filename

## Create synthetic read coordinates
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(
	count_mat_filename=count_mat_filename, samtools_directory=samtools_directory, INPUT_bamfile=INPUT_bamfile, ref_peakfile=outdirectory + "/" + ref_peakfile, directory_cellnumber=directory_cellnumber, outdirectory=outdirectory, BED_filename=BED_filename_pre, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file)
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(
	count_mat_filename=count_mat_comple_filename, samtools_directory=samtools_directory, INPUT_bamfile=INPUT_bamfile, ref_peakfile=outdirectory + "/" + ref_comple_peakfile, directory_cellnumber=directory_cellnumber, outdirectory=outdirectory, BED_filename=BED_COMPLE_filename_pre, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file)

# Combine foreground and background bed file
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory, BED_filename_pre, BED_COMPLE_filename_pre, BED_filename_combined_pre)
```

### Convert BED files to FASTQ files
Use function `scATAC_BED2FASTQ` to convert BED file to FASTQ file. This function takes the following arguments:
- `bedtools_directory`: Path to software bedtools.
- `seqtk_directory`: Path to software seqtk.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Output directory of the synthteic bed file and its corresponding cell barcodes file.
- `BED_filename_combined`: Base name of the combined bed file output by function `scATAC_CombineBED`.
- `synthetic_fastq_prename`: Specify the base name of the output FASTQ files.
- `sort_FASTQ`: (Optional, default: True) Set `True` to sort the output FASTQ file.
This function will output paired-end reads in FASTQ files named as `synthetic_fastq_prename`.read1.bed2fa.fq, `synthetic_fastq_prename`.read2.bed2fa.fq to directory `outdirectory`.

```{code-block} python3
referenceGenome_name = "chr1"
referenceGenome_dir = "~/Projects/scATAC_Simulator/package_development/package_data" 
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)
synthetic_fastq_prename = BED_filename_combined_pre
output_BAM_pre = "%s.syntheticBAM.CBincluded" % filename

scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=BED_filename_combined_pre, synthetic_fastq_prename=synthetic_fastq_prename, sort_FASTQ = True)
```

### Convert FASTQ files to BAM file
Use function `AlignSyntheticBam_Pair` to align FASTQ files onto reference genome. It takes the following arguments:
- `bowtie2_directory`: Path to software bowtie2.
- `samtools_directory`: Path to software samtools.
- `outdirectory`: Specify the output directory of the synthteic BAM file.
- `referenceGenome_name`: Base name of the reference genome FASTA file. For example, you should input "chr1" for file "chr1.fa".
- `referenceGenome_dir`: Path to the reference genome FASTA file.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scATAC_BED2FASTQ`.
- `output_BAM_pre`: Specify the base name of the output BAM file.

**Important** Note that before using function `AlignSyntheticBam_Pair`, the reference gemome FASTA file should be indexed by bowtie2 through `bowtie2-build ${referenceGenome_name}.fa referenceGenome_name` command and make sure the output index files are within the same directory to `referenceGenome_name`.fa.

```{code-block} python3
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename, output_BAM_pre=output_BAM_pre)
```

### Introduce Error to synthetic data
Use function `scATAC_ErrorBase` to introduce random error to synthetic reads. It takes the following arguments:
- `fgbio_jarfile`: Path to software fgbio jar script.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Specify the output directory of the synthteic FASTQ file with random errors.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scATAC_BED2FASTQ`.
This function will output synthetic reads with random errors in FASTQ files named as `synthetic_fastq_prename`.ErrorIncluded.read1.bed2fa.fq, `synthetic_fastq_prename`.ErrorIncluded.read2.bed2fa.fq to directory `outdirectory`.

**Important** Note that before using function `scATAC_ErrorBase`, please create the reference dictionary with function `CreateSequenceDictionary` using software Picard and make sure the output .dict files are within the same directory to `referenceGenome_name`.fa.

```{code-block} python3
scATAC_GenerateBAM.scATAC_ErrorBase(fgbio_jarfile=fgbio_jarfile, INPUT_bamfile=INPUT_bamfile, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, synthetic_fastq_prename=synthetic_fastq_prename)
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename + ".ErrorIncluded" , output_BAM_pre=output_BAM_pre+ ".ErrorIncluded")
```


