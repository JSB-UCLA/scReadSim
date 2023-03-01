# scReadSim on 10x scATAC-seq with user-input chromatin regions

Import modules.

```{code-block} python3
import sys, os
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scATAC_GenerateBAM as scATAC_GenerateBAM
import pkg_resources
```

## Required softwares for scReadSim
scReadSim requires users to pre-install the following softwares:
- [MACS3](https://github.com/macs3-project/MACS)
- [samtools](http://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [seqtk](https://github.com/lh3/seqtk)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [fgbio](http://fulcrumgenomics.github.io/fgbio/)


## Step 1: Download test sample
The example deploys scReadSim on the [10x single cell ATAC-seq](https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0) dataset. The demo BAM file and its corresponding cell barcode file could be accessed through the following chunk. This BAM file uses mm10 as reference genome, the required chromosome size file is also embedded within the package.  

```{code-block} python3
INPUT_cells_barcode_file = pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv') 
filename = "10X_ATAC_chr1_4194444_4399104"
INPUT_bamfile = pkg_resources.resource_filename("scReadSim", 'data/%s.bam' % filename)
INPUT_genome_size_file = pkg_resources.resource_filename("scReadSim", 'data/mm10.chrom.sizes')
```

Use the following chunk to download other required files for this example, inlcuding the reference genome FASTA file (indexed by bowtie2) and annotation gtf file. 

```{code-block} console
$ mkdir example/refgenome_dir 
$ cd example/refgenome_dir
$ wget http://compbio10data.stat.ucla.edu/repository/gayan/Projects/scReadSim/reference.genome.chr1.tar.gz # 292 MB
$ wget http://compbio10data.stat.ucla.edu/repository/gayan/Projects/scReadSim/gencode.vM10.annotation.gtf # 765 MB
$ tar -xf reference.genome.chr1.tar.gz
```

### Pre-process input BAM file
Note: Input BAM file for scReadSim needs pre-processing to add the cell barcode in front of the read name. For example, in 10x sequencing data, cell barcode `TGGACCGGTTCACCCA-1` is stored in the field `CB:Z:TGGACCGGTTCACCCA-1`. 

```{code-block} console
$ samtools view 10X_ATAC_chr1_4194444_4399104_unprocess.bam | head -n 1
A00836:472:HTNW5DMXX:1:1372:16260:18129      83      chr1    4194410 60      50M     =       4193976 -484    TGCCTTGCTACAGCAGCTCAGGAAATGTCTTTGTGCCCACAGTCTGTGGT   :FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF      NM:i:0  MD:Z:50 AS:i:50 XS:i:0  CR:Z:TCCGGGACAGCTAACA   CY:Z:FFFFFFFFFFFFFFF:   CB:Z:TGGACCGGTTCACCCA-1 BC:Z:AAACTCAT        QT:Z::FFFFFFF   RG:Z:e18_mouse_brain_fresh_5k:MissingLibrary:1:HTNW5DMXX:1
```

The following code chunk adds the cell barcodes in front of the read names.

```{code-block} console
$ # extract the header file
$ mkdir tmp
$ samtools view 10X_ATAC_chr1_4194444_4399104_unprocess.bam -H > tmp/10X_ATAC_chr1_4194444_4399104.header.sam

$ # create a bam file with the barcode embedded into the read name
$ time(cat <( cat tmp/10X_ATAC_chr1_4194444_4399104.header.sam ) \
 <( samtools view 10X_ATAC_chr1_4194444_4399104_unprocess.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
 | samtools view -bS - > 10X_ATAC_chr1_4194444_4399104.bam) 
$ rm -dr tmp

$ samtools view 10X_ATAC_chr1_4194444_4399104.bam | head -n 1
TGGACCGGTTCACCCA-1:A00836:472:HTNW5DMXX:1:1372:16260:18129      83      chr1    4194410 60      50M     =       4193976 -484    TGCCTTGCTACAGCAGCTCAGGAAATGTCTTTGTGCCCACAGTCTGTGGT   :FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF      NM:i:0  MD:Z:50 AS:i:50 XS:i:0  CR:Z:TCCGGGACAGCTAACA   CY:Z:FFFFFFFFFFFFFFF:   CB:Z:TGGACCGGTTCACCCA-1 BC:Z:AAACTCAT        QT:Z::FFFFFFF   RG:Z:e18_mouse_brain_fresh_5k:MissingLibrary:1:HTNW5DMXX:1
```

## Step 2: Feature space construction
scReadSim allows users to specify open chromatin regions (referred to as "output peaks") and then generate synthetic scATAC-seq reads accordingly. When users take this option, scReadSim requires users to input the BAM file, trustworthy peaks and non-peaks (i.e., input peaks and input non-peaks. Alternatively, if users do not specify input peaks and non-peaks, scReadSim by default uses [MACS3](https://github.com/macs3-project/MACS) with stringent criteria to call trustworthy peaks (q-value `0.01`) and non-peaks (q-value `0.1`) from the input BAM file) of the BAM file, and a list of output peaks. Given the specified output peaks, scReadSim takes the inter-output-peak as the output non-peaks. In summary, scReadSim defines two sets of peaks and non-peaks: the "input peak and input non-peak" set based on the user-specified (or scReadSim-generated) input peaks and input non-peak and the "output peak and output non-peak" set based on the user-specified output peaks. The following chunks show how to prepare the features for scReadSim with user-spepcified open chromatin regions when anlayzing scATAC-seq data.

### Specify output directory
Specify the absolute path of output directory. Create output directory if it does not exist.

```{code-block} python3
outdirectory = "/home/users/example/outputs" # use absolute path
os.mkdir(outdirectory)
```

### Prepare features
To generate bed files for input peak and input non-peak, and output peak and output non-peak, use function `Utility.scATAC_CreateFeatureSets` with following arguments:
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `samtools_directory`: Directory of software samtools.
- `bedtools_directory`: Directory of software bedtools.
- `outdirectory`: Output directory of the prepared features.
- `genome_size_file`: Directory of Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicates the size.  
- `macs3_directory`: Path to software MACS3.
- `INPUT_peakfile`: (Optional, default: None) Directory of user-specified input peak file.
- `INPUT_nonpeakfile`: (Optional, default: None) Directory of user-specified input non-peak file.
- `OUTPUT_peakfile`: (Optional, default: None) Directory of user-specified output peak file. Synthetic scATAC-seq reads will be generated taking `OUTPUT_peakfile` as ground truth peaks. Note that `OUTPUT_peakfile` does not name the generated feature files by function `scATAC_CreateFeatureSets`.

#### **Case 1** without user-specified input peaks and input non-peaks
If users do not specify peaks and non-peaks (by setting options `INPUT_peakfile` and `INPUT_nonpeakfile` as the default values None) for the input BAM file, scReadSim by default uses MACS3 to determine the peaks and non-peaks. This function will generate the following three bed files into directory `outdirectory` for following analysis:

- peak bed file: *scReadSim.MACS3.peak.bed*
- non-peak bed file: *scReadSim.MACS3.nonpeak.bed*
- gray area bed file: *scReadSim.grayareas.bed*

```{code-block} python3
# Specify the path to output peak bed file
OUTPUT_peakfile = pkg_resources.resource_filename("scReadSim", 'data/10x_ATAC_chr1_4194444_4399104.output.peaks.bed') 

# Prepare features without user-specified peaks and non-peaks
Utility.scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, INPUT_genome_size_file, macs3_directory, INPUT_peakfile=None, INPUT_nonpeakfile=None, OUTPUT_peakfile=OUTPUT_peakfile)
```


#### **Case 2** with user-specified input peaks and input non-peaks
If users specify peaks and non-peaks (by setting options `INPUT_peakfile` and `INPUT_nonpeakfile` as the path to the ground truth peak and non-peak bed files) for the input BAM file, scReadSim further preprocesses the bed files and generate the following three bed files into directory `outdirectory` for following analysis:

- peak bed file: *scReadSim.UserInput.peak.bed*
- non-peak bed file: *scReadSim.UserInput.nonpeak.bed*
- gray area bed file: *scReadSim.grayareas.bed*

```{code-block} python3
# Specify the path to peak and non-peak bed files
INPUT_peakfile = pkg_resources.resource_filename("scReadSim", 'data/10x_ATAC_chr1_4194444_4399104.input.peak.bed')
INPUT_nonpeakfile = pkg_resources.resource_filename("scReadSim", 'data/10x_ATAC_chr1_4194444_4399104.input.nonpeak.bed')
# Specify the path to output peak bed file
OUTPUT_peakfile = pkg_resources.resource_filename("scReadSim", 'data/10x_ATAC_chr1_4194444_4399104.output.peaks.bed') 

# Prepare features with user-specified peaks and non-peaks
Utility.scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, INPUT_genome_size_file, macs3_directory, INPUT_peakfile, INPUT_nonpeakfile, OUTPUT_peakfile=OUTPUT_peakfile)
```


## Step 3: Count matrix construction
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


## Step 4: Synthetic count matrix simulation
The current version of scReadSim implements [scDesign2](https://github.com/JSB-UCLA/scDesign2) to generate synthetic count matrix based on the constructed count matrix from the input BAM file. Use function `GenerateSyntheticCount.scATAC_GenerateSyntheticCount` to generate synthetic count matrix with following paramters

- `bed_file`: Features' bed file to generate the count matrix (Generated by function `Utility.scATAC_CreateFeatureSets`).
- `count_mat_filename`: Base name of the count matrix output by function `Utility.scATAC_bam2countmat_paral`.
- `directory`: Path to the count matrix.
- `outdirectory`: Specify the output directory of the synthetic count matrix file.
- `n_cell_new`: (Optional, default: 'None') Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
- `total_count_new`: (Optional, default: 'None') Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
- `celllabel_file`: (Optional, default: 'None') Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.

Given the input count matrix *`count_mat_filename`.txt*, scReadSim generates the syntheitic count matrix file to `outdirectory` for following analysis:

- **`count_mat_filename`.scDesign2Simulated.txt**: Synthetic count matrix.


```{code-block} python3
# Generate synthetic count matrix for peak-by-cell count matrix
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_peak_filename, directory=outdirectory, outdirectory=outdirectory)
# Generate synthetic count matrix for nonpeak-by-cell count matrix
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_nonpeak_filename, directory=outdirectory, outdirectory=outdirectory)
```


## Step 5: Synthetic BAM file generation

### Generate synthetic reads in BED format
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

This function will output two bed files *`read_bedfile_prename:`.read1.bed* and *`read_bedfile_prename:`.read2.bed* storing the coordinates information of synthetic reads and its cell barcode file `OUTPUT_cells_barcode_file` in directory `outdirectory`.

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
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(bed_file=peak_bedfile, count_mat_file=outdirectory + "/" + synthetic_countmat_peak_file, synthetic_cell_label_file=synthetic_cell_label_file, read_bedfile_prename=peak_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50)

# Create synthetic read bed file for non-peaks
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(bed_file=nonpeak_bedfile, count_mat_file=outdirectory + "/" + synthetic_countmat_nonpeak_file, synthetic_cell_label_file=synthetic_cell_label_file, read_bedfile_prename=nonpeak_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50,  GrayAreaModeling=True))

# Combine bed files
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory=outdirectory, peak_read_bedfile_prename=peak_read_bedfile_prename, nonpeak_read_bedfile_prename=nonpeak_read_bedfile_prename, BED_filename_combined_pre=BED_filename_combined_pre)
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

```{code-block} python3
referenceGenome_name = "chr1"
referenceGenome_dir = "/home/users/example/refgenome_dir" 
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)
synthetic_fastq_prename = BED_filename_combined_pre

# Convert combined bed file into FASTQ files
scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=BED_filename_combined_pre, synthetic_fastq_prename=synthetic_fastq_prename)
```

### Convert FASTQ files to BAM file (optional)
Use function `AlignSyntheticBam_Pair` to align FASTQ files onto reference genome. It takes the following arguments:
- `bowtie2_directory`: Path to software bowtie2.
- `samtools_directory`: Path to software samtools.
- `outdirectory`: Specify the output directory of the synthteic BAM file.
- `referenceGenome_name`: Base name of the reference genome FASTA file. For example, you should input "chr1" for file "chr1.fa".
- `referenceGenome_dir`: Path to the reference genome FASTA file.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scATAC_BED2FASTQ`.
- `output_BAM_pre`: Specify the base name of the output BAM file.

> **Important** Note that before using function `AlignSyntheticBam_Pair`, the reference gemome FASTA file should be indexed by bowtie2 through following chunk and make sure the output index files are within the same directory to *`referenceGenome_name`.fa*.

```{code-block} console
$ cd example/refgenome_dir # change to directory where your reference genome file is
$ bowtie2-build chr1.fa chr1
```

In the demo data, We have indexed chr1.fa stored in reference.genome.chr1.tar.gz. Now align the synthetic reads on to the reference genome with bowtie2.

```{code-block} python3
output_BAM_pre = "%s.syntheticBAM.CBincluded" % filename

# Convert FASTQ files to BAM file
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename, output_BAM_pre=output_BAM_pre)
```

### Introduce Error to synthetic data 
Use function `scATAC_ErrorBase` to introduce random error to synthetic reads. It takes the following arguments:
- `fgbio_jarfile`: Path to software fgbio jar script.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Specify the output directory of the synthteic FASTQ file with random errors.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scATAC_BED2FASTQ`.

This function will output synthetic reads with random errors in FASTQ files named as *`synthetic_fastq_prename`.ErrorIncluded.read1.bed2fa.fq*, *`synthetic_fastq_prename`.ErrorIncluded.read2.bed2fa.fq* to directory `outdirectory`.

> **Important** Note that before using function `scATAC_ErrorBase`, please create the reference dictionary for the reference genome with function `CreateSequenceDictionary` using software Picard and make sure that the output *.dict* files are within the same directory to *`referenceGenome_name`.fa*.

```{code-block} console
$ cd example/refgenome_dir # change to directory where your reference genome file is
$ java -jar /home/users/picard/build/libs/picard.jar CreateSequenceDictionary \
$       -R chr1.fa \
$       -O chr1.fa.dict
```

In the demo data, We have built the dictionary file chr1.fa.dict for chr1.fa stored in reference.genome.chr1.tar.gz. 

```{code-block} python3
# Generate reads with errors in FASTQs
scATAC_GenerateBAM.scATAC_ErrorBase(fgbio_jarfile=fgbio_jarfile, INPUT_bamfile=INPUT_bamfile, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, synthetic_fastq_prename=synthetic_fastq_prename)
# Reads alignment (optional)
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename + ".ErrorIncluded" , output_BAM_pre=output_BAM_pre+ ".ErrorIncluded")
```


