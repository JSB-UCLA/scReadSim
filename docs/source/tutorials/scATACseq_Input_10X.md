# scReadSim on 10x scATAC-seq with user-input chromatin regions

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
INPUT_peakfile = pkg_resources.resource_filename("scReadSim", 'data/%s.INPUT.peaks.bed' % filename) 
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
$ rm -d tmp

$ samtools view 10X_ATAC_chr1_4194444_4399104.bam | head -n 1
TGGACCGGTTCACCCA-1:A00836:472:HTNW5DMXX:1:1372:16260:18129      83      chr1    4194410 60      50M     =       4193976 -484    TGCCTTGCTACAGCAGCTCAGGAAATGTCTTTGTGCCCACAGTCTGTGGT   :FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF      NM:i:0  MD:Z:50 AS:i:50 XS:i:0  CR:Z:TCCGGGACAGCTAACA   CY:Z:FFFFFFFFFFFFFFF:   CB:Z:TGGACCGGTTCACCCA-1 BC:Z:AAACTCAT        QT:Z::FFFFFFF   RG:Z:e18_mouse_brain_fresh_5k:MissingLibrary:1:HTNW5DMXX:1
```

## Step 2: Feature space construction
When users specify features for their synthetic scATAC-seq data, scReadSim takes the user input features as the foreground features, and uses the copmlementary regions along the reference genome as the background features. The deom input chromatin open regions used in this example are truncated from the transcription starting sites of genes with length randomly chosen from $[250, 550]$ bp. 

Meanwhile, scReadSim also needs to construct another pair of foreground and background features by using the real data. Specifically, scReadSim uses the chromatin open regions (peaks) identified by [MACS3](https://github.com/macs3-project/MACS) as the feature space. Specifically, scReasSim takes peak regions as foreground features and the copmlementary regions along the reference genome as the background features. 

### Specify input parameters 
Specify the absolute path of output directory. Create output directory if it does not exist.

```{code-block} python3
outdirectory = "/home/users/example/outputs" # use absolute path
os.mkdir(outdirectory)
```

### Prepare user input features
To generate background features for input features, use function `Utility.ComplementFeature` with following arguments:

- `feature_file`: Input feature set.
- `comple_feature_peakfile`: Specify the base name of output background feature bed file.    
- `genome_size_file`: Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicating the size.
- `outdirectory`: Output directory.
- `bedtools_directory`: Path to software bedtools.

```{code-block} python3
input_comple_peakfile = "%s.INPUT.COMPLE.peaks.bed" % filename 

# Generate background features for user input features
Utility.ComplementFeature(feature_file=INPUT_peakfile, comple_feature_peakfile=input_comple_peakfile, genome_size_file=INPUT_genome_size_file, outdirectory=outdirectory, bedtools_directory=bedtools_directory)
```

### Prepare real features

To identify chromatin open regions for scATAC-seq, scReadSim utilizes MACS3 through function `Utility.CallPeak` with following arguments
- `macs3_directory`: Path to software MACS3.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Output directory of peak calling.
- `MACS3_peakname_pre`: Base name of peak calling results for MACS3.

The peak calling results by MACS3 would be output into directory `outdirectory`.

```{code-block} python3
MACS3_peakname_pre = filename + ".MACS3"

# Peak calling
Utility.CallPeak(macs3_directory=macs3_directory, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, MACS3_peakname_pre=MACS3_peakname_pre)
```

Use function `Utility.scATAC_CreateFeatureSets` to generate features. This function needs user to specify

- `INPUT_bamfile`: Input BAM file for anlaysis.
- `samtools_directory`: Path to software *samtools*.
- `bedtools_directory`: Path to software *bedtools*.
- `outdirectory`: Specify the output directory of the features files.
- `genome_size_file`: Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicating the size.
- `ref_peakfile`: Specify the name of output foreground feature bed file.
- `ref_comple_peakfile`: Specify the name of output background feature bed file.    
- `MACS3_peakname_pre`: Base name of peak calling results for MACS3.

```{code-block} python3
ref_peakfile = "%s_peaks.bed" % MACS3_peakname_pre
ref_comple_peakfile = "%s_peaks.COMPLE.bed" % MACS3_peakname_pre

# Generate real features 
Utility.scATAC_CreateFeatureSets(INPUT_bamfile=INPUT_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_size_file=INPUT_genome_size_file, ref_peakfile=ref_peakfile, ref_comple_peakfile=ref_comple_peakfile, MACS3_peakname_pre=MACS3_peakname_pre)
```


## Step 3: Count matrix construction
To consturct count matrix for user-input foreground features and its corresponding background counterparts, scReadSim first uses function `Utility.match_peak` to assign a real foreground(or background) feature to each user-input foreground (or background) feature based on the feature length. 

Function `Utility.match_peak` takes following input arguments:

- `input_peakfile`: User input foreground(or background) features bed file.
- `real_peakfile`: Real foreground(or background) features bed file. 
- `outdirectory`: Output directory of the features assingment file.
- `assignment_file`: Specify the name of features assignment file.

The assignment file `assignment_file` will be output into directory `outdirectory`.

```{code-block} python3
assignment_file = filename + ".assigned.peaks.txt"
assignment_comple_file=filename + ".COMPLE.assigned.peaks.txt"

# Construct mappings for foreground features
Utility.match_peak(true_peakfile=INPUT_peakfile, ref_peakfile=outdirectory + "/" + ref_peakfile, outdirectory=outdirectory, assignment_file=assignment_file)
# Construct mappings for background features
Utility.match_peak(true_peakfile=outdirectory + "/" + input_comple_peakfile, ref_peakfile=outdirectory + "/" + ref_comple_peakfile, outdirectory=outdirectory, assignment_file=assignment_comple_file)
```

Based on the feature sets output in **Step 2**, scReasSim constructs the count matrices for foreground and background features through function `Utility.bam2countmat_INPUT`. This function needs user to specify

- `cells_barcode_file`: Cell barcode file corresponding to the input BAM file.
- `assignment_file`: Features assignment file output by function `match_peak`.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Specify the output directory of the count matrix file.
- `count_mat_filename`: Specify the base name of output count matrix.

For the user specified `count_mat_filename`, scReadSim will generate a count matrix named *`count_mat_filename`.txt* to directory `outdirectory`.

```{code-block} python3
count_mat_filename = "%s.assigned.countmatrix" % filename
count_mat_comple_filename = "%s.assigned.COMPLE.countmatrix" % filename

# Construct count matrix for foregroud features
Utility.bam2countmat_INPUT(cells_barcode_file=INPUT_cells_barcode_file, assignment_file=outdirectory + "/" + assignment_file, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_filename)
# Construct count matrix for background features
Utility.bam2countmat_INPUT(cells_barcode_file=INPUT_cells_barcode_file, assignment_file=outdirectory + "/" + assignment_comple_file, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_comple_filename)
```


## Step 4: Synthetic count matrix simulation
The current version of scReadSim implement [scDesign2](https://github.com/JSB-UCLA/scDesign2) to generate synthetic count matrix based on the constructed count matrix from the input BAM file. Use function `GenerateSyntheticCount.scATAC_GenerateSyntheticCount` to generate synthetic count matrix with following paramters

- `count_mat_filename`: Base name of the count matrix output by function bam2countmat().
- `directory`: Path to the count matrix.
- `outdirectory`: Output directory of coordinate files.
- `n_cell_new`: (Optional, default: 'None') Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
- `total_count_new`: (Optional, default: 'None') Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
- `celllabel_file`: (Optional, default: 'None') Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.


Given the input count matrix *`count_mat_filename`.txt*, scReadSim generates two files to `outdirectory` for following analysis:

- **`count_mat_filename`.scDesign2Simulated.txt**: Synthetic count matrix.
- **`count_mat_filename`.scDesign2Simulated.nReadRegionmargional.txt**: The per-feature summation of counts for synthetic count matrix.

```{code-block} python3
# Generate synthetic count matrix for foregroud features
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_filename, directory=outdirectory, outdirectory=outdirectory)
# Generate synthetic count matrix for backgroud features
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_comple_filename, directory=outdirectory, outdirectory=outdirectory)
```


## Step 5: Synthetic BAM file generation

### Generate synthetic reads in BED format
Based on the synthetic count matrix, scReadSim generates synthetic reads by randomly sampling from the real BAM file input by users. First use function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord_INPUT` to create the synthetic reads and output in BED file storing the coordinates information. Function `scATAC_GenerateBAM.scATAC_GenerateBAMCoord_INPUT` takes following input arguments:
- `count_mat_filename`: The base name of output count matrix in bam2countmat.
- `samtools_directory`: Path to software samtools.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `assignment_file`: Features assignment file output by function `match_peak`.
- `directory_cellnumber`: Directory of the marginal synthetic count vector file output in scATAC_GenerateSyntheticCount step.
- `outdirectory`: Specify the output directory for synthetic reads bed file.
- `BED_filename`: Specify the base name of output bed file.
- `OUTPUT_cells_barcode_file`: Specify the file name storing the synthetic cell barcodes.
- `read_len`: (Optional, default: '50') Specify the length of synthetic reads. Default value is 50 bp.
- `jitter_size`: (Optional, default: '5') Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.

This function will output a bed file *`BED_filename`.bed* storing the coordinates information of synthetic reads and its cell barcode file `OUTPUT_cells_barcode_file` in directory `outdirectory`.

After generation of synthetic reads for both foreground and background features, combine the two bed files using function `scATAC_GenerateBAM.scATAC_CombineBED`, which takes following input arguments:
- `outdirectory`: Directory of `BED_filename_pre`.txt and `BED_COMPLE_filename_pre`.txt.
- `BED_filename_pre`: File prename of foreground synthetic reads bed file.
- `BED_COMPLE_filename_pre`: File prename of background synthetic reads bed file.
- `BED_filename_combined_pre`: Specify the combined syntehtic reads bed file prename. The combined bed file will be output to `outdirectory`.

```{code-block} python3
directory_cellnumber = outdirectory
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
BED_filename_pre = "%s.syntheticBAM" % filename
BED_COMPLE_filename_pre = "%s.syntheticBAM.COMPLE" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined" % filename

# Create synthetic read coordinates for foregroud features
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_INPUT(
	count_mat_filename=count_mat_filename, samtools_directory=samtools_directory, INPUT_bamfile=INPUT_bamfile, assignment_file=outdirectory + "/" + assignment_file, directory_cellnumber=directory_cellnumber, outdirectory=outdirectory, BED_filename=BED_filename_pre, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file)
# Create synthetic read coordinates for backgroud features
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_INPUT(
	count_mat_filename=count_mat_comple_filename, samtools_directory=samtools_directory, INPUT_bamfile=INPUT_bamfile, assignment_file=outdirectory + "/" + assignment_comple_file, directory_cellnumber=directory_cellnumber, outdirectory=outdirectory, BED_filename=BED_COMPLE_filename_pre, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file)

# Combine foreground and background bed file
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory=outdirectory, BED_filename_pre=BED_filename_pre, BED_COMPLE_filename_pre=BED_COMPLE_filename_pre, BED_filename_combined_pre=BED_filename_combined_pre)
```


### Convert BED files to FASTQ files
Use function `scATAC_GenerateBAM.scATAC_BED2FASTQ` to convert BED file to FASTQ file. This function takes the following arguments:
- `bedtools_directory`: Path to software bedtools.
- `seqtk_directory`: Path to software seqtk.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Output directory of the synthteic bed file and its corresponding cell barcodes file.
- `BED_filename_combined`: Base name of the combined bed file output by function `scATAC_CombineBED`.
- `synthetic_fastq_prename`: Specify the base name of the output FASTQ files.
- `sort_FASTQ`: (Optional, default: True) Set `True` to sort the output FASTQ file.

This function will output paired-end reads in FASTQ files named as *`synthetic_fastq_prename`.read1.bed2fa.fq*, *`synthetic_fastq_prename`.read2.bed2fa.fq* to directory `outdirectory`.

```{code-block} python3
referenceGenome_name = "chr1"
referenceGenome_dir = "/home/users/example/refgenome_dir" # Use absolut path
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)

# Convert combined bed file into FASTQ files
scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=BED_filename_combined_pre, synthetic_fastq_prename=synthetic_fastq_prename, sort_FASTQ = True)
```

### Convert FASTQ files to BAM file (optional)
Use function `scATAC_GenerateBAM.AlignSyntheticBam_Pair` to align FASTQ files onto reference genome. It takes the following arguments:
- `bowtie2_directory`: Path to software bowtie2.
- `samtools_directory`: Path to software samtools.
- `outdirectory`: Specify the output directory of the synthteic BAM file.
- `referenceGenome_name`: Base name of the reference genome FASTA file. For example, you should input "chr1" for file "chr1.fa".
- `referenceGenome_dir`: Path to the reference genome FASTA file.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scATAC_BED2FASTQ`.
- `output_BAM_pre`: Specify the base name of the output BAM file.

> **Important** Note that before using function `scATAC_GenerateBAM.AlignSyntheticBam_Pair`, the reference gemome FASTA file should be indexed by bowtie2 through following chunk and make sure the output index files are within the same directory to *`referenceGenome_name`.fa*.

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
Use function `scATAC_GenerateBAM.scATAC_ErrorBase` to introduce random error to synthetic reads. It takes the following arguments:
- `fgbio_jarfile`: Path to software fgbio jar script.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Specify the output directory of the synthteic FASTQ file with random errors.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scATAC_BED2FASTQ`.

This function will output synthetic reads with random errors in FASTQ files named as *`synthetic_fastq_prename`.ErrorIncluded.read1.bed2fa.fq*, *`synthetic_fastq_prename`.ErrorIncluded.read2.bed2fa.fq* to directory `outdirectory`.

> **Important** Note that before using function `scATAC_GenerateBAM.scATAC_ErrorBase`, please create the reference dictionary with function `CreateSequenceDictionary` using software Picard and make sure the output *.dict* files are within the same directory to *`referenceGenome_name`.fa*.

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
