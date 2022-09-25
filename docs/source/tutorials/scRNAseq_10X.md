# scReadSim on 10x scRNA-seq 

Import modules.

```{code-block} python3
import sys, os
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scRNA_GenerateBAM as scRNA_GenerateBAM
import pkg_resources
```


## Step 1: Download test sample
The example deploys scReadSim on the [10x single cell RNA-seq](https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0) dataset. The demo BAM file and its corresponding cell barcode file could be accessed through the following chunk. This BAM file uses mm10 as reference genome, the required chromosome size file is also pacakged for this example.  

```{code-block} python3
INPUT_cells_barcode_file = pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv') 
filename = "10X_RNA_chr1_3073253_4526737"
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
Note: Input BAM file for scReadSim needs pre-processing to add the cell barcode in front of the read name. For example, in 10x sequencing data, cell barcode `AACTTAGTCACAAGCT-1` is stored in the field `CB:Z:AACTTAGTCACAAGCT-1`. 

```{code-block} console
$ samtools view 10X_RNA_chr1_3073253_4526737_unprocess.bam | head -n 1
A00984:207:HGWCKDSXY:2:2306:14253:36886      16      chr1    3013015 255     17M186701N73M   *       0       0       TTTTTTTTTTTTTTGTTTTAAAATGACCACAGTGTACTTTATTTAATGATTTTTGTACTTTGTGTTGCAATAAAATAAAAAAAAAATCTA   ::F::FFFF,F:,FFF:FFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFF:::FFF:FFFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFF      NH:i:1  HI:i:1  AS:i:78      nM:i:0  RG:Z:e18_mouse_brain_fresh_5k:0:1:HGWCKDSXY:2   RE:A:I  xf:i:0  CR:Z:AACTTAGTCACAAGCT   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:AACTTAGTCACAAGCT-1 UR:Z:TAAGGTTCGACA   UY:Z:FFFFFFFFFFFF        UB:Z:TAAGGTTCGACA
```

The following code chunk adds the cell barcodes in front of the read names.

```{code-block} console
$ # extract the header file
$ mkdir tmp
$ samtools view 10X_RNA_chr1_3073253_4526737_unprocess.bam -H > tmp/10X_RNA_chr1_3073253_4526737.header.sam

$ # create a bam file with the barcode embedded into the read name
$ time(cat <( cat tmp/10X_RNA_chr1_3073253_4526737.header.sam ) \
 <( samtools view 10X_RNA_chr1_3073253_4526737_unprocess.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
 | samtools view -bS - > 10X_RNA_chr1_3073253_4526737.bam) 
$ rm -d tmp

$ samtools view 10X_RNA_chr1_3073253_4526737.bam | head -n 1
AACTTAGTCACAAGCT-1:A00984:207:HGWCKDSXY:2:2306:14253:36886      16      chr1    3013015 255     17M186701N73M   *       0       0       TTTTTTTTTTTTTTGTTTTAAAATGACCACAGTGTACTTTATTTAATGATTTTTGTACTTTGTGTTGCAATAAAATAAAAAAAAAATCTA   ::F::FFFF,F:,FFF:FFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFF:::FFF:FFFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFF      NH:i:1  HI:i:1  AS:i:78      nM:i:0  RG:Z:e18_mouse_brain_fresh_5k:0:1:HGWCKDSXY:2   RE:A:I  xf:i:0  CR:Z:AACTTAGTCACAAGCT   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:AACTTAGTCACAAGCT-1 UR:Z:TAAGGTTCGACA   UY:Z:FFFFFFFFFFFF        UB:Z:TAAGGTTCGACA
```

## Step 2: Feature space construction
For scRNA-seq, scReadSim autamatically uses gene transcript regions as the feature space. Specifically, scReasSim takes gene regions as foreground features and the copmlementary regions along the reference genome as the background features. 

### Specify input parameters 
Specify the absolute path of output directory. Create output directory if it does not exist.

```{code-block} python3
outdirectory = "/home/users/example/outputs" # use absolute path
os.mkdir(outdirectory)
```
### Generate feature sets 
Use function `scRNA_CreateFeatureSets` to generate features. This function needs user to specify

- `INPUT_bamfile`: Input BAM file for anlaysis.
- `samtools_directory`: Path to software *samtools*.
- `bedtools_directory`: Path to software *bedtools*.
- `outdirectory`: Specify the output directory of the features files.
- `genome_annotation`: Genome annotation file for the reference genome that the input BAM aligned on or the synthetic BAM should align on.
- `genome_size_file`: Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicating the size.
- `ref_peakfile`: Specify the name of output foreground feature bed file.
- `ref_comple_peakfile`: Specify the name of output background feature bed file.    

```{code-block} python3
INPUT_genome_annotation = "example/refgenome_dir/gencode.vM10.annotation.gtf" # Change the path
ref_peakfile = "%s_peaks.bed" % filename
ref_comple_peakfile = "%s_peaks.COMPLE.bed" % filename

# Generate features
Utility.scRNA_CreateFeatureSets(INPUT_bamfile=INPUT_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_annotation=INPUT_genome_annotation, genome_size_file=INPUT_genome_size_file, ref_peakfile=ref_peakfile, ref_comple_peakfile=ref_comple_peakfile)
```


## Step 3: Count matrix construction
Based on the feature sets output in **Step 2**, scReasSim constructs the count matrices for both foreground feautures and background features through function `Utility.scRNA_bam2countmat_paral`. This function needs user to specify

- `cells_barcode_file`: Cell barcode file corresponding to the input BAM file.
- `bed_file`: Features bed file to generate the count matrix.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `outdirectory`: Specify the output directory of the count matrix file.
- `count_mat_filename`: Specify the base name of output count matrix.
- `UMI_modeling`: (Optional, default: False) Specify whether scReadSim should also model UMI count of the input BAM file.
- `UMI_count_mat_filename`: (Optional, default: 'UMI_countmat') If UMI_modeling is set to True, specify the base name of output UMI count matrix.
- `n_cores`: (Optional, default: '1') Specify the number of cores for parallel computing when generating count matrix.

For the user specified `count_mat_filename`, scReadSim will generate a count matrix named *`count_mat_filename`.txt* to directory `outdirectory`. In case modeling both UMI and read count of the input BAM file (set `UMI_modeling` to be True), scReadSim generate two count matrices named *`count_mat_filename`.txt* and *`UMI_count_mat_filename`.txt*to directory `outdirectory`.

```{code-block} python3
count_mat_filename = "%s.countmatrix" % filename
count_mat_comple_filename = "%s.COMPLE.countmatrix" % filename
UMI_count_mat_filename = "%s.UMI.countmatrix" % filename

# Construct count matrix for foregroud features
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=outdirectory + "/" + ref_peakfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_filename, UMI_modeling=True, UMI_count_mat_filename = UMI_count_mat_filename, n_cores=1)
# Construct count matrix for background features
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=outdirectory + "/" + ref_comple_peakfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=count_mat_comple_filename, n_cores=1))
```

## Step 4: Synthetic count matrix simulation
The current version of scReadSim implement scDesign2 (reference) to generate synthetic count matrix based on the constructed count matrix from the input BAM file. Use function `GenerateSyntheticCount.scRNA_GenerateSyntheticCount` to generate synthetic count matrix with following paramters

- `count_mat_filename`: Base name of the count matrix output by function bam2countmat().
- `directory`: Path to the count matrix.
- `outdirectory`: Output directory of coordinate files.
- `UMI_modeling`: (Optional, default: False) Specify whether scReadSim should also model UMI count of the input BAM file.
- `UMI_count_mat_filename`: (Optional, default: 'UMI_countmat') Base name of the UMI count matrix output by function `scRNA_bam2countmat()` with option UMI_modeling setting to Ture.
- `n_cell_new`: (Optional, default: None) Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
- `total_count_new`: (Optional, default: None) Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
- `celllabel_file`: (Optional, default: None) Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.


Given the input count matrix *`count_mat_filename`.txt*, scReadSim generates two files to `outdirectory` for following analysis:

- **`count_mat_filename`.scDesign2Simulated.txt**: Synthetic count matrix.
- **`count_mat_filename`.scDesign2Simulated.nReadRegionmargional.txt**: The per-feature summation of counts for synthetic count matrix.

```{code-block} python3
# Generate synthetic count matrix for foregroud features
GenerateSyntheticCount.scRNA_GenerateSyntheticCount(count_mat_filename=count_mat_filename, directory=outdirectory, outdirectory=outdirectory, UMI_modeling=True, UMI_count_mat_filename = UMI_count_mat_filename)
# Generate synthetic count matrix for backgroud features
GenerateSyntheticCount.scRNA_GenerateSyntheticCount(count_mat_filename=count_mat_comple_filename, directory=outdirectory, outdirectory=outdirectory)
```

## Step 5: Synthetic BAM file generation

### Generate synthetic reads in BED format
Based on the synthetic count matrix, scReadSim generates synthetic reads by randomly sampling from the real BAM file input by users. First use function `scRNA_GenerateBAMCoord_paral` to create the synthetic reads and output in BED file storing the coordinates information. Function `scRNA_GenerateBAMCoord_paral` takes following input arguments:
- `count_mat_filename`: The base name of output count matrix in bam2countmat.
- `samtools_directory`: Path to software samtools.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `ref_peakfile`: Features bed file.
- `directory_cellnumber`: Directory of the marginal synthetic count vector file output in scATAC_GenerateSyntheticCount step.
- `outdirectory`: Specify the output directory for synthetic reads bed file.
- `BED_filename`: Specify the base name of output bed file.
- `OUTPUT_cells_barcode_file`: Specify the file name storing the synthetic cell barcodes.
- `read_len`: (Optional, default: '90') Specify the length of synthetic reads. Default value is 90 bp.
- `jitter_size`: (Optional, default: '5') Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
- `CB_len`: (Optional, default: '16') Specify the length of cell barcode. Default value is 16 bp.
- `UMI_modeling`: (Optional, default: False) Specify whether scReadSim should also model UMI count of the input BAM file.
- `UMI_count_mat_filename`: (Optional, default: 'UMI_countmat') Base name of the UMI count matrix output by function `scRNA_bam2countmat()` with option UMI_modeling setting to Ture.
- `UB_len`: (Optional, default: '10') Specify the length of UMI barcode. Default value is 10 bp.
- `n_cores`: (Optional, default: '1') Specify the number of cores for parallel computing.

This function will output a bed file *`BED_filename`.bed* storing the coordinates information of synthetic reads and its cell barcode file `OUTPUT_cells_barcode_file` in directory `outdirectory`.

After generation of synthetic reads for both foreground and background features, combine the two bed files using function `scRNA_GenerateBAM.scRNA_CombineBED`, which takes following input arguments:
- `outdirectory`: Directory of `BED_filename_pre`.txt and `BED_COMPLE_filename_pre`.txt.
- `BED_filename_pre`: File prename of foreground synthetic reads bed file.
- `BED_COMPLE_filename_pre`: File prename of background synthetic reads bed file.
- `BED_filename_combined_pre`: Specify the combined syntehtic reads bed file prename. The combined bed file will be output to `outdirectory`.

```{code-block} python3
directory_cellnumber = outdirectory
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
BED_filename_pre = "%s.syntheticBAM.CBincluded" % filename
BED_COMPLE_filename_pre = "%s.syntheticBAM.COMPLE.CBincluded" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined.CBincluded" % filename

# Create synthetic read coordinates for foregroud features
scRNA_GenerateBAM.scRNA_GenerateBAMCoord_paral(
	count_mat_filename=count_mat_filename, samtools_directory=samtools_directory, INPUT_bamfile=INPUT_bamfile, ref_peakfile=outdirectory + "/" + ref_peakfile, directory_cellnumber=directory_cellnumber, outdirectory=outdirectory, BED_filename=BED_filename_pre, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, read_len=90, jitter_size=5, CB_len=16, UMI_modeling=True, UMI_count_mat_file=UMI_count_mat_filename, UB_len=10, n_cores=1)
# Create synthetic read coordinates for backgroud features
scRNA_GenerateBAM.scRNA_GenerateBAMCoord_paral(
	count_mat_filename=count_mat_comple_filename, samtools_directory=samtools_directory, INPUT_bamfile=INPUT_bamfile, ref_peakfile=outdirectory + "/" + ref_comple_peakfile, directory_cellnumber=directory_cellnumber, outdirectory=outdirectory, BED_filename=BED_COMPLE_filename_pre, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, read_len=90, jitter_size=5, CB_len=16, UMI_modeling=True, UMI_count_mat_file=UMI_count_mat_filename, UB_len=10, n_cores=1)

# Combine foreground and background bed file
scRNA_GenerateBAM.scRNA_CombineBED(outdirectory=outdirectory, BED_filename_pre=BED_filename_pre, BED_COMPLE_filename_pre=BED_COMPLE_filename_pre, BED_filename_combined_pre=BED_filename_combined_pre)
```

### Convert BED files to FASTQ files
Use function `scRNA_BED2FASTQ` to convert BED file to FASTQ file. This function takes the following arguments:
- `bedtools_directory`: Path to software bedtools.
- `seqtk_directory`: Path to software seqtk.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Output directory of the synthteic bed file and its corresponding cell barcodes file.
- `BED_filename_combined`: Base name of the combined bed file output by function `scRNA_CombineBED`.
- `synthetic_fastq_prename`: Specify the base name of the output FASTQ files.
- `sort_FASTQ`: (Optional, default: True) Set `True` to sort the output FASTQ file.

This function will output paired-end reads in FASTQ files named as *`BED_filename_combined`.read1.bed2fa.fq*, *`BED_filename_combined`.read2.bed2fa.fq* to directory `outdirectory`.

```{code-block} python3
referenceGenome_name = "chr1"
referenceGenome_dir = "/home/users/example/refgenome_dir" 
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)
synthetic_fastq_prename = BED_filename_combined_pre

# Convert combined bed file into FASTQ files
scRNA_GenerateBAM.scRNA_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=BED_filename_combined_pre, synthetic_fastq_prename=synthetic_fastq_prename, sort_FASTQ = True)
```

### Convert FASTQ files to BAM file (optional)
Use function `AlignSyntheticBam_Pair` to align FASTQ files onto reference genome. It takes the following arguments:
- `bowtie2_directory`: Path to software bowtie2.
- `samtools_directory`: Path to software samtools.
- `outdirectory`: Specify the output directory of the synthteic BAM file.
- `referenceGenome_name`: Base name of the reference genome FASTA file. For example, you should input "chr1" for file "chr1.fa".
- `referenceGenome_dir`: Path to the reference genome FASTA file.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scRNA_BED2FASTQ`.
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
scRNA_GenerateBAM.AlignSyntheticBam_Single(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename, output_BAM_pre=output_BAM_pre)
```

### Introduce Error to synthetic data
Use function `scRNA_ErrorBase` to introduce random error to synthetic reads. It takes the following arguments:
- `fgbio_jarfile`: Path to software fgbio jar script.
- `INPUT_bamfile`: Input BAM file for anlaysis.
- `referenceGenome_file`: Reference genome FASTA file that the synthteic reads should align.
- `outdirectory`: Specify the output directory of the synthteic FASTQ file with random errors.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scRNA_BED2FASTQ`.

This function will output synthetic reads with random errors in FASTQ files named as *`synthetic_fastq_prename`.ErrorIncluded.read1.bed2fa.fq*, *`synthetic_fastq_prename`.ErrorIncluded.read2.bed2fa.fq* to directory `outdirectory`.

> **Important** Note that before using function `scRNA_ErrorBase`, please create the reference dictionary with function `CreateSequenceDictionary` using software Picard and make sure the output *.dict* files are within the same directory to *`referenceGenome_name`.fa*.

```{code-block} console
$ cd example/refgenome_dir # change to directory where your reference genome file is
$ java -jar /home/users/picard/build/libs/picard.jar CreateSequenceDictionary \
$       -R chr1.fa \
$       -O chr1.fa.dict
```

In the demo data, We have built the dictionary file chr1.fa.dict for chr1.fa stored in reference.genome.chr1.tar.gz. 

```{code-block} python3
# Generate reads with errors in FASTQs
scRNA_GenerateBAM.scRNA_ErrorBase(fgbio_jarfile=fgbio_jarfile, INPUT_bamfile=INPUT_bamfile, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, synthetic_fastq_prename=synthetic_fastq_prename)
# Reads alignment (optional)
scRNA_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename + ".ErrorIncluded" , output_BAM_pre=output_BAM_pre+ ".ErrorIncluded")
```








