# scReadSim on 10x scRNA-seq 

Import modules.

```{code-block} python3
import sys, os
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scRNA_GenerateBAM as scRNA_GenerateBAM
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
To pre-process real scRNA-seq data for training, scReadSim requires a BAM file (containing scRNA-seq reads in cells) and a gene annotation file (in GTF format). Based on the gene coordinates in the annotation file, scReadSim segregates the reference genome into two sets of features: genes and inter-genes.

### Specify output directory
Specify the absolute path of output directory. Create output directory if it does not exist.

```{code-block} python3
outdirectory = "/home/users/example/outputs" # use absolute path
os.mkdir(outdirectory)
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

```{code-block} python3
INPUT_genome_annotation = "example/refgenome_dir/gencode.vM10.annotation.gtf" # Change the path
gene_bedfile = "scReadSim.Gene.bed"
intergene_bedfile = "scReadSim.InterGene.bed"

# Generate features
Utility.scRNA_CreateFeatureSets(INPUT_bamfile=INPUT_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_annotation=INPUT_genome_annotation, genome_size_file=INPUT_genome_size_file)
```


## Step 3: UMI Count matrix construction
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
UMI_gene_count_mat_filename = "%s.gene.countmatrix" % filename
UMI_intergene_count_mat_filename = "%s.intergene.countmatrix" % filename

# Construct count matrix for genes
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=outdirectory + "/" + gene_bedfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=UMI_gene_count_mat_filename, UMI_modeling=True, UMI_tag = "UB:Z", n_cores=8)
# Construct count matrix for inter-genes
Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file, bed_file=outdirectory + "/" + intergene_bedfile, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, count_mat_filename=UMI_intergene_count_mat_filename, UMI_modeling=True, UMI_tag = "UB:Z", n_cores=8)
```


## Step 4: Synthetic count matrix simulation
The current version of scReadSim implement scDesign2 (reference) to generate synthetic count matrix based on the constructed count matrix from the input BAM file. Use function `GenerateSyntheticCount.scRNA_GenerateSyntheticCount` to generate synthetic count matrix with following paramters

- `count_mat_filename`: Base name of the count matrix output by function `Utility.scRNA_bam2countmat_paral`.
- `directory`: Path to the count matrix.
- `outdirectory`: Output directory of coordinate files.
- `n_cell_new`: (Optional, default: None) Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
- `total_count_new`: (Optional, default: None) Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
- `celllabel_file`: (Optional, default: None) Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.

Given the input count matrix *`count_mat_filename`.txt*, scReadSim generates the syntheitic count matrix file to `outdirectory` for following analysis:

- **`count_mat_filename`.scDesign2Simulated.txt**: Synthetic count matrix.


```{code-block} python3
# Generate synthetic count matrix for gene-by-cell count matrix
GenerateSyntheticCount.scRNA_GenerateSyntheticCount(count_mat_filename=UMI_gene_count_mat_filename, directory=outdirectory, outdirectory=outdirectory)
# Generate synthetic count matrix for non-gene-by-cell count matrix
GenerateSyntheticCount.scRNA_GenerateSyntheticCount(count_mat_filename=UMI_intergene_count_mat_filename, directory=outdirectory, outdirectory=outdirectory)
```

## Step 5: Synthetic BAM file generation

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
        bed_file=outdirectory + "/" + gene_bedfile, UMI_count_mat_file=outdirectory + "/" + synthetic_countmat_gene_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=gene_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=90)
# Create synthetic read coordinates for intergenes
scRNA_GenerateBAM.scRNA_GenerateBAMCoord(
        bed_file=outdirectory + "/" + intergene_bedfile, UMI_count_mat_file=outdirectory + "/" + synthetic_countmat_intergene_file, synthetic_cell_label_file=outdirectory + "/" + synthetic_cell_label_file, read_bedfile_prename=intergene_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=90)

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

```{code-block} python3
referenceGenome_name = "chr1"
referenceGenome_dir = "/home/users/example/refgenome_dir" 
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)
synthetic_fastq_prename = BED_filename_combined_pre

# Convert combined bed file into FASTQ files
scRNA_GenerateBAM.scRNA_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory, BED_filename_combined=BED_filename_combined_pre, synthetic_fastq_prename=synthetic_fastq_prename)
```

### Convert FASTQ files to BAM file (optional)
Use function `scRNA_GenerateBAM.AlignSyntheticBam_Single` to align FASTQ files onto reference genome. It takes the following arguments:
- `bowtie2_directory`: Path to software bowtie2.
- `samtools_directory`: Path to software samtools.
- `outdirectory`: Specify the output directory of the synthteic BAM file.
- `referenceGenome_name`: Base name of the reference genome FASTA file. For example, you should input "chr1" for file "chr1.fa".
- `referenceGenome_dir`: Path to the reference genome FASTA file.
- `synthetic_fastq_prename`: Base name of the synthetic FASTQ files output by function `scRNA_BED2FASTQ`.
- `output_BAM_pre`: Specify the base name of the output BAM file.

> **Important** Note that before using function `scRNA_GenerateBAM.AlignSyntheticBam_Single`, the reference gemome FASTA file should be indexed by bowtie2 through following chunk and make sure the output index files are within the same directory to *`referenceGenome_name`.fa*.

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

> **Important** Note that before using function `scRNA_ErrorBase`, please create the reference dictionary for the reference genome with function `CreateSequenceDictionary` using software Picard and make sure the output *.dict* files are within the same directory to *`referenceGenome_name`.fa*.

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
scRNA_GenerateBAM.AlignSyntheticBam_Single(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename + ".ErrorIncluded" , output_BAM_pre=output_BAM_pre+ ".ErrorIncluded")
```








