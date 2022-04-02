# scReadSim on 10X scRNA-seq 

Import modules.

```{code-block} python3
import sys, os
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scRNA_GenerateBAM as scRNA_GenerateBAM
import pkg_resources
```


## Step 1: Download test sample. 
The example deploys scReadSim on the 10X single cell RNA-seq (give url to original website). The demo BAM file and its corresponding cell barcode file could be accessed through the following chunk. This BAM file uses mm10 as reference genome, the required chromosome size file is also embedded within the package.  

```{code-block} python3
INPUT_cells_barcode_file = pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv') 
filename = "10X_RNA_chr1_3073253_4526737"
INPUT_bamfile = pkg_resources.resource_filename("scReadSim", 'data/%s.bam' % filename) # this one work, still need to test if R script works
INPUT_genome_size_file = pkg_resources.resource_filename("scReadSim", 'data/mm10.chrom.sizes')
```

Other required files for this example inlcuding the reference genome FASTA file and annotation gtf file are downloadable through the following chunk.  

```{code-block} bash
wget http://compbio10data.stat.ucla.edu/repository/gayan/Projects/scReadSim/reference.genmoe.tar.gz # 292 MB
wget http://compbio10data.stat.ucla.edu/repository/gayan/Projects/scReadSim/gencode.vM10.annotation.gtf # 765MB
```

## Step 2: Feature space construction
For scRNA-seq, scReadSim autamatically uses gene transcript regions as the feature space. Specifically, scReasSim takes gene regions as foreground features and the copmlementary regions along the reference genome as the background features. 

Use function `scRNA_CreateFeatureSets` to generate features. This function needs user to specify

- `INPUT_bamfile`: Input BAM file for anlaysis.
- `samtools_directory`: Path to software *samtools*.
- `bedtools_directory`: Path to software *bedtools*.
- `outdirectory`: Specify the output directory of the features files.
- `genome_annotation`: Genome annotation file for the reference genome that the input BAM aligned on or the synthetic BAM should align on.
- `genome_size_file`: Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicating the size.
- `ref_peakfile`: Specify the name of output foreground feature bed file.
- `ref_comple_peakfile`: Specify the name of output background feature bed file.    


### Specify input parameters 

```{code-block} python3
INPUT_genome_annotation = "gencode.vM10.annotation.gtf" # Change the path
outdirectory = "example/outputs"
ref_peakfile = "%s_peaks.bed" % filename
ref_comple_peakfile = "%s_peaks.COMPLE.bed" % filename
```
### Generate feature sets 

```{code-block} python3
####################### Generate Feature Set ######################## 
Utility.scRNA_CreateFeatureSets(INPUT_bamfile=INPUT_bamfile, samtools_directory=samtools_directory, bedtools_directory=bedtools_directory, outdirectory=outdirectory, genome_annotation=INPUT_genome_annotation, genome_size_file=INPUT_genome_size_file, ref_peakfile=ref_peakfile, ref_comple_peakfile=ref_comple_peakfile)

```

```
referenceGenome_name = "chr1"
referenceGenome_dir = "~/Projects/scATAC_Simulator/package_development/package_data" 
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)

samtools_directory="/home/gayan/Tools/samtools/bin"
macs3_directory="/home/gayan/.local/bin"
bedtools_directory="/home/gayan/Tools/bedtools/bedtools2/bin"
seqtk_directory = "~/Tools/seqtk/seqtk"
bowtie2_directory = "/usr/bin"

OUTPUT_cells_barcode_file = outdirectory + "/synthetic_cell_barcode.txt"

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
The current version of scReadSim implement scDesign2 (reference) to generate synthetic count matrix based on the constructed count matrix from the input BAM file. Use function `GenerateSyntheticCount.scRNA_SampleSyntheticReads` to generate synthetic count matrix with following paramters

- `count_mat_filename`: Base name of the count matrix output by function bam2countmat().
- `directory`: Path of the count matrix.
- `outdirectory`: Output directory of coordinate files.
- `cluster_prestep`: Set `cluster_prestep=True` to perform a Louvain clustering before implementing scDesign2.

Given the input count matrix `count_mat_filename`.txt, scReadSim generates two files to `outdirectory` for following analysis:

- **`count_mat_filename`.scDesign2Simulated.txt**: Synthetic count matrix.
- **`count_mat_filename`.scDesign2Simulated.nReadRegionmargional.txt**: The per-feature summation of counts for synthetic count matrix.

```{code-block} python3

GenerateSyntheticCount.scRNA_SampleSyntheticReads(count_mat_filename=count_mat_filename, directory=outdirectory, outdirectory=outdirectory, cluster_prestep=True)
GenerateSyntheticCount.scRNA_SampleSyntheticReads(count_mat_filename=count_mat_comple_filename, directory=outdirectory, outdirectory=outdirectory, cluster_prestep=True)

```

## STEP4: Synthetic BAM file generation
```{code-block} python3
cellnumberfile = "%s/%s.scDesign2Simulated.nReadRegionmargional.txt" % (outdirectory, count_mat_filename)
cellnumberfile_comple = "%s/%s.scDesign2Simulated.nReadRegionmargional.txt" % (outdirectory, count_mat_comple_filename)
synthetic_countmat_file = "%s.scDesign2Simulated.%s" % (count_mat_filename, count_mat_format)
synthetic_countmat_file_comple = "%s.scDesign2Simulated.%s" % (count_mat_comple_filename, count_mat_format)

######################## Generate Synthetic FASTQ file ######################## 
coordinate_file = "BAMfile_coordinates.txt"
coordinate_COMPLE_file = "BAMfile_COMPLE_coordinates.txt"
BED_filename_pre = "%s.syntheticBAM.CBincluded" % filename
BED_COMPLE_filename_pre = "%s.syntheticBAM.COMPLE.CBincluded" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined.CBincluded" % filename

## Parsing bam files according to referenced features, modify the position according to true features
scRNA_GenerateBAM.scRNA_GenerateSyntheticReads(coordinate_file, samtools_directory, INPUT_bamfile, outdirectory, outdirectory + "/" + ref_peakfile, cellnumberfile)
scRNA_GenerateBAM.scRNA_GenerateSyntheticReads(coordinate_COMPLE_file, samtools_directory, INPUT_bamfile, outdirectory, outdirectory + "/" + ref_comple_peakfile, cellnumberfile_comple)
## Create synthetic read coordinates
scRNA_GenerateBAM.scRNA_GenerateBAMCoord(outdirectory, coordinate_file, ref_peakfile, synthetic_countmat_file, cellnumberfile, BED_filename_pre, OUTPUT_cells_barcode_file)
scRNA_GenerateBAM.scRNA_GenerateBAMCoord(outdirectory, coordinate_COMPLE_file, ref_comple_peakfile, synthetic_countmat_file_comple, cellnumberfile_comple, BED_COMPLE_filename_pre, OUTPUT_cells_barcode_file)

## Combine peak and comple.peak 
scRNA_GenerateBAM.scRNA_CombineBED(outdirectory, BED_filename_pre, BED_COMPLE_filename_pre, BED_filename_combined_pre)


## Convert bed files to FASTQ files
output_BAM_pre = "%s.syntheticBAM.CBincluded" % filename
scRNA_GenerateBAM.scRNA_BED2FASTQ(bedtools_directory, seqtk_directory, referenceGenome_file, outdirectory, BED_filename_combined_pre, sort_FASTQ = True)

######################## Generate BAM file ######################## 
scRNA_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory, samtools_directory, outdirectory, referenceGenome_name, referenceGenome_dir, BED_filename_combined_pre, output_BAM_pre, doIndex = False)

```