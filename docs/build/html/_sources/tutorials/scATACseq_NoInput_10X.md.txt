# scReadSim on 10X scATAC-seq 

Import modules.

```{code-block} python3
import sys, os
import scReadSim.scATAC_GenerateBAM as scATAC_GenerateBAM
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import pkg_resources
```

## STEP1: Download test example.
```{code-block} bash
wget http://renlab.sdsc.edu/r3fang/share/SnapTools/snaptools_test.tar.gz
tar -xf snaptools_test.tar.gz
cd snaptools_test/
gunzip mm10.fa.gz
```

## STEP1: Feature space construction
```{code-block} python3
#####################################################################
############################## User Input #########################
#####################################################################
samtools_directory="/home/gayan/Tools/samtools/bin"
macs3_directory="/home/gayan/.local/bin"
bedtools_directory="/home/gayan/Tools/bedtools/bedtools2/bin"
seqtk_directory = "~/Tools/seqtk/seqtk"
bowtie2_directory = "/usr/bin"

filename = "10X_ATAC_chr1_4194444_4399104"
INPUT_cells_barcode_file = pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv') 
INPUT_bamfile = pkg_resources.resource_filename("scReadSim", 'data/%s.bam' % filename) # this one work, still need to test if R script works
outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220330_10X_scATACseq_NONINPUT"

referenceGenome_name = "chr1"
referenceGenome_dir = "~/Projects/scATAC_Simulator/package_development/package_data" # Put these into server FileGator for tutorial downloading 
# referenceGenome_name = "GRCm38.primary_assembly.genome"
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)


#####################################################################
############################ Main function #########################
#####################################################################
os.mkdir(outdirectory)
command = "cd %s" % outdirectory
os.system(command)

# INPUT_bamfile = "%s/%s.bam" % (directory, filename)
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"

######################## MACS Call Peaks ######################## 
MACS3_peakname_pre = filename + ".MACS3"
Utility.CallPeak(macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre)

######################## Generate Feature Set ######################## 
ref_peakfile = "%s_peaks.bed" % MACS3_peakname_pre
ref_comple_peakfile = "%s_peaks.COMPLE.bed" % MACS3_peakname_pre
# genome_size_file = "/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/mm10.chrom.sizes.removed"
genome_size_file = pkg_resources.resource_filename("scReadSim", 'data/mm10.chrom.sizes')
Utility.scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_size_file, ref_peakfile, ref_comple_peakfile, MACS3_peakname_pre)
 
```

## STEP2: Count matrix construction
```{code-block} python3
######################## Generate Count matrix ########################
count_mat_filename = "%s.countmatrix" % filename
count_mat_comple_filename = "%s.COMPLE.countmatrix" % filename
Utility.bam2countmat(INPUT_cells_barcode_file, outdirectory + "/" + ref_peakfile, INPUT_bamfile, outdirectory, count_mat_filename + ".txt")
Utility.bam2countmat(INPUT_cells_barcode_file, outdirectory + "/" + ref_comple_peakfile, INPUT_bamfile, outdirectory, count_mat_comple_filename + ".txt")
```

## STEP3: Synthetic count matrix simulation
```{code-block} python3
######################## Synthetic Matrix Training ######################## 
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename, count_mat_format, outdirectory, outdirectory)
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_comple_filename, count_mat_format, outdirectory, outdirectory)

cellnumberfile = "%s/%s.scDesign2Simulated.nPairsRegionmargional.txt" % (outdirectory, count_mat_filename)
cellnumberfile_comple = "%s/%s.scDesign2Simulated.nPairsRegionmargional.txt" % (outdirectory, count_mat_comple_filename)
synthetic_countmat_file = "%s.scDesign2Simulated.%s" % (count_mat_filename)
synthetic_countmat_file_comple = "%s.scDesign2Simulated.txt" % (count_mat_comple_filename)
```

## STEP4: Synthetic BAM file generation
```{code-block} python3
######################## Generate Synthetic FASTQ file ######################## 
coordinate_file = "BAMfile_halfsampled_coordinates.txt"
coordinate_COMPLE_file = "BAMfile_halfsampled_COMPLE_coordinates.txt"
BED_filename_pre = "%s.syntheticBAM.CBincluded" % filename
BED_COMPLE_filename_pre = "%s.syntheticBAM.COMPLE.CBincluded" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined.CBincluded" % filename
## Parsing bam files according to referenced features, modify the position according to true features
scATAC_GenerateBAM.scATAC_SampleSyntheticReads(coordinate_file, samtools_directory, INPUT_bamfile, outdirectory, outdirectory + "/" + ref_peakfile, cellnumberfile)
scATAC_GenerateBAM.scATAC_SampleSyntheticReads(coordinate_COMPLE_file, samtools_directory, INPUT_bamfile, outdirectory, outdirectory + "/" + ref_comple_peakfile, cellnumberfile_comple)
## Create synthetic read coordinates
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(outdirectory, outdirectory + "/" + coordinate_file, outdirectory + "/" + ref_peakfile, outdirectory + "/" + synthetic_countmat_file, outdirectory + "/" + cellnumberfile, BED_filename_pre, OUTPUT_cells_barcode_file)
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(outdirectory, outdirectory + "/" + coordinate_COMPLE_file, outdirectory + "/" + ref_comple_peakfile, outdirectory + "/" + synthetic_countmat_file_comple, outdirectory + "/" + cellnumberfile_comple, BED_COMPLE_filename_pre, OUTPUT_cells_barcode_file)

## Combine peak and comple.peak 
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory, BED_filename_pre, BED_COMPLE_filename_pre, BED_filename_combined_pre)

## Convert bed files to FASTQ files
output_BAM_pre = "%s.syntheticBAM.CBincluded" % filename
scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory, seqtk_directory, referenceGenome_file, outdirectory, BED_filename_combined_pre, sort_FASTQ = True)

######################## Generate BAM file ######################## 
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory, samtools_directory, outdirectory, referenceGenome_name, referenceGenome_dir, BED_filename_combined_pre, output_BAM_pre)

```