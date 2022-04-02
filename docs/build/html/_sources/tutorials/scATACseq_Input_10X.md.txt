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

# Reference genome
referenceGenome_name = "chr1"
referenceGenome_dir = "~/Projects/scATAC_Simulator/package_development/package_data" # Put these into server FileGator for tutorial downloading 
# referenceGenome_name = "GRCm38.primary_assembly.genome"
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)

# Genome size file
# genome_size_file = "/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/mm10.chrom.sizes"
genome_size_file = pkg_resources.resource_filename("scReadSim", 'data/mm10.chrom.sizes')

# Specify parameters
filename = "10X_ATAC_chr1_4194444_4399104"
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220125_%s_NONINPUT_withCluster" % filename
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220116_%s_NONINPUT_withCluster" % filename
outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220330_10X_scATACseq_INPUT"
INPUT_bamfile = pkg_resources.resource_filename("scReadSim", 'data/%s.bam' % filename) # this one work, still need to test if R script works
input_peakfile = pkg_resources.resource_filename("scReadSim", 'data/%s.INPUT.peaks.bed' % filename) # this one work, still need to test if R script works
INPUT_cells_barcode_file = pkg_resources.resource_filename("scReadSim", 'data/barcodes.tsv') # this one work, still need to test if R script works

#####################################################################
############################ Main function #########################
#####################################################################
os.mkdir(outdirectory)
command = "cd %s" % outdirectory
os.system(command)

# directory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_data"
# filename = "10X_ATAC_chr1_4194444_4399104"
# INPUT_bamfile = "%s/%s.bam" % (directory, filename)
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
# INPUT_cells_barcode_file = "/home/gayan/Projects/scATAC_Simulator/package_development/package_data/barcodes.tsv"

######################## Prepare Input Feature Set ######################## 
input_comple_peakfile = "%s.INPUT.COMPLE.peaks.bed" % filename 
# input_peakfile = "%s.INPUT.peaks.bed" % filename
## Generate INPUT peak ## Fill in 
# Take complementary regions and stored as true-non-peak.bed
Utility.ComplementFeature(input_peakfile, input_comple_peakfile, genome_size_file, outdirectory, bedtools_directory)

######################## Reference feature set ######################## 
######################## MACS Call Peaks ######################## 
MACS3_peakname_pre = filename + ".MACS3"
ref_peakfile = "%s_peaks.bed" % MACS3_peakname_pre
ref_comple_peakfile = "%s_peaks.COMPLE.bed" % MACS3_peakname_pre
Utility.CallPeak(macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre)
######################## Generate Reference Feature Set ######################## 
Utility.scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_size_file, ref_peakfile, ref_comple_peakfile, MACS3_peakname_pre)
```


## STEP2: Count matrix construction
```{code-block} python3

######################## Generate Count matrix ######################## 
assignment_file = filename + ".assigned.peaks.txt"
assignment_comple_file=filename + ".COMPLE.assigned.peaks.txt"
Utility.match_peak(input_peakfile, outdirectory + "/" + ref_peakfile, outdirectory, assignment_file)
Utility.match_peak(outdirectory + "/" + input_comple_peakfile, outdirectory + "/" + ref_comple_peakfile, outdirectory, assignment_comple_file)
count_mat_filename = "%s.assigned.countmatrix" % filename
count_mat_comple_filename = "%s.assigned.COMPLE.countmatrix" % filename
count_mat_format = "txt"
Utility.bam2countmat_INPUT(INPUT_cells_barcode_file, outdirectory + "/" + assignment_file, INPUT_bamfile, outdirectory, count_mat_filename + "." + count_mat_format)
Utility.bam2countmat_INPUT(INPUT_cells_barcode_file, outdirectory + "/" + assignment_comple_file, INPUT_bamfile, outdirectory, count_mat_comple_filename + "." + count_mat_format)

```

## STEP3: Synthetic count matrix simulation
```{code-block} python3

######################## Synthetic Matrix Training ######################## 
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename, count_mat_format, outdirectory, outdirectory)
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_comple_filename, count_mat_format, outdirectory, outdirectory)

cellnumberfile = "%s/%s.scDesign2Simulated.nPairsRegionmargional.txt" % (outdirectory, count_mat_filename)
cellnumberfile_comple = "%s/%s.scDesign2Simulated.nPairsRegionmargional.txt" % (outdirectory, count_mat_comple_filename)
synthetic_countmat_file = "%s.scDesign2Simulated.%s" % (count_mat_filename, count_mat_format)
synthetic_countmat_file_comple = "%s.scDesign2Simulated.%s" % (count_mat_comple_filename, count_mat_format)

```

## STEP4: Synthetic BAM file generation
```{code-block} python3

######################## Generate Synthetic FASTQ file ######################## 
coordinate_file = "BAMfile_halfsampled_coordinates.txt"
coordinate_COMPLE_file = "BAMfile_halfsampled_COMPLE_coordinates.txt"
BED_filename_pre = "%s.syntheticBAM" % filename
BED_COMPLE_filename_pre = "%s.syntheticBAM.COMPLE" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined" % filename
## Parsing bam files according to referenced features, modify the position according to true features
scATAC_GenerateBAM.scATAC_SampleSyntheticReads_INPUT(samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, outdirectory + "/" + assignment_file, cellnumberfile)
scATAC_GenerateBAM.scATAC_SampleSyntheticReads_INPUT(samtools_directory, INPUT_bamfile, outdirectory, coordinate_COMPLE_file, outdirectory + "/" + assignment_comple_file, cellnumberfile_comple)
## Create synthetic read coordinates
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_INPUT(outdirectory, outdirectory + "/" + coordinate_file, outdirectory + "/" + assignment_file, outdirectory + "/" + synthetic_countmat_file, outdirectory + "/" + cellnumberfile, BED_filename_pre, OUTPUT_cells_barcode_file)
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_INPUT(outdirectory, outdirectory + "/" + coordinate_COMPLE_file, outdirectory + "/" + assignment_comple_file, outdirectory + "/" + synthetic_countmat_file_comple, outdirectory + "/" + cellnumberfile_comple, BED_COMPLE_filename_pre, OUTPUT_cells_barcode_file)

## Combine peak and comple.peak 
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory, BED_filename_pre, BED_COMPLE_filename_pre, BED_filename_combined_pre)

## Convert bed files to FASTQ files
output_BAM_pre = "%s.syntheticBAM.CBincluded" % filename
scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory, seqtk_directory, referenceGenome_file, outdirectory, BED_filename_combined_pre, sort_FASTQ = True)

######################## Generate BAM file ######################## 
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory, samtools_directory, outdirectory, referenceGenome_name, referenceGenome_dir, BED_filename_combined_pre, output_BAM_pre, doIndex = False)

```