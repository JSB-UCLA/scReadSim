import sys, os
import subprocess
import tempfile
import scReadSim.scATAC_GenerateBAM as scATAC_GenerateBAM
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount

#####################################################################
############################## User Input #########################
#####################################################################
samtools_directory="/home/gayan/Tools/samtools/bin"
macs3_directory="/home/gayan/.local/bin"
bedtools_directory="/home/gayan/Tools/bedtools/bedtools2/bin"
seqtk_directory = "~/Tools/seqtk/seqtk"
bowtie2_directory = "/usr/bin"

directory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_data"
filename = "10X_ATAC_chr1_4194444_4399104"
INPUT_cells_barcode_file = "/home/gayan/Projects/scATAC_Simulator/package_development/package_data/barcodes.tsv"
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220125_%s_NONINPUT_withCluster" % filename
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220116_%s_NONINPUT_withCluster" % filename
outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220311_10X_scATACseq_INPUT"

#####################################################################
############################ Main function #########################
#####################################################################
os.mkdir(outdirectory)
command = "cd %s" % outdirectory
os.system(command)

INPUT_bamfile = "%s/%s.bam" % (directory, filename)
genome_size_file = "/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/mm10.chrom.sizes.removed"
OUTPUT_cells_barcode_file = outdirectory + "/synthetic_cell_barcode.txt"


######################## Prepare Input Feature Set ######################## 
input_peakfile = "%s.INPUT.peaks.bed" % filename
input_comple_peakfile = "%s.INPUT.COMPLE.peaks.bed" % filename 

## Generate INPUT peak ## Fill in 
# Take complementary regions and stored as true-non-peak.bed
Utility_INPUT.ComplementFeature(outdirectory + "/" + input_peakfile, input_comple_peakfile, genome_size_file, outdirectory, bedtools_directory)

######################## Reference feature set ######################## 
######################## MACS Call Peaks ######################## 
MACS3_peakname_pre = filename + ".MACS3"
ref_peakfile = "%s_peaks.bed" % MACS3_peakname_pre
ref_comple_peakfile = "%s_peaks.COMPLE.bed" % MACS3_peakname_pre
Utility.CallPeak(macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre)
######################## Generate Reference Feature Set ######################## 
Utility.scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_size_file, ref_peakfile, ref_comple_peakfile, MACS3_peakname_pre)

######################## Generate Count matrix ######################## 
assignment_file = filename + ".assigned.peaks.txt"
assignment_comple_file=filename + ".COMPLE.assigned.peaks.txt"
Utility_INPUT.match_peak(outdirectory + "/" + input_peakfile, outdirectory + "/" + ref_peakfile, outdirectory, assignment_file)
Utility_INPUT.match_peak(outdirectory + "/" + input_comple_peakfile, outdirectory + "/" + ref_comple_peakfile, outdirectory, assignment_comple_file)
count_mat_filename = "%s.assigned.countmatrix" % filename
count_mat_comple_filename = "%s.assigned.COMPLE.countmatrix" % filename
count_mat_format = "txt"
Utility_INPUT.bam2countmat(INPUT_cells_barcode_file, outdirectory, assignment_file, INPUT_bamfile, outdirectory, count_mat_filename + "." + count_mat_format)
Utility_INPUT.bam2countmat(INPUT_cells_barcode_file, outdirectory, assignment_comple_file, INPUT_bamfile, outdirectory, count_mat_comple_filename + "." + count_mat_format)

######################## Synthetic Matrix Training ######################## 
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename, count_mat_format, outdirectory, outdirectory)
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_comple_filename, count_mat_format, outdirectory, outdirectory)

cellnumberfile = "%s/%s.scDesign2Simulated.nPairsRegionmargional.txt" % (outdirectory, count_mat_filename)
cellnumberfile_comple = "%s/%s.scDesign2Simulated.nPairsRegionmargional.txt" % (outdirectory, count_mat_comple_filename)
synthetic_countmat_file = "%s.scDesign2Simulated.%s" % (count_mat_filename, count_mat_format)
synthetic_countmat_file_comple = "%s.scDesign2Simulated.%s" % (count_mat_comple_filename, count_mat_format)

######################## Generate Synthetic FASTQ file ######################## 
coordinate_file = "BAMfile_halfsampled_coordinates.txt"
coordinate_COMPLE_file = "BAMfile_halfsampled_COMPLE_coordinates.txt"
BED_filename_pre = "%s.syntheticBAM" % filename
BED_COMPLE_filename_pre = "%s.syntheticBAM.COMPLE" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined" % filename
## Parsing bam files according to referenced features, modify the position according to true features
scATAC_GenerateBAM.INPUT_GenerateSyntheticReads(samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, assignment_file, cellnumberfile)
scATAC_GenerateBAM.INPUT_GenerateSyntheticReads(samtools_directory, INPUT_bamfile, outdirectory, coordinate_COMPLE_file, assignment_comple_file, cellnumberfile_comple)
## Create synthetic read coordinates
scATAC_GenerateBAM.scATAC_INPUT_GenerateBAMCoord(outdirectory, coordinate_file, assignment_file, synthetic_countmat_file, cellnumberfile, BED_filename_pre, OUTPUT_cells_barcode_file)
scATAC_GenerateBAM.scATAC_INPUT_GenerateBAMCoord(outdirectory, coordinate_COMPLE_file, assignment_comple_file, synthetic_countmat_file_comple, cellnumberfile_comple, BED_COMPLE_filename_pre, OUTPUT_cells_barcode_file)

## Combine peak and comple.peak 
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory, BED_filename_pre, BED_COMPLE_filename_pre, BED_filename_combined_pre)

## Convert bed files to FASTQ files
referenceGenome_name = "GRCm38.primary_assembly.genome"
referenceGenome_dir = "/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE"
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)
output_BAM_pre = "%s.syntheticBAM.CBincluded" % filename
scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory, seqtk_directory, referenceGenome_file, outdirectory, BED_filename_combined_pre, sort_FASTQ = True)

######################## Generate BAM file ######################## 
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory, samtools_directory, outdirectory, referenceGenome_name, referenceGenome_dir, BED_filename_combined_pre, output_BAM_pre, doIndex = False)






































