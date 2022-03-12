import sys, os
import subprocess
import tempfile
import scReadSim.scATAC_GenerateBAM as scATAC_GenerateBAM
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount

#####################################################################
############################ Prepare the data #######################
#####################################################################
# full_filename=e18_mouse_brain_fresh_5k_atac_possorted_bam
# full_directory=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA
# directory=${full_directory}/split.${full_filename}
########################  Prestep (by user): Attach barcodes to BAM file ########################
# mkdir ${directory}
# mkdir ${full_directory}/tmp
# ${samtools_directory}/samtools view ${full_directory}/${full_filename}.bam -H > ${full_directory}/tmp/${full_filename}.header.sam # extract the header file
# time(cat <( cat ${full_directory}/tmp/${full_filename}.header.sam ) \
# <( ${samtools_directory}/samtools view ${full_directory}/${full_filename}.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
# | ${samtools_directory}/samtools view -bS - > ${full_directory}/${full_filename}.scReadSim.bam) # create a bam file with the barcode embedded into the read name 241m40.507s
# rm -d ${full_directory}/tmp/
########################  Prepare Genome sizes ########################
# wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes
# cat mm10.chrom.sizes | grep -v "_" > mm10.chrom.sizes.removed 
######################## Split, Index and Sort BAM file ######################## 
# ${samtools_directory}/samtools sort ${full_directory}/${full_filename}.scReadSim.bam -o ${full_directory}/${full_filename}.scReadSim.sorted.bam
# ${samtools_directory}/samtools index ${full_directory}/${full_filename}.scReadSim.sorted.bam ${full_directory}/${full_filename}.scReadSim.sorted.bam.bai 
# ## Separate bam files
# for (( counter=1; counter<20; counter++ ))
# do
# time(${samtools_directory}/samtools view -b ${full_directory}/${full_filename}.scReadSim.sorted.bam chr$counter > ${directory}/${full_filename}_chr${counter}.bam)
# time(${samtools_directory}/samtools index ${directory}/${full_filename}_chr${counter}.bam)
# printf "chr$counter\tDone\n"
# done
# printf "Done\n"
# time(${samtools_directory}/samtools view -b ${full_directory}/${full_filename}.scReadSim.sorted.bam chrX > "${directory}/${full_filename}_chrX.bam")
# time(${samtools_directory}/samtools index "${directory}/${full_filename}_chrX.bam")
# printf "chrX\tDone\n"
# time(${samtools_directory}/samtools view -b ${full_directory}/${full_filename}.scReadSim.sorted.bam chrY > "${directory}/${full_filename}_chrY.bam")
# time(${samtools_directory}/samtools index "${directory}/${full_filename}_chrY.bam")
# printf "chrY\tDone\n"
# time(${samtools_directory}/samtools view -b ${full_directory}/${full_filename}.scReadSim.sorted.bam chrM > "${directory}/${full_filename}_chrM.bam")
# time(${samtools_directory}/samtools index "${directory}/${full_filename}_chrM.bam")
# printf "chrM\tDone\n"



#####################################################################
############################## User Input #########################
#####################################################################
samtools_directory="/home/gayan/Tools/samtools/bin"
macs3_directory="/home/gayan/.local/bin"
bedtools_directory="/home/gayan/Tools/bedtools/bedtools2/bin"
seqtk_directory = "~/Tools/seqtk/seqtk"
bowtie2_directory = "/usr/bin"

# path = os.path.abspath()
# new_path = os.xxx($path/xxx)

directory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_data"
filename = "10X_ATAC_chr1_4194444_4399104"
INPUT_cells_barcode_file = "/home/gayan/Projects/scATAC_Simulator/package_development/package_data/barcodes.tsv"
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220125_%s_NONINPUT_withCluster" % filename
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220116_%s_NONINPUT_withCluster" % filename
outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220310_10X_scATACseq_NONINPUT"

#####################################################################
############################ Main function #########################
#####################################################################
os.mkdir(outdirectory)
command = "cd %s" % outdirectory
os.system(command)

INPUT_bamfile = "%s/%s.bam" % (directory, filename)
OUTPUT_cells_barcode_file = outdirectory + "/synthetic_cell_barcode.txt"

######################## MACS Call Peaks ######################## 
MACS3_peakname_pre = filename + ".MACS3"
Utility.CallPeak(macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre)

######################## Generate Feature Set ######################## 
ref_peakfile = "%s_peaks.bed" % MACS3_peakname_pre
ref_comple_peakfile = "%s_peaks.COMPLE.bed" % MACS3_peakname_pre
genome_size_file = "/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/mm10.chrom.sizes.removed"
Utility.scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_size_file, ref_peakfile, ref_comple_peakfile, MACS3_peakname_pre)
  
######################## Generate Count matrix ########################
count_mat_filename = "%s.countmatrix" % filename
count_mat_comple_filename = "%s.COMPLE.countmatrix" % filename
count_mat_format = "txt"
Utility.bam2countmat(INPUT_cells_barcode_file, outdirectory, ref_peakfile, INPUT_bamfile, outdirectory, "%s.%s" % (count_mat_filename, count_mat_format))
Utility.bam2countmat(INPUT_cells_barcode_file, outdirectory, ref_comple_peakfile, INPUT_bamfile, outdirectory, "%s.%s" % (count_mat_comple_filename, count_mat_format))

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
BED_filename_pre = "%s.syntheticBAM.CBincluded" % filename
BED_COMPLE_filename_pre = "%s.syntheticBAM.COMPLE.CBincluded" % filename
BED_filename_combined_pre = "%s.syntheticBAM.combined.CBincluded" % filename
## Parsing bam files according to referenced features, modify the position according to true features
scATAC_GenerateBAM.GenerateSyntheticReads(samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, ref_peakfile, cellnumberfile)
scATAC_GenerateBAM.GenerateSyntheticReads(samtools_directory, INPUT_bamfile, outdirectory, coordinate_COMPLE_file, ref_comple_peakfile, cellnumberfile_comple)
## Create synthetic read coordinates
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(outdirectory, coordinate_file, ref_peakfile, synthetic_countmat_file, cellnumberfile, BED_filename_pre, OUTPUT_cells_barcode_file)
scATAC_GenerateBAM.scATAC_GenerateBAMCoord(outdirectory, coordinate_COMPLE_file, ref_comple_peakfile, synthetic_countmat_file_comple, cellnumberfile_comple, BED_COMPLE_filename_pre, OUTPUT_cells_barcode_file)

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










