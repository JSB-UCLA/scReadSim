import pandas as pd
import numpy as np
import csv
import collections
import time
import sys
import os
pd.options.mode.chained_assignment = None  # default='warn'
import string
import random
import subprocess
from tqdm import tqdm

def flatten(x):
    """Flatten a nested list.

    """
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]


def cellbarcode_generator(length, size=10):
	"""Generate random cellbarcode.

    Parameters
    ----------
    length: `int`
        Number of cells.
    size: `int` (default: '10')
        Size of cell barcode. Default value is 10 bp.

    Return
    ------
    cb_list: `list`
        List of randomly generated cell barcodes.
    """
	chars = 'ACGT'
	cb_list = [''.join(random.choice(chars) for _ in range(size)) for cell in range(length)]
	return cb_list



def scRNA_SampleSyntheticReads(count_mat_filename, samtools_directory, INPUT_bamfile, outdirectory, ref_peakfile, directory_cellnumber):
	"""Sample Synthetic reads. The sampled reads coordinates are stored as `coordinate_file` in `outdirectory`.

	Parameters
	----------
	count_mat_filename: `str`
		The base name of output count matrix in bam2countmat.
	samtools_directory: `str`
		Path to software samtools.
	INPUT_bamfile: `str`
		Path to input BAM file.
	outdirectory: `str`
		Output directory of coordinate files.
	ref_peakfile: `str`
		Features bed file.
	directory_cellnumber: `str`
		Directory of the marginal synthetic count vector file output in scATAC_GenerateSyntheticCount step.
	"""
	coordinate_file = count_mat_filename + ".BAMfile_coordinates.txt"
	cellnumberfile = "%s/%s.scDesign2Simulated.nReadRegionmargional.txt" % (directory_cellnumber, count_mat_filename)
	rm_coor_command = "rm %s/%s" % (outdirectory, coordinate_file)
	os.system(rm_coor_command)
	create_coor_command = "touch %s/%s" % (outdirectory, coordinate_file)
	os.system(create_coor_command)
	cmd = "while true; do read -r region <&3 || break;  read -r ncell <&4 || break; region=$(echo ${region} | cut -f 1,2,3 | perl -lane 'print \"$F[0]:$F[1]-$F[2]\"');  paste -d\"\t\" <(awk -v nsample=${ncell} -v region=${region} 'BEGIN{for(c=0;c<nsample;c++) print region}') <(%s/samtools view %s ${region} | shuf -r -n ${ncell} | cut -f3,4) >> %s/%s;  done 3<%s 4<%s" % (samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, ref_peakfile, cellnumberfile)
	# Testing using copied directory
	output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to generate synthetic reads:\n', error.decode())



def scRNA_PerTruePeakEdition(peak_record, count_vec, read_lines, random_cellbarcode_list, read_len=80, jitter_size=5):
	"""Formulate Synthetic reads.

	Parameters
	----------
	peak_record: `numpy.array`
		Coordinates of a query peak.
	count_vec: `numpy.array`
    	Count vector of the query peak.
	read_lines: `pandas.dataframe`
		Coordinates of synthetic reads sampled from the input BAM file.
	random_cellbarcode_list: `list`
    	List of cell barcodes randomly generated using [cellbarcode_generator().
	read_len: `int` (default: '50')
		Specify the length of synthetic reads. Default value is 50 bp.
	jitter_size: `int` (default: '5')
		Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.

	Return
	------
	read_1_df: `pandas.dataframe`
		Coordinates of synthetic read 1 generated for the query peak.
	read_2_df: `pandas.dataframe`
		Coordinates of synthetic read 2 generated for the query peak.
	"""
	# peak_record = peaks_assignments.loc[1,] # Input
	true_peak_concat = peak_record[0] + ":" + str(peak_record[1]) + "-" + str(peak_record[2])
	reads_cur = read_lines[read_lines['peak_name'] == true_peak_concat] # Input
	nread_cur= np.sum(count_vec).astype(int)
	# Add cell information
	nonempty_cell_ind = np.where(count_vec != 0)[0]
	random_umi_list = cellbarcode_generator(nread_cur, size=10)
	read_code_simu_cur = [random_cellbarcode_list[nonempty_cell_ind[ind]] + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(true_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_vec[nonempty_cell_ind[ind]])]
	read_code_simu_cur_withUMI = [read_code_simu_cur[read_id].split(":",1)[0] + random_umi_list[read_id] + ":" + read_code_simu_cur[read_id].split(":",1)[1] for read_id in range(len(read_code_simu_cur))]
	# read_code_simu_cur = ["CellType1" + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(true_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_vec[nonempty_cell_ind[ind]])]
	# start = time.time()
	jitter_value_vec = np.random.random_integers(-jitter_size,jitter_size,size=np.shape(reads_cur)[0])  # nrow(reads_cur) should equal to nfrag_cur
	reads_cur['r1_start_shifted'] = reads_cur['r1_start'] + jitter_value_vec
	reads_cur['r1_end_shifted'] = reads_cur['r1_start'] + read_len + jitter_value_vec
	read_1_df = reads_cur[['chr','r1_start_shifted', 'r1_end_shifted']]
	read_1_df['read_name'] = read_code_simu_cur_withUMI
	read_1_df['read_length'] = read_len
	read_1_df['strand'] = '+'
	return read_1_df


# def scRNA_SampleSyntheticReads(coordinate_file, samtools_directory, INPUT_bamfile, outdirectory, ref_peakfile, cellnumberfile):
# 	"""Sample Synthetic reads. The sampled reads coordinates are stored as `coordinate_file` in `outdirectory`.

# 	Parameters
# 	----------
# 	coordinate_file: `str`
# 		Specify name of the coordinates file storing synthetic reads. This is not a final output of scReadSim.
# 	samtools_directory: `str`
# 		Path to software samtools.
# 	INPUT_bamfile: `str`
# 		Path to input BAM file.
# 	outdirectory: `str`
# 		Output directory of coordinate files.
# 	ref_peakfile: `str`
# 		Directory of the features bed file.
# 	cellnumberfile: `str`
# 		Directory of the marginal synthetic count vector file output in scATAC_GenerateSyntheticCount step.
# 	"""
# 	rm_coor_command = "rm %s/%s" % (outdirectory, coordinate_file)
# 	os.system(rm_coor_command)
# 	create_coor_command = "touch %s/%s" % (outdirectory, coordinate_file)
# 	os.system(create_coor_command)
# 	cmd = "while true; do read -r region <&3 || break;  read -r ncell <&4 || break; region=$(echo ${region} | cut -f 1,2,3 | perl -lane 'print \"$F[0]:$F[1]-$F[2]\"');  paste -d\"\t\" <(awk -v nsample=${ncell} -v region=${region} 'BEGIN{for(c=0;c<nsample;c++) print region}') <(%s/samtools view %s ${region} | shuf -r -n ${ncell} | cut -f3,4) >> %s/%s;  done 3<%s/%s 4<%s" % (samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, outdirectory, ref_peakfile, cellnumberfile)
# 	# Testing using copied directory
# 	output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
# 	if error:
# 	     print('[ERROR] Fail to generate synthetic reads:\n', error.decode())


# def scRNA_GenerateBAMCoord(outdirectory, coordinate_file, assignment_file, count_mat_file, cellnumberfile, BED_filename, OUTPUT_cells_barcode_file):
# 	"""Generate Synthetic reads in BED format. 

# 	Parameters
# 	----------
# 	outdirectory: `str`
# 		Output directory of the synthteic bed file and its corresponding cell barcodes file.
# 	coordinate_file: `str`
# 		Directory of the coordinates file storing synthetic reads output by scATAC_SampleSyntheticReads.
# 	ref_peakfile: `str`
# 		Directory of the features bed file.
# 	count_mat_file: `str`
# 		Directory of synthetic count matrix.
# 	cellnumberfile: `str`
# 		Directory of the marginal synthetic count vector file output in scATAC_GenerateSyntheticCount step.
# 	BED_filename: `str`
# 		Specify the base name of output bed file.
# 	OUTPUT_cells_barcode_file: `str`
# 		Specify the file name storing the synthetic cell barcodes.
# 	read_len: `int` (default: '50')
# 		Specify the length of synthetic reads. Default value is 50 bp.
# 	jitter_size: `int` (default: '5')
# 		Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
# 	"""
# 	random.seed(2022)
# 	read_lines = pd.read_csv("%s/%s" % (outdirectory, coordinate_file), delimiter="\t",  names=['peak_name', 'chr', 'r1_start'])
# 	peaks = pd.read_csv("%s/%s" % (outdirectory, assignment_file), delimiter="\t",  names=['chr', 'start', 'end']).to_numpy()
# 	count_mat = pd.read_csv("%s/%s" % (outdirectory, count_mat_file), header=0, delimiter="\t").to_numpy()
# 	marginal_cell_number = pd.read_csv("%s" %  cellnumberfile, header=None, delimiter="\t").to_numpy()
# 	n_cell = np.shape(count_mat)[1]
# 	random_cellbarcode_list = cellbarcode_generator(n_cell, size=16)
# 	read_len = 80
# 	jitter_size = 5
# 	read1_bedfile="%s.read.bed" % BED_filename
# 	start = time.time()
# 	with open(OUTPUT_cells_barcode_file, 'w') as f:
# 		for item in random_cellbarcode_list:
# 		    f.write("%s\n" % item)
# 	with open("%s/%s" % (outdirectory, read1_bedfile), 'w') as fp:
# 		pass
# 	peak_nonzero_id = np.nonzero(marginal_cell_number)[0]
# 	for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
# 		peak_ind = peak_nonzero_id[relative_peak_ind]
# 		# peak_ind = 192
# 		peak_record = peaks[peak_ind]
# 		count_vec = count_mat[peak_ind,:] # Input
# 		print(peak_ind)
# 		read_1_df = scRNA_PerTruePeakEdition(peak_record, count_vec, read_lines, read_len, jitter_size, random_cellbarcode_list)
# 		read_1_df.to_csv("%s/%s" % (outdirectory, read1_bedfile), header=None, index=None, sep='\t', mode='a')
# 	# for peak_ind in tqdm(range(len(peaks))):
# 	# 	# peak_ind = 192
# 	# 	peak_record = peaks[peak_ind]
# 	# 	count_vec = count_mat[peak_ind,:] # Input
# 	# 	print(peak_ind)
# 	# 	read_1_df = scRNA_PerTruePeakEdition(peak_record, count_vec, read_lines, read_len, jitter_size, random_cellbarcode_list)
# 	# 	read_1_df.to_csv("%s/%s" % (outdirectory, read1_bedfile), header=None, index=None, sep='\t', mode='a')
# 	end = time.time()
# 	print(end - start)

def scRNA_GenerateBAMCoord(count_mat_filename, samtools_directory, INPUT_bamfile, ref_peakfile, directory_cellnumber, outdirectory, BED_filename, OUTPUT_cells_barcode_file, read_len=80, jitter_size=5):
	"""Generate Synthetic reads in BED format. 

	Parameters
	----------
	count_mat_filename: `str`
		The base name of output count matrix in bam2countmat.
	samtools_directory: `str`
		Path to software samtools.
	INPUT_bamfile: `str`
		Input BAM file for anlaysis.
	ref_peakfile: `str`
		Features bed file.
	directory_cellnumber: `str`
		Directory of the marginal synthetic count vector file output in scATAC_GenerateSyntheticCount step.
	outdirectory: `str`
		Specify the output directory for synthetic reads bed file.
	BED_filename: `str`
		Specify the base name of output bed file.
	OUTPUT_cells_barcode_file: `str`
		Specify the file name storing the synthetic cell barcodes.
	read_len: `int` (default: '80')
		Specify the length of synthetic reads. Default value is 50 bp.
	jitter_size: `int` (default: '5')
		Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
	"""
	# Generating BAMfile_coordinates file
	print('scReadSim scRNA_GenerateBAMCoord Running...')
	print('\t- Sampling synthetic reads...')
	coordinate_file = count_mat_filename + ".BAMfile_coordinates.txt"
	cellnumberfile = "%s/%s.scDesign2Simulated.nReadRegionmargional.txt" % (directory_cellnumber, count_mat_filename)
	rm_coor_command = "rm %s/%s" % (outdirectory, coordinate_file)
	os.system(rm_coor_command)
	create_coor_command = "touch %s/%s" % (outdirectory, coordinate_file)
	os.system(create_coor_command)
	cmd = "while true; do read -r region <&3 || break;  read -r ncell <&4 || break; region=$(echo ${region} | cut -f 1,2,3 | perl -lane 'print \"$F[0]:$F[1]-$F[2]\"');  paste -d\"\t\" <(awk -v nsample=${ncell} -v region=${region} 'BEGIN{for(c=0;c<nsample;c++) print region}') <(%s/samtools view %s ${region} | shuf -r -n ${ncell} | cut -f3,4) >> %s/%s;  done 3<%s 4<%s" % (samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, ref_peakfile, cellnumberfile)
	output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to generate synthetic reads:\n', error.decode())
	print('\t- Sampling synthetic reads done!')
	# Generating BED file
	print('\t- Writing synthetic reads...')
	count_mat_file = "%s.scDesign2Simulated.txt" % count_mat_filename
	random.seed(2022)
	read_lines = pd.read_csv("%s/%s" % (outdirectory, coordinate_file), delimiter="\t",  names=['peak_name', 'chr', 'r1_start'])
	peaks = pd.read_csv(ref_peakfile, delimiter="\t",  names=['chr', 'start', 'end']).to_numpy()
	count_mat = pd.read_csv("%s/%s" % (directory_cellnumber, count_mat_file), header=0, delimiter="\t").to_numpy()
	marginal_cell_number = pd.read_csv(cellnumberfile, header=None, delimiter="\t").to_numpy()
	n_cell = np.shape(count_mat)[1]
	random_cellbarcode_list = cellbarcode_generator(n_cell, size=16)
	read1_bedfile="%s.read.bed" % BED_filename
	with open(outdirectory + "/" + OUTPUT_cells_barcode_file, 'w') as f:
		for item in random_cellbarcode_list:
		    f.write("%s\n" % item)
	with open("%s/%s" % (outdirectory, read1_bedfile), 'w') as fp:
		pass
	peak_nonzero_id = np.nonzero(marginal_cell_number)[0]
	for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
		peak_ind = peak_nonzero_id[relative_peak_ind]
		# peak_ind = 192
		peak_record = peaks[peak_ind]
		count_vec = count_mat[peak_ind,:] # Input
		print(peak_ind)
		read_1_df = scRNA_PerTruePeakEdition(peak_record, count_vec, read_lines, random_cellbarcode_list, read_len, jitter_size)
		read_1_df.to_csv("%s/%s" % (outdirectory, read1_bedfile), header=None, index=None, sep='\t', mode='a')
	# for peak_ind in tqdm(range(len(peaks))):
	# 	# peak_ind = 192
	# 	peak_record = peaks[peak_ind]
	# 	count_vec = count_mat[peak_ind,:] # Input
	# 	print(peak_ind)
	# 	read_1_df = scRNA_PerTruePeakEdition(peak_record, count_vec, read_lines, read_len, jitter_size, random_cellbarcode_list)
	# 	read_1_df.to_csv("%s/%s" % (outdirectory, read1_bedfile), header=None, index=None, sep='\t', mode='a')
	print('\t- Writing synthetic reads done.')
	print('Done!\n')


def scRNA_CombineBED(outdirectory, BED_filename_pre, BED_COMPLE_filename_pre, BED_filename_combined_pre):
	"""Combine the bed files of foreground and background feature sets into one bed file.

	"""
	combine_read_cmd = "cat %s/%s.read.bed %s/%s.read.bed | sort -k1,1 -k2,2n | cut -f1-5 > %s/%s.read.bed" % (outdirectory, BED_filename_pre, outdirectory, BED_COMPLE_filename_pre, outdirectory, BED_filename_combined_pre)
	output, error = subprocess.Popen(combine_read_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
 		print('[ERROR] Fail to create combine synthetic read bed files:\n', error.decode())

def scRNA_BED2FASTQ(bedtools_directory, seqtk_directory, referenceGenome_file, outdirectory, BED_filename_combined, synthetic_fastq_prename, sort_FASTQ = True):
	"""Convert Synthetic reads from BED to FASTQ. 

	Parameters
	----------
	bedtools_directory: `str`
		Path to software bedtools.
	seqtk_directory: `str`
		Path to software seqtk.
	referenceGenome_file: `str`
		Reference genome FASTA file that the synthteic reads should align.
	outdirectory: `str`
		Output directory of the synthteic bed file and its corresponding cell barcodes file.
	BED_filename_combined: `str`
		Base name of the combined bed file output by function `scRNA_CombineBED`.
	synthetic_fastq_prename
		Specify the base name of the output FASTQ files.
	sort_FASTQ: `bool` (Default: True)
		Set `True` to sort the output FASTQ file.
	"""
	# Create FASTA
	print('scReadSim BED2FASTQ_Pair Running...')
	print('\t- Creating FASTA files...')
	fasta_read2_cmd = "%s/bedtools getfasta -fi %s -bed %s/%s.read.bed -fo %s/%s.read2.bed2fa.fa -nameOnly" % (bedtools_directory, referenceGenome_file, outdirectory, BED_filename_combined, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(fasta_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print(error.decode())
	fasta_read1_cmd = "awk 'NR%%2==0 {print substr(p,2,26);} NR%%2 {p=$0;print p;}' %s/%s.read2.bed2fa.fa > %s/%s.read1.bed2fa.fa" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(fasta_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print(error.decode())
	print('\t- Converting FASTA files to FASTQ files...')
	# FASTA to FASTQ
	fastq_read1_cmd = "%s/seqtk seq -F 'F' %s/%s.read1.bed2fa.fa > %s/%s.read1.bed2fa.fq" % (seqtk_directory, outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(fastq_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to convert read1 synthetic fasta file to fastq file:', error.decode())
	fastq_read2_cmd = "%s/seqtk seq -F 'F' %s/%s.read2.bed2fa.fa > %s/%s.read2.bed2fa.fq" % (seqtk_directory, outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(fastq_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to convert read2 synthetic fasta file to fastq file:', error.decode())
	print('\t- FASTQ files %s.read1.bed2fa.fq, %s.read2.bed2fa.fq stored in %s.' % (synthetic_fastq_prename, synthetic_fastq_prename, outdirectory))
	if sort_FASTQ == True:
		print('\t- Sorting FASTQ files...')
		sort_fastq_read1_cmd = "cat %s/%s.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.read1.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
		output, error = subprocess.Popen(sort_fastq_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if error:
		     print('[ERROR] Fail to sort read1 synthetic fastq file:', error.decode())
		sort_fastq_read2_cmd = "cat %s/%s.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.read2.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
		output, error = subprocess.Popen(sort_fastq_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if error:
		     print('[ERROR] Fail to sort read2 synthetic fastq file:', error.decode())
		print('\t- Sorted FASTQ files %s.read1.bed2fa.sorted.fq, %s.read2.bed2fa.sorted.fq stored in %s.' % (synthetic_fastq_prename, synthetic_fastq_prename, outdirectory))
	print('Done!')


def AlignSyntheticBam_Pair(bowtie2_directory, samtools_directory, outdirectory, referenceGenome_name, referenceGenome_dir, synthetic_fastq_prename, output_BAM_pre):
	"""Convert Synthetic reads from FASTQ to BAM.

	Parameters
	----------
	bowtie2_directory: `str`
		Path to software bowtie2.
	samtools_directory: `str`
		Path to software samtools.
	outdirectory: `str`
		Specify the output directory of the synthteic BAM file.
	referenceGenome_name: `str`
		Base name of the eference genome FASTA file. For example, you should input "chr1" for file "chr1.fa".
	referenceGenome_dir: `str`
		Path to the reference genome FASTA file.
	synthetic_fastq_prename: `str`
		Base name of the synthetic FASTQ files output by function `scRNA_BED2FASTQ`.
	output_BAM_pre: `str`
		Specify the base name of the output BAM file.
	"""
	# print('scReadSim AlignSyntheticBam_Pair Running...')
	# if doIndex == True:
	# 	print('\t- Indexing reference genome file...')
	# 	index_ref_cmd = "%s/bowtie2-build %s/%s.fa %s" % (bowtie2_directory, referenceGenome_dir, referenceGenome_name, referenceGenome_name)
	# 	output, error = subprocess.Popen(index_ref_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	# 	if error:
	# 		print(error.decode())
	# 		# print('[ERROR] Fail to index the reference genome:\nPlease index the reference genome by setting \'doIndex=True\'\n', error.decode())
	# Align using bwa
	print('\t- Aligning FASTQ files onto reference genome files...')
	alignment_cmd = "%s/bowtie2 -x %s/%s -1 %s/%s.read1.bed2fa.sorted.fq -2 %s/%s.read2.bed2fa.sorted.fq | %s/samtools view -bS - > %s/%s.synthetic.noCB.bam" % (bowtie2_directory, referenceGenome_dir, referenceGenome_name,  outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename, samtools_directory, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(alignment_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print(error.decode())
	print('\t- Alignment Done.')
	print('\t- Generating cell barcode tag...')
	addBC2BAM_header_cmd = "%s/samtools view %s/%s.synthetic.noCB.bam -H > %s/%s.synthetic.noCB.header.sam" % (samtools_directory, outdirectory, output_BAM_pre, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(addBC2BAM_header_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	addBC2BAM_cmd = "cat <( cat %s/%s.synthetic.noCB.header.sam ) <( paste <(%s/samtools view %s/%s.synthetic.noCB.bam ) <(%s/samtools view %s/%s.synthetic.noCB.bam | cut -f1 | cut -d':' -f1 | sed -e 's/^/CB:Z:/')) | %s/samtools view -bS - > %s/%s.synthetic.bam" % (outdirectory, output_BAM_pre, samtools_directory, outdirectory, output_BAM_pre, samtools_directory, outdirectory, output_BAM_pre, samtools_directory, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(addBC2BAM_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to add BC tag to synthetic BAM file:', error.decode())
	print('\t- Generating cell barcode tag done.')
	print('\t- Sorting synthetic BAM file...')
	sortBAMcmd = "%s/samtools sort %s/%s.synthetic.bam > %s/%s.synthetic.sorted.bam" % (samtools_directory, outdirectory, output_BAM_pre, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(sortBAMcmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to sort synthetic BAM file:', error.decode())
	print('\t- Indexing synthetic BAM file...')
	indexBAMcmd = "%s/samtools index %s/%s.synthetic.sorted.bam" % (samtools_directory, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(indexBAMcmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to index synthetic BAM file:', error.decode())
	print('Done!\n')


## Error rate
def ErrorBase(base, prop, base_call_ref):
	"""Sample random errors. 

	"""
	err_base_call_id = np.random.choice(a=[0, 1, 2], size=1, p=prop)[0]      
	err_base_call = base_call_ref[base][err_base_call_id]
	return err_base_call


def ErroneousRead(real_error_rate_read, read_df, output_fq_file):
	"""Generate random errors according to input real data error rates. 

	"""
	n_read = int(np.shape(read_df)[0]/4)
	## Prepare Error rate
	real_error_rate_read_A = real_error_rate_read[['a_to_c_error_rate', 'a_to_g_error_rate', 'a_to_t_error_rate']].to_numpy() 
	real_error_rate_read_A_prop = real_error_rate_read_A/real_error_rate_read_A.sum(axis=1,keepdims=1)
	real_error_rate_read_A_prop[np.isnan(real_error_rate_read_A_prop).any(axis=1),:] = 1/3
	real_error_rate_read_C = real_error_rate_read[['c_to_a_error_rate', 'c_to_g_error_rate', 'c_to_t_error_rate']].to_numpy() 
	real_error_rate_read_C_prop = real_error_rate_read_C/real_error_rate_read_C.sum(axis=1,keepdims=1)
	real_error_rate_read_C_prop[np.isnan(real_error_rate_read_C_prop).any(axis=1),:] = 1/3
	real_error_rate_read_G = real_error_rate_read[['g_to_a_error_rate', 'g_to_c_error_rate', 'g_to_t_error_rate']].to_numpy() 
	real_error_rate_read_G_prop = real_error_rate_read_G/real_error_rate_read_G.sum(axis=1,keepdims=1)
	real_error_rate_read_G_prop[np.isnan(real_error_rate_read_G_prop).any(axis=1),:] = 1/3
	real_error_rate_read_T = real_error_rate_read[['t_to_a_error_rate', 't_to_c_error_rate', 't_to_g_error_rate']].to_numpy() 
	real_error_rate_read_T_prop = real_error_rate_read_T/real_error_rate_read_T.sum(axis=1,keepdims=1)
	real_error_rate_read_T_prop[np.isnan(real_error_rate_read_T_prop).any(axis=1),:] = 1/3
	# Base decision matrix
	real_error_rate_read_prop_dict = {'A': real_error_rate_read_A_prop, 'C': real_error_rate_read_C_prop, 'G': real_error_rate_read_G_prop, 'T': real_error_rate_read_T_prop}
	# Error decision vector
	real_error_rate_read_perbase = real_error_rate_read['error_rate'].to_numpy()
	## Decide whether error occurs for each read
	read_length = real_error_rate_read.shape[0]
	error_read_perbase_indicator = np.zeros((n_read, read_length), dtype=int)
	random.seed(1)
	for base_id in range(read_length):
		error_read_perbase_indicator[:,base_id] = np.random.binomial(n=1, p=real_error_rate_read_perbase[base_id], size=n_read)
	erroneous_read_id = np.where(np.sum(error_read_perbase_indicator, axis=1) > 0)[0]
	## For erroneous reads, generate erroneous base based on the probability matrix
	base_call_ref = {'A': ['C', 'G', 'T'], 'C': ['A', 'G', 'T'], 'G': ['A', 'C', 'T'], 'T': ['A', 'C', 'G']}
	random.seed(2021)
	read_df_witherror = read_df
	for read_id_tqdm in tqdm(range(len(erroneous_read_id))):
		read_id = erroneous_read_id[read_id_tqdm]
		read_cur = read_df[(read_id*4) : (read_id*4 + 4)]
		bases = list(read_cur[1][0].upper())
		Qscores = list(read_cur[3][0])
		for errorneous_base_id in np.where(error_read_perbase_indicator[read_id,:] > 0)[0]:
			if errorneous_base_id < len(bases):
				base_cur = bases[errorneous_base_id]
				if base_cur in real_error_rate_read_prop_dict:
					prop = real_error_rate_read_prop_dict[base_cur][errorneous_base_id]
					# Decide error base
					err_base_call = ErrorBase(base_cur, prop, base_call_ref)
					bases[errorneous_base_id] = err_base_call
					# Decide Q score
					# Use 9 (Phred score 24) for erroneous base
					Qscores[errorneous_base_id] = '9'
		read_df_witherror[read_id*4+1] = ''.join(bases)
		read_df_witherror[read_id*4+3] = ''.join(Qscores)
	## Write out
	np.savetxt(output_fq_file, read_df_witherror, fmt='%s')


def SubstiError(real_error_rate_file, outdirectory, synthetic_fastq_prename):
	"""Generate random errors for single-end sequencing reads according to input real data error rates. 

	Parameters
	----------
	real_error_rate_file: `str`
		Path to software fgbio jar script.
	outdirectory: `str`
		Specify the output directory of the synthteic FASTQ file with random errors.
	synthetic_fastq_prename: `str`
		Specify the base name of the synthetic erroneous reads' FASTQ files.
	"""
	# Read in real error rates
	real_error_rate_dir = real_error_rate_file
	real_error_rate = pd.read_csv(real_error_rate_dir, header=0, delimiter="\t")
	# Read in perfect reads
	read2_fq = outdirectory  + "/" + synthetic_fastq_prename + ".read2.bed2fa.fq"
	read2_df = pd.read_csv(read2_fq, header=None).to_numpy()
	# Real data quality score
	# Error rate to Qscore
	# err_rate_qscore_read1 = -10 * np.log10(real_error_rate_read1['error_rate'])
	# fastqc_result = "/home/gayan/Projects/scATAC_Simulator/results/20220408_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_1_5000000_NONINPUT/VerifyQuality/fastqc_Real/10X_ATAC_chr1_1_5000000_fastqc/fastqc_data.txt"
	# with open(fastqc_result) as f:
	#     fastqc_result_contents = f.read().splitlines()
	# fastqc_result_contents = np.array(fastqc_result_contents)
	# real_quality_range = np.where(fastqc_result_contents == '>>END_MODULE')[0][0:2]
	# real_quality = list(fastqc_result_contents[(real_quality_range[0]+2) : real_quality_range[1]])
	# real_quality_array = np.array([i.split('\t') for i in real_quality])
	ErroneousRead(real_error_rate, read2_df, outdirectory + "/" + synthetic_fastq_prename + ".ErrorIncluded.read2.bed2fa.fq") 


def scRNA_ErrorBase(fgbio_jarfile, INPUT_bamfile, referenceGenome_file, outdirectory, synthetic_fastq_prename):
	"""Introduce random substitution errors into synthetic reads according to real data error rates.

	Parameters
	----------
	fgbio_jarfile: `str`
		Path to software fgbio jar script.
	INPUT_bamfile: `str`
		Input BAM file for anlaysis.
	referenceGenome_file: 'str'
		Reference genome FASTA file that the synthteic reads should align.
	outdirectory: `str`
		Specify the output directory of the synthteic FASTQ file with random errors.
	synthetic_fastq_prename: `str`
		Base name of the synthetic FASTQ files output by function `scATAC_BED2FASTQ`.
	"""
	combine_read1_cmd = "java -jar %s ErrorRateByReadPosition -i %s -r %s -o %s/Real --collapse false" % (fgbio_jarfile, INPUT_bamfile, referenceGenome_file, outdirectory)
	output, error = subprocess.Popen(combine_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[Messages] Running fgbio on real bam file:\n', error.decode())
	# Generate Errors into fastq files
	real_error_rate_file = outdirectory + "/" + "Real.error_rate_by_read_position.txt"
	SubstiError(real_error_rate_file, outdirectory, synthetic_fastq_prename)
	# Combine FASTQs
	print('\t- Sorting FASTQ files...')
	sort_fastq_read1_cmd = "cat %s/%s.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.read1.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(sort_fastq_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to sort read1 synthetic fastq file:', error.decode())
	sort_fastq_read2_cmd = "cat %s/%s.ErrorIncluded.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.ErrorIncluded.read2.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(sort_fastq_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to sort read2 synthetic fastq file:', error.decode())
	print('\t- Sorted FASTQ files %s.read1.bed2fa.sorted.fq, %s.ErrorIncluded.read2.bed2fa.sorted.fq stored in %s.' % (synthetic_fastq_prename, synthetic_fastq_prename, outdirectory))



