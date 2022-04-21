import pandas as pd
import pickle
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


def cellbarode_generator(length, size=10):
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


def scATAC_INPUT_PerTruePeakEdition(peak_record, count_vec, read_lines, random_cellbarcode_list, read_len=50, jitter_size=5):
	"""Formulate Synthetic reads for task with user input features set.

    Parameters
    ----------
	peak_record: `numpy.array`
        Coordinates of an input foreground(/background) peak (columns 1:3) and its assigned reference peak coordinates (columns 4:6).
	count_vec: `numpy.array`
		Count vector of the input foreground(/background) peak.
	read_lines: `pandas.dataframe`
		Coordinates of synthetic reads sampled from the input BAM file.
	random_cellbarcode_list: `list`
		List of cell barcodes randomly generated using cellbarode_generator().
	read_len: `int` (default: '50')
		Specify the length of synthetic reads. Default value is 50 bp.
	jitter_size: `int` (default: '50')
		Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.

    Return
    ------
    read_1_df: `pandas.dataframe`
        Coordinates of synthetic read 1 generated for the input foreground(/background) peak.
    read_2_df: `pandas.dataframe`
        Coordinates of synthetic read 2 generated for the input foreground(/background) peak.
    """
	ref_peak_concat = peak_record[3] + ":" + str(peak_record[4]) + "-" + str(peak_record[5])
	true_peak_concat = peak_record[0] + ":" + str(peak_record[1]) + "-" + str(peak_record[2])
	reads_cur = read_lines[read_lines['true_peak_name'] == true_peak_concat] # Input
	nread_cur= np.sum(count_vec)
	count_frag_vec = np.ceil(count_vec/2).astype(int)
	nfrag_cur= np.sum(count_frag_vec).astype(int) # nrow(reads_cur) should equal to nfrag_cur
	shift_number = peak_record[1] - peak_record[4] # ref + shift_number = true position
	# Add cell information
	nonempty_cell_ind = np.where(count_frag_vec != 0)[0]
	read_code_simu_cur = [random_cellbarcode_list[nonempty_cell_ind[ind]] + ":CellType1" + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(true_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_frag_vec[nonempty_cell_ind[ind]])]
	jitter_value_vec = np.random.random_integers(-jitter_size,jitter_size,size=np.shape(reads_cur)[0])  # nrow(reads_cur) should equal to nfrag_cur
	contain_read_indicator = reads_cur['r1_start'] == reads_cur['r2_start']
	reads_cur['read_length'] = read_len
	reads_cur['read_length'][contain_read_indicator] = abs(reads_cur['length'].astype(int)[contain_read_indicator])
	reads_cur['r1_start_shifted'] = reads_cur['r1_start'].astype(int)  + shift_number + jitter_value_vec
	reads_cur['r2_start_shifted'] = reads_cur['r2_start'].astype(int)  + shift_number + jitter_value_vec
	reads_cur['r1_end_shifted'] = reads_cur['r1_start'].astype(int) + reads_cur['read_length'].astype(int) + shift_number + jitter_value_vec
	reads_cur['r2_end_shifted'] = reads_cur['r2_start'].astype(int) + reads_cur['read_length'].astype(int) + shift_number + jitter_value_vec
	# Create dataframes for paired end reads
	read_1_df = pd.concat([reads_cur.loc[reads_cur['length'] >= 0, ['chr','r1_start_shifted', 'r1_end_shifted', 'read_length']], reads_cur.loc[reads_cur['length'] < 0, ['chr','r2_start_shifted', 'r2_end_shifted', 'read_length']].rename(columns={'r2_start_shifted':'r1_start_shifted', 'r2_end_shifted':'r1_end_shifted'})], ignore_index=True)
	read_2_df = pd.concat([reads_cur.loc[reads_cur['length'] >= 0, ['chr','r2_start_shifted', 'r2_end_shifted', 'read_length']], reads_cur.loc[reads_cur['length'] < 0, ['chr','r1_start_shifted', 'r1_end_shifted', 'read_length']].rename(columns={'r1_start_shifted':'r2_start_shifted', 'r1_end_shifted':'r2_end_shifted'})], ignore_index=True)
	# Fill in other information
	read_1_df['read_name'] = read_code_simu_cur
	read_2_df['read_name'] = read_code_simu_cur
	# read_1_df['read_length'] = read_len
	# read_2_df['read_length'] = read_len
	read_1_df['strand'] = '+'
	read_2_df['strand'] = '-'
	read_1_df_order = read_1_df[['chr','r1_start_shifted', 'r1_end_shifted', 'read_name', 'read_length', 'strand']]
	read_2_df_order = read_2_df[['chr','r2_start_shifted', 'r2_end_shifted', 'read_name', 'read_length', 'strand']]
	return read_1_df_order, read_2_df_order

def scATAC_PerTruePeakEdition(peak_record, count_vec, read_lines, random_cellbarcode_list, read_len=50, jitter_size=5):
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
    	List of cell barcodes randomly generated using cellbarode_generator().
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
	true_peak_concat = peak_record[0] + ":" + str(peak_record[1]) + "-" + str(peak_record[2])
	reads_cur = read_lines[read_lines['peak_name'] == true_peak_concat] # Input
	nread_cur= np.sum(count_vec)
	count_frag_vec = np.ceil(count_vec/2).astype(int)
	nfrag_cur= np.sum(count_frag_vec).astype(int) # nrow(reads_cur) should equal to nfrag_cur
	# Add cell information
	nonempty_cell_ind = np.where(count_frag_vec != 0)[0]
	read_code_simu_cur = [random_cellbarcode_list[nonempty_cell_ind[ind]] + ":CellType1" + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(true_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_frag_vec[nonempty_cell_ind[ind]])]
	# start = time.time()
	jitter_value_vec = np.random.random_integers(-jitter_size,jitter_size,size=np.shape(reads_cur)[0])  # nrow(reads_cur) should equal to nfrag_cur
	contain_read_indicator = reads_cur['r1_start'] == reads_cur['r2_start']
	reads_cur['read_length'] = read_len
	reads_cur['read_length'][contain_read_indicator] = abs(reads_cur['length'].astype(int)[contain_read_indicator])
	reads_cur['r1_start_shifted'] = reads_cur['r1_start'].astype(int)  + jitter_value_vec
	reads_cur['r2_start_shifted'] = reads_cur['r2_start'].astype(int)  + jitter_value_vec
	reads_cur['r1_end_shifted'] = reads_cur['r1_start'].astype(int) + reads_cur['read_length'].astype(int) + jitter_value_vec
	reads_cur['r2_end_shifted'] = reads_cur['r2_start'].astype(int) + reads_cur['read_length'].astype(int) + jitter_value_vec
	# Create dataframes for paired end reads
	read_1_df = pd.concat([reads_cur.loc[reads_cur['length'] >= 0, ['chr','r1_start_shifted', 'r1_end_shifted', 'read_length']], reads_cur.loc[reads_cur['length'] < 0, ['chr','r2_start_shifted', 'r2_end_shifted', 'read_length']].rename(columns={'r2_start_shifted':'r1_start_shifted', 'r2_end_shifted':'r1_end_shifted'})], ignore_index=True)
	read_2_df = pd.concat([reads_cur.loc[reads_cur['length'] >= 0, ['chr','r2_start_shifted', 'r2_end_shifted', 'read_length']], reads_cur.loc[reads_cur['length'] < 0, ['chr','r1_start_shifted', 'r1_end_shifted', 'read_length']].rename(columns={'r1_start_shifted':'r2_start_shifted', 'r1_end_shifted':'r2_end_shifted'})], ignore_index=True)
	# Fill in other information
	read_1_df['read_name'] = read_code_simu_cur
	read_2_df['read_name'] = read_code_simu_cur
	# read_1_df['read_length'] = read_len
	# read_2_df['read_length'] = read_len
	read_1_df['strand'] = '+'
	read_2_df['strand'] = '-'
	read_1_df_order = read_1_df[['chr','r1_start_shifted', 'r1_end_shifted', 'read_name', 'read_length', 'strand']]
	read_2_df_order = read_2_df[['chr','r2_start_shifted', 'r2_end_shifted', 'read_name', 'read_length', 'strand']]
	return read_1_df_order, read_2_df_order

# def scATAC_SampleSyntheticReads(coordinate_file, samtools_directory, INPUT_bamfile, outdirectory, ref_peakfile, cellnumberfile):
# 	"""Sample Synthetic reads. The sampled reads coordinates are stored as `coordinate_file` in `outdirectory`.

#     Parameters
#     ----------
# 	coordinate_file: `str`
# 		Specify name of the coordinates file storing synthetic reads. This is not a final output of scReadSim.
# 	samtools_directory: `str`
# 		Directory of software samtools.
# 	INPUT_bamfile: `str`
# 		Directory of input BAM file.
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
# 	cmd = "while true; do read -r region <&3 || break;  read -r ncell <&4 || break; region=$(echo ${region} | cut -f 1,2,3 | perl -lane 'print \"$F[0]:$F[1]-$F[2]\"');  paste -d\"\t\" <(awk -v nsample=${ncell} -v region=${region} 'BEGIN{for(c=0;c<nsample;c++) print region}') <(%s/samtools view %s ${region} | shuf -r -n ${ncell} | cut -f3,4,8,9) >> %s/%s;  done 3<%s 4<%s" % (samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, ref_peakfile, cellnumberfile)
# 	# Testing using copied directory
# 	output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
# 	if error:
# 	     print('[ERROR] Fail to generate synthetic reads:\n', error.decode())


# def scATAC_GenerateBAMCoord(outdirectory, coordinate_file, ref_peakfile, count_mat_file, cellnumberfile, BED_filename, OUTPUT_cells_barcode_file, read_len = 50, jitter_size = 5):
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
# 	read_lines = pd.read_csv(coordinate_file, delimiter="\t",  names=['peak_name', 'chr', 'r1_start', 'r2_start', 'length'])
# 	peaks_assignments = pd.read_csv(ref_peakfile, delimiter="\t",  names=['chr', 'start', 'end']).to_numpy()
# 	count_mat = pd.read_csv(count_mat_file, header=0, delimiter="\t").to_numpy()
# 	marginal_cell_number = pd.read_csv(cellnumberfile, header=None, delimiter="\t").to_numpy()
# 	n_cell = np.shape(count_mat)[1]
# 	random_cellbarcode_list = cellbarode_generator(n_cell, size=16)
# 	read1_bedfile="%s.read1.bed" % BED_filename
# 	read2_bedfile="%s.read2.bed" % BED_filename
# 	start = time.time()
# 	with open(outdirectory + "/" + OUTPUT_cells_barcode_file, 'w') as f:
# 		for item in random_cellbarcode_list:
# 		    f.write("%s\n" % item)
# 	with open("%s/%s" % (outdirectory, read1_bedfile), 'w') as fp:
# 		pass
# 	with open("%s/%s" % (outdirectory, read2_bedfile), 'w') as fp:
# 		pass
# 	peak_nonzero_id = np.nonzero(marginal_cell_number)[0]
# 	for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
# 		peak_ind = peak_nonzero_id[relative_peak_ind]
# 		# peak_ind = 192
# 		peak_record = peaks_assignments[peak_ind]
# 		count_vec = count_mat[peak_ind,:] # Input
# 		print(peak_ind)
# 		read_1_df, read_2_df = scATAC_PerTruePeakEdition(peak_record, count_vec, read_lines, random_cellbarcode_list, read_len, jitter_size)
# 		read_1_df.to_csv("%s/%s" % (outdirectory, read1_bedfile), header=None, index=None, sep='\t', mode='a')
# 		read_2_df.to_csv("%s/%s" % (outdirectory, read2_bedfile), header=None, index=None, sep='\t', mode='a')
# 	end = time.time()
# 	print(end - start)

def scATAC_GenerateBAMCoord(count_mat_filename, samtools_directory, INPUT_bamfile, ref_peakfile, directory_cellnumber, outdirectory, BED_filename, OUTPUT_cells_barcode_file, read_len = 50, jitter_size = 5):
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
	read_len: `int` (default: '50')
		Specify the length of synthetic reads. Default value is 50 bp.
	jitter_size: `int` (default: '5')
		Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
	"""
	print('scReadSim scRNA_GenerateBAMCoord Running...')
	print('\t- Sampling synthetic reads...')
	coordinate_file = count_mat_filename + ".BAMfile_halfsampled_coordinates.txt"
	cellnumberfile = "%s/%s.scDesign2Simulated.nPairsRegionmargional.txt" % (directory_cellnumber, count_mat_filename)
	rm_coor_command = "rm %s/%s" % (outdirectory, coordinate_file)
	os.system(rm_coor_command)
	create_coor_command = "touch %s/%s" % (outdirectory, coordinate_file)
	os.system(create_coor_command)
	cmd = "while true; do read -r region <&3 || break;  read -r ncell <&4 || break; region=$(echo ${region} | cut -f 1,2,3 | perl -lane 'print \"$F[0]:$F[1]-$F[2]\"');  paste -d\"\t\" <(awk -v nsample=${ncell} -v region=${region} 'BEGIN{for(c=0;c<nsample;c++) print region}') <(%s/samtools view %s ${region} | shuf -r -n ${ncell} | cut -f3,4,8,9) >> %s/%s;  done 3<%s 4<%s" % (samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, ref_peakfile, cellnumberfile)
	# Testing using copied directory
	output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to generate synthetic reads:\n', error.decode())
	print('\t- Sampling synthetic reads done!')
	count_mat_file = "%s.scDesign2Simulated.txt" % count_mat_filename
	random.seed(2022)
	read_lines = pd.read_csv("%s/%s" % (outdirectory, coordinate_file), delimiter="\t",  names=['peak_name', 'chr', 'r1_start', 'r2_start', 'length'])
	peaks_assignments = pd.read_csv(ref_peakfile, delimiter="\t",  names=['chr', 'start', 'end']).to_numpy()
	count_mat = pd.read_csv("%s/%s" % (directory_cellnumber, count_mat_file), header=0, delimiter="\t").to_numpy()
	marginal_cell_number = pd.read_csv(cellnumberfile, header=None, delimiter="\t").to_numpy()
	n_cell = np.shape(count_mat)[1]
	random_cellbarcode_list = cellbarode_generator(n_cell, size=16)
	read1_bedfile="%s.read1.bed" % BED_filename
	read2_bedfile="%s.read2.bed" % BED_filename
	start = time.time()
	with open(outdirectory + "/" + OUTPUT_cells_barcode_file, 'w') as f:
		for item in random_cellbarcode_list:
		    f.write("%s\n" % item)
	with open("%s/%s" % (outdirectory, read1_bedfile), 'w') as fp:
		pass
	with open("%s/%s" % (outdirectory, read2_bedfile), 'w') as fp:
		pass
	peak_nonzero_id = np.nonzero(marginal_cell_number)[0]
	for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
		peak_ind = peak_nonzero_id[relative_peak_ind]
		# peak_ind = 192
		peak_record = peaks_assignments[peak_ind]
		count_vec = count_mat[peak_ind,:] # Input
		print(peak_ind)
		read_1_df, read_2_df = scATAC_PerTruePeakEdition(peak_record, count_vec, read_lines, random_cellbarcode_list, read_len, jitter_size)
		read_1_df.to_csv("%s/%s" % (outdirectory, read1_bedfile), header=None, index=None, sep='\t', mode='a')
		read_2_df.to_csv("%s/%s" % (outdirectory, read2_bedfile), header=None, index=None, sep='\t', mode='a')
	end = time.time()
	print(end - start)
 

# def scATAC_SampleSyntheticReads_INPUT(coordinate_file, samtools_directory, INPUT_bamfile, outdirectory, assignment_file, cellnumberfile):
# 	"""Sample Synthetic reads. The sampled reads coordinates are stored as `coordinate_file` in `outdirectory`.

# 	Parameters
# 	----------
# 	coordinate_file: `str`
# 		Specify name of the coordinates file storing synthetic reads. This is not a final output of scReadSim.
# 	samtools_directory: `str`
# 		Directory of software samtools.
# 	INPUT_bamfile: `str`
# 		Directory of input BAM file.
# 	outdirectory: `str`
# 		Output directory of coordinate files.
# 	assignment_file: `str`
# 		Directory of the peak assignment file generated in match_peak.
# 	cellnumberfile: `str`
# 		Directory of the marginal synthetic count vector file output in scATAC_GenerateSyntheticCount step.
# 	"""
# 	rm_coor_command = "rm %s/%s" % (outdirectory, coordinate_file)
# 	os.system(rm_coor_command)
# 	create_coor_command = "touch %s/%s" % (outdirectory, coordinate_file)
# 	os.system(create_coor_command)
# 	cmd = "while true; do read -r region <&3 || break;  read -r ncell <&4 || break; true_region=$(echo ${region} | cut -f 1,2,3 | perl -lane 'print \"$F[0]:$F[1]-$F[2]\"'); ref_region=$(echo ${region} | cut -f 4,5,6 | perl -lane 'print \"$F[0]:$F[1]-$F[2]\"'); paste -d\"\t\" <(awk -v nsample=${ncell} -v region=${true_region} 'BEGIN{for(c=0;c<nsample;c++) print region}') <(%s/samtools view %s ${ref_region} | shuf -r -n ${ncell} | cut -f3,4,8,9) >> %s/%s;  done 3<%s 4<%s" % (samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, assignment_file, cellnumberfile)
# 	# Testing using copied directory
# 	output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
# 	if error:
# 	     print('[ERROR] Fail to generate synthetic reads:\n', error.decode())
# 	print('Done!')


# def scATAC_GenerateBAMCoord_INPUT(outdirectory, coordinate_file, assignment_file, count_mat_file, cellnumberfile, BED_filename, OUTPUT_cells_barcode_file, read_len=50, jitter_size=5):
# 	"""Generate Synthetic reads in BED format. 

# 	Parameters
# 	----------
# 	outdirectory: `str`
# 		Output directory of the synthteic bed file and its corresponding cell barcodes file.
# 	coordinate_file: `str`
# 		Directory of the coordinates file storing synthetic reads output by scATAC_SampleSyntheticReads.
# 	assignment_file: `str`
# 		Directory of the peak assignment file generated in match_peak.
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
# 	jitter_size: `int` (default: '50')
# 		Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
# 	"""
# 	random.seed(2022)
# 	read_lines = pd.read_csv(coordinate_file, delimiter="\t",  names=['true_peak_name', 'chr', 'r1_start', 'r2_start', 'length'])
# 	peaks_assignments = pd.read_csv(assignment_file, delimiter="\t",  names=['true_chr', 'true_start', 'true_end', 'ref_chr', 'ref_start', 'ref_end']).to_numpy()
# 	count_mat = pd.read_csv(count_mat_file, header=0, delimiter="\t").to_numpy()
# 	marginal_cell_number = pd.read_csv(cellnumberfile, header=None, delimiter="\t").to_numpy()
# 	n_cell = np.shape(count_mat)[1]
# 	random_cellbarcode_list = cellbarode_generator(n_cell, size=16)
# 	read1_bedfile="%s.read1.bed" % BED_filename
# 	read2_bedfile="%s.read2.bed" % BED_filename
# 	start = time.time()
# 	with open(outdirectory + "/" + OUTPUT_cells_barcode_file, 'w') as f:
# 		for item in random_cellbarcode_list:
# 		    f.write("%s\n" % item)
# 	with open("%s/%s" % (outdirectory, read1_bedfile), 'w') as fp:
# 		pass
# 	with open("%s/%s" % (outdirectory, read2_bedfile), 'w') as fp:
# 		pass
# 	peak_nonzero_id = np.nonzero(marginal_cell_number)[0]
# 	for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
# 		peak_ind = peak_nonzero_id[relative_peak_ind]
# 		# peak_ind = 192
# 		peak_record = peaks_assignments[peak_ind]
# 		count_vec = count_mat[peak_ind,:] # Input
# 		print(peak_ind)
# 		read_1_df, read_2_df = scATAC_INPUT_PerTruePeakEdition(peak_record, count_vec, read_lines, random_cellbarcode_list, read_len, jitter_size)
# 		read_1_df.to_csv("%s/%s" % (outdirectory, read1_bedfile), header=None, index=None, sep='\t', mode='a')
# 		read_2_df.to_csv("%s/%s" % (outdirectory, read2_bedfile), header=None, index=None, sep='\t', mode='a')
# 	end = time.time()
# 	print(end - start)

def scATAC_GenerateBAMCoord_INPUT(count_mat_filename, samtools_directory, INPUT_bamfile, assignment_file, directory_cellnumber, outdirectory, BED_filename, OUTPUT_cells_barcode_file, read_len = 50, jitter_size = 5):
	"""Generate Synthetic reads in BED format. 

	Parameters
	----------
	count_mat_filename: `str`
		The base name of output count matrix in bam2countmat.
	samtools_directory: `str`
		Path to software samtools.
	INPUT_bamfile: `str`
		Input BAM file for anlaysis.
	assignment_file: `str`
		Features assignment file output by function `match_peak`.
	directory_cellnumber: `str`
		Directory of the marginal synthetic count vector file output in scATAC_GenerateSyntheticCount step.
	outdirectory: `str`
		Specify the output directory for synthetic reads bed file.
	BED_filename: `str`
		Specify the base name of output bed file.
	OUTPUT_cells_barcode_file: `str`
		Specify the file name storing the synthetic cell barcodes.
	read_len: `int` (default: '50')
		Specify the length of synthetic reads. Default value is 50 bp.
	jitter_size: `int` (default: '5')
		Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
	"""
	print('scReadSim scRNA_GenerateBAMCoord Running...')
	print('\t- Sampling synthetic reads...')
	coordinate_file = count_mat_filename + ".BAMfile_halfsampled_coordinates.txt"
	cellnumberfile = "%s/%s.scDesign2Simulated.nPairsRegionmargional.txt" % (directory_cellnumber, count_mat_filename)
	rm_coor_command = "rm %s/%s" % (outdirectory, coordinate_file)
	os.system(rm_coor_command)
	create_coor_command = "touch %s/%s" % (outdirectory, coordinate_file)
	os.system(create_coor_command)
	cmd = "while true; do read -r region <&3 || break;  read -r ncell <&4 || break; true_region=$(echo ${region} | cut -f 1,2,3 | perl -lane 'print \"$F[0]:$F[1]-$F[2]\"'); ref_region=$(echo ${region} | cut -f 4,5,6 | perl -lane 'print \"$F[0]:$F[1]-$F[2]\"'); paste -d\"\t\" <(awk -v nsample=${ncell} -v region=${true_region} 'BEGIN{for(c=0;c<nsample;c++) print region}') <(%s/samtools view %s ${ref_region} | shuf -r -n ${ncell} | cut -f3,4,8,9) >> %s/%s;  done 3<%s 4<%s" % (samtools_directory, INPUT_bamfile, outdirectory, coordinate_file, assignment_file, cellnumberfile)
	# Testing using copied directory
	output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to generate synthetic reads:\n', error.decode())
	print('Done!')
	print('\t- Sampling synthetic reads done!')
	print('\t- Writing synthetic reads...')
	count_mat_file = "%s.scDesign2Simulated.txt" % count_mat_filename
	random.seed(2022)
	read_lines = pd.read_csv(outdirectory + "/" + coordinate_file, delimiter="\t",  names=['true_peak_name', 'chr', 'r1_start', 'r2_start', 'length'])
	peaks_assignments = pd.read_csv(assignment_file, delimiter="\t",  names=['true_chr', 'true_start', 'true_end', 'ref_chr', 'ref_start', 'ref_end']).to_numpy()
	count_mat = pd.read_csv(outdirectory + "/" + count_mat_file, header=0, delimiter="\t").to_numpy()
	marginal_cell_number = pd.read_csv(cellnumberfile, header=None, delimiter="\t").to_numpy()
	n_cell = np.shape(count_mat)[1]
	random_cellbarcode_list = cellbarode_generator(n_cell, size=16)
	read1_bedfile="%s.read1.bed" % BED_filename
	read2_bedfile="%s.read2.bed" % BED_filename
	with open(outdirectory + "/" + OUTPUT_cells_barcode_file, 'w') as f:
		for item in random_cellbarcode_list:
		    f.write("%s\n" % item)
	with open("%s/%s" % (outdirectory, read1_bedfile), 'w') as fp:
		pass
	with open("%s/%s" % (outdirectory, read2_bedfile), 'w') as fp:
		pass
	peak_nonzero_id = np.nonzero(marginal_cell_number)[0]
	for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
		peak_ind = peak_nonzero_id[relative_peak_ind]
		# peak_ind = 192
		peak_record = peaks_assignments[peak_ind]
		count_vec = count_mat[peak_ind,:] # Input
		print(peak_ind)
		read_1_df, read_2_df = scATAC_INPUT_PerTruePeakEdition(peak_record, count_vec, read_lines, random_cellbarcode_list, read_len, jitter_size)
		read_1_df.to_csv("%s/%s" % (outdirectory, read1_bedfile), header=None, index=None, sep='\t', mode='a')
		read_2_df.to_csv("%s/%s" % (outdirectory, read2_bedfile), header=None, index=None, sep='\t', mode='a')
	print('\t- Writing synthetic reads done.')
	print('Done!\n')




def scATAC_CombineBED(outdirectory, BED_filename_pre, BED_COMPLE_filename_pre, BED_filename_combined_pre):
	"""Combine the bed files of foreground and background feature sets into one bed file.

	"""
	combine_read1_cmd = "cat %s/%s.read1.bed %s/%s.read1.bed | sort -k1,1 -k2,2n > %s/%s.read1.bed" % (outdirectory, BED_filename_pre, outdirectory, BED_COMPLE_filename_pre, outdirectory, BED_filename_combined_pre)
	output, error = subprocess.Popen(combine_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to create combine synthetic read1 bed files:\n', error.decode())
	combine_read2_cmd = "cat %s/%s.read2.bed %s/%s.read2.bed | sort -k1,1 -k2,2n > %s/%s.read2.bed" % (outdirectory, BED_filename_pre, outdirectory, BED_COMPLE_filename_pre, outdirectory, BED_filename_combined_pre)
	output, error = subprocess.Popen(combine_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to create combine synthetic read2 bed files:\n', error.decode())
		# sys.exit('[ERROR] Fail to create combine synthetic read2 bed files:\n', error.decode())


def scATAC_BED2FASTQ(bedtools_directory, seqtk_directory, referenceGenome_file, outdirectory, BED_filename_combined, sort_FASTQ = True):
	"""Convert Synthetic reads from BED to FASTQ. 

	Parameters
	----------
	bedtools_directory: `str`
		Directory of software bedtools.
	seqtk_directory: `str`
		Directory of software seqtk.
	referenceGenome_file: `str`
		Directory of the reference genome FASTA file that the synthteic reads should align.
	outdirectory: `str`
		Output directory of the synthteic bed file and its corresponding cell barcodes file.
	BED_filename_combined: `str`
		Specify the base name of output bed file.
	sort_FASTQ: `bool` (Default: True)
		Set `True` to sort the output FASTQ file.
	"""
	# Create FASTA
	print('scReadSim BED2FASTQ_Pair Running...')
	print('\t- Creating FASTA files...')
	fasta_read1_cmd = "%s/bedtools getfasta -s -fi %s -bed %s/%s.read1.bed -fo %s/%s.read1.bed2fa.strand.fa -nameOnly" % (bedtools_directory, referenceGenome_file, outdirectory, BED_filename_combined, outdirectory, BED_filename_combined)
	output, error = subprocess.Popen(fasta_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print(error.decode())
	fasta_read2_cmd = "%s/bedtools getfasta -s -fi %s -bed %s/%s.read2.bed -fo %s/%s.read2.bed2fa.strand.fa -nameOnly" % (bedtools_directory, referenceGenome_file, outdirectory, BED_filename_combined, outdirectory, BED_filename_combined)
	output, error = subprocess.Popen(fasta_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print(error.decode())
	print('\t- Converting FASTA files to FASTQ files...')
 	# remove (-) or (+)
	org_fasta_read1_cmd = "sed '/^>/s/.\{3\}$//' %s/%s.read1.bed2fa.strand.fa > %s/%s.read1.bed2fa.fa" % (outdirectory, BED_filename_combined, outdirectory, BED_filename_combined)
	output, error = subprocess.Popen(org_fasta_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to remove strand infomormation from synthetic read1 fasta file:', error.decode())
	org_fasta_read2_cmd = "sed '/^>/s/.\{3\}$//' %s/%s.read2.bed2fa.strand.fa > %s/%s.read2.bed2fa.fa" % (outdirectory, BED_filename_combined, outdirectory, BED_filename_combined)
	output, error = subprocess.Popen(org_fasta_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to remove strand infomormation from synthetic read2 fasta file:', error.decode())
	# FASTA to FASTQ
	fastq_read1_cmd = "%s/seqtk seq -F 'F' %s/%s.read1.bed2fa.fa > %s/%s.read1.bed2fa.fq" % (seqtk_directory, outdirectory, BED_filename_combined, outdirectory, BED_filename_combined)
	output, error = subprocess.Popen(fastq_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to convert read1 synthetic fasta file to fastq file:', error.decode())
	fastq_read2_cmd = "%s/seqtk seq -F 'F' %s/%s.read2.bed2fa.fa > %s/%s.read2.bed2fa.fq" % (seqtk_directory, outdirectory, BED_filename_combined, outdirectory, BED_filename_combined)
	output, error = subprocess.Popen(fastq_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	     print('[ERROR] Fail to convert read2 synthetic fasta file to fastq file:', error.decode())
	print('\t- FASTQ files %s.read1.bed2fa.fq, %s.read2.bed2fa.fq stored in %s.' % (BED_filename_combined, BED_filename_combined, outdirectory))
	if sort_FASTQ == True:
		print('\t- Sorting FASTQ files...')
		sort_fastq_read1_cmd = "cat %s/%s.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.read1.bed2fa.sorted.fq" % (outdirectory, BED_filename_combined, outdirectory, BED_filename_combined)
		output, error = subprocess.Popen(sort_fastq_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if error:
		     print('[ERROR] Fail to sort read1 synthetic fastq file:', error.decode())
		sort_fastq_read2_cmd = "cat %s/%s.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.read2.bed2fa.sorted.fq" % (outdirectory, BED_filename_combined, outdirectory, BED_filename_combined)
		output, error = subprocess.Popen(sort_fastq_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if error:
		     print('[ERROR] Fail to sort read2 synthetic fastq file:', error.decode())
		print('\t- Sorted FASTQ files %s.read1.bed2fa.sorted.fq, %s.read2.bed2fa.sorted.fq stored in %s.' % (BED_filename_combined, BED_filename_combined, outdirectory))
	print('Done!\n')


def AlignSyntheticBam_Pair(bowtie2_directory, samtools_directory, outdirectory, referenceGenome_name, referenceGenome_dir, BED_filename_combined, output_BAM_pre):
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
	BED_filename_combined: `str`
		Base name of the combined bed file output by function `scRNA_CombineBED`.
	output_BAM_pre: `str`
		Specify the base name of the output BAM file.
	"""
	print('scReadSim AlignSyntheticBam_Pair Running...')
	# if doIndex == True:
	# 	print('\t- Indexing reference genome file...')
	# 	index_ref_cmd = "%s/bowtie2-build %s/%s.fa %s" % (bowtie2_directory, referenceGenome_dir, referenceGenome_name, referenceGenome_name)
	# 	output, error = subprocess.Popen(index_ref_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	# 	if error:
	# 	     print('[ERROR] Fail to index the reference genome:\nPlease index the reference genome by setting \'doIndex=True\'\n', error.decode())
	# Align using bwa
	print('\t- Aligning FASTQ files onto reference genome files...')
	alignment_cmd = "%s/bowtie2 -x %s/%s -1 %s/%s.read1.bed2fa.sorted.fq -2 %s/%s.read2.bed2fa.sorted.fq | %s/samtools view -bS - > %s/%s.synthetic.noCB.bam" % (bowtie2_directory, referenceGenome_dir, referenceGenome_name,  outdirectory, BED_filename_combined, outdirectory, BED_filename_combined, samtools_directory, outdirectory, output_BAM_pre)
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
	sortBAMcmd = "%s/samtools sort %s/%s.synthetic.bam > %s/%s.synthetic.sorted.bam" % (samtools_directory, outdirectory, output_BAM_pre, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(sortBAMcmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to sort synthetic BAM file:', error.decode())
	indexBAMcmd = "%s/samtools index %s/%s.synthetic.sorted.bam" % (samtools_directory, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(indexBAMcmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to index synthetic BAM file:', error.decode())
	print('Done!\n')














