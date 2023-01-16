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
from pathlib import Path
from joblib import Parallel, delayed
import pysam
from collections import defaultdict


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


def scRNA_GenerateBAMCoord(bed_file,  UMI_count_mat_file, read_bedfile_prename, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file, jitter_size=5, read_len=90):
	"""Generate Synthetic reads in BED format. 

    Parameters
    ----------
    bed_file: `str`
        Input the bed file of features.
    UMI_count_mat_file: `str`
        The path to synthetic UMI count matrix.
    read_bedfile_prename: `str`
        Specify the base name of output bed file.
    INPUT_bamfile: `str`
        Input BAM file for anlaysis.
    outdirectory: `str`
        Specify the output directory for synthetic reads bed file.
    OUTPUT_cells_barcode_file: `str`
        Specify the file name storing the synthetic cell barcodes.
    jitter_size: `int` (default: '5')
        Specify the range of random shift to avoid replicate synthetic reads. Default value is 5 bp.
    read_len: `int` (default: '90')
        Specify the length of synthetic reads. Default value is 50 bp.
    """
	UMI_count_mat_df = pd.read_csv("%s" % (UMI_count_mat_file), header=0, delimiter="\t")
	UMI_count_mat = UMI_count_mat_df.to_numpy()
	UMI_count_mat_cluster = UMI_count_mat_df.columns.to_numpy()
	n_cell = np.shape(UMI_count_mat)[1]
	samfile = pysam.AlignmentFile(INPUT_bamfile, "rb")
	with open(bed_file) as open_peak:
		reader = csv.reader(open_peak, delimiter="\t")
		open_peak = np.asarray(list(reader))
	peak_nonzero_id = np.nonzero(UMI_count_mat.sum(axis=1))[0]
	random.seed(2022)
	random_cellbarcode_list = cellbarcode_generator(n_cell, size=16)
	with open(OUTPUT_cells_barcode_file, 'w') as f:
			for item in random_cellbarcode_list:
				f.write(item + "\n")
	cellbarcode_list_withclusters = np.vstack((random_cellbarcode_list, UMI_count_mat_cluster)).transpose()
	with open(OUTPUT_cells_barcode_file + ".withSynthCluster", 'w') as f:
			for item in cellbarcode_list_withclusters:
				f.write("\t".join(item) + "\n")
	with open("%s/%s.read.bed" % (outdirectory, read_bedfile_prename), 'w') as fp:
		pass
	print("[scReadSim] Generating Synthetic Reads for Feature Set: %s" % (bed_file))
	for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
		peak_ind = peak_nonzero_id[relative_peak_ind]
		# peak_ind = 1
		rec = open_peak[peak_ind]
		rec_name = '_'.join(rec)
		UMI_read_dict_pergene = defaultdict(lambda: [])
		reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
		# TODO: create a dictionary recording UMI(key) to reads starting positions(values)
		for read in reads:
			if read.has_tag('UB:Z'):
				UMI = read.get_tag('UB:Z')
				start = read.reference_start
				if read.is_reverse==1:
					strand = "-" 
				else:
					strand = "+"
				# strand = read.is_reverse
				UMI_read_dict_pergene[UMI].append(str(start) + ":" + str(strand))
		## Real data: Obtain the empirical distribution of # reads of UMI in the gene 
		UMI_real_vec = list(UMI_read_dict_pergene.keys())
		nread_perUMI = [len(x) for x in UMI_read_dict_pergene.values()] # return a vector with the length of # UMIs in the gene
		# nread_perUMI_prob = nread_perUMI / np.sum(nread_perUMI)
		## Synthetic data
		# Sample syntheitc UMI's read number and assign real UMI to synthetic UMI (only for sampling reads)
		UMI_count_vec = UMI_count_mat[peak_ind,:] # Synthetic umi count
		nonzer_UMI_zero_id = np.nonzero(UMI_count_vec)[0]
		synthetic_nread_perUMI_percell = [np.random.choice(nread_perUMI, size=UMI_count_vec[i], replace=True) for i in nonzer_UMI_zero_id] # return a vector with the length of non-zero count cells, each entry idnicates the number of reads for each UMI
		nread_synthetic = sum(np.sum(x) for x in synthetic_nread_perUMI_percell) # total number of synthetic reads
		# Assign real UMI to synthetic UMI based on read number per UMI
		# Construct a dic with keys as number of reasd per UMI and values as the corresponding UMI strings
		nread_UMI_dic = defaultdict(lambda: [])
		for i,item in enumerate(nread_perUMI):
			nread_UMI_dic[item].append(UMI_real_vec[i])
		synthetic_nread_perUMI_percell_unlist = [xs for x in synthetic_nread_perUMI_percell for xs in list(x)]
		synthetic_realUMI_assign_array = np.array(synthetic_nread_perUMI_percell_unlist).astype(str)
		synthetic_nread_perUMI_percell_unlist_array = np.array(synthetic_nread_perUMI_percell_unlist)
		for i in set(synthetic_nread_perUMI_percell_unlist):
			tmp_ind = np.where(synthetic_nread_perUMI_percell_unlist_array == i)[0]
			synthetic_realUMI_assign_array[tmp_ind] = np.random.choice(nread_UMI_dic[i], size=len(tmp_ind), replace=True)
		synthetic_realUMI_assign_vec = list(synthetic_realUMI_assign_array) # a vectir sane length to synthetic_nread_perUMI_percell_unlist, synthetic UMI's number 
		# synthetic_realUMI_assign_vec = [np.random.choice(nread_UMI_dic[nread], size=1, replace=True) for nread in synthetic_nread_perUMI_percell_unlist]
		# 
		# Sample reads for each synthetic UMI
		synthetic_read_list = [np.random.choice(UMI_read_dict_pergene[synthetic_realUMI_assign_vec[id]], size=synthetic_nread_perUMI_percell_unlist[id], replace=True) for id in range(len(synthetic_realUMI_assign_vec))]
		synthetic_read_unlist = [xs for x in synthetic_read_list for xs in list(x)]
		synthetic_read_start_list = [int(xs.split(':', 1)[0]) - 1  for xs in synthetic_read_unlist]  # -1 accounts for the BAM coordinates is +1 compared with getfasta bedtools
		synthetic_read_strand_list = [str(xs.split(':', 1)[1]) for xs in synthetic_read_unlist]
		# Creat synthetic UMI string
		UMI_synthetic_uniq = cellbarcode_generator(np.sum(UMI_count_vec), size=10)
		UMI_synthetic_list = np.repeat(np.array(UMI_synthetic_uniq), synthetic_nread_perUMI_percell_unlist)
		# Create syntheitc CB string
		CB_synthetic_repeat = [np.sum(x) for x in synthetic_nread_perUMI_percell] # same length to nozero cells
		CB_synthetic_list = np.repeat(np.array(random_cellbarcode_list)[nonzer_UMI_zero_id], CB_synthetic_repeat)
		# Create syntheitc read name string
		read_name_remaining_list = ["CellNo" + str(nonzer_UMI_zero_id[ind] + 1) + ":" + str(rec_name) + "#" + str(count).zfill(4) for ind in range(len(CB_synthetic_repeat)) for count in range(CB_synthetic_repeat[ind])]
		read_name_list = [CB_synthetic_list[id] + UMI_synthetic_list[id] + ":" + read_name_remaining_list[id] for id in range(nread_synthetic)]
		# Create output table
		jitter_value_vec = np.random.randint(-jitter_size,jitter_size,size=nread_synthetic)  # nrow(reads_cur) should equal to nfrag_cur
		reads_cur = pd.DataFrame({
			'chr': rec[0],
			'r1_start_shifted': synthetic_read_start_list  + jitter_value_vec,
			'r1_end_shifted': synthetic_read_start_list + jitter_value_vec + read_len,
			'read_name': read_name_list,
			'read_length': read_len,
			'strand': synthetic_read_strand_list
			})
		reads_cur.to_csv("%s/%s.read.bed" % (outdirectory, read_bedfile_prename), header=None, index=None, sep='\t', mode='a')
	print("\n[scReadSim] Created:")
	print("[scReadSim] Read bed file: %s/%s.read.bed" % (outdirectory, read_bedfile_prename))


def scRNA_CombineBED(outdirectory, gene_read_bedfile_prename, intergene_read_bedfile_prename, BED_filename_combined_pre):
	"""Combine the bed files of foreground and background feature sets into one bed file.

	Parameters
	----------
	outdirectory: `str`
		Directory of `gene_read_bedfile_prename`.txt and `intergene_read_bedfile_prename`.txt.
	gene_read_bedfile_prename: 'str'
		File prename of foreground synthetic reads bed file.
	intergene_read_bedfile_prename: 'str'
		File prename of background synthetic reads bed file.
	BED_filename_combined_pre: 'str'
		Specify the combined syntehtic reads bed file prename. The combined bed file will be output to `outdirectory`.
	"""
	print("[scReadSim] Combining Synthetic Read Bed Files from Genes and InterGenes.")
	combine_read_cmd = "cat %s/%s.read.bed %s/%s.read.bed | sort -k1,1 -k2,2n > %s/%s.read.bed" % (outdirectory, gene_read_bedfile_prename, outdirectory, intergene_read_bedfile_prename, outdirectory, BED_filename_combined_pre)
	output, error = subprocess.Popen(combine_read_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to create combine synthetic read bed files:\n', error.decode())
	print("\n[scReadSim] Created:")
	print("[scReadSim] Combined Read Bed File: %s/%s.read.bed" % (outdirectory, BED_filename_combined_pre))
	print("[scReadSim] Done.")


def scRNA_BED2FASTQ(bedtools_directory, seqtk_directory, referenceGenome_file, outdirectory, BED_filename_combined, synthetic_fastq_prename):
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
	"""
	# Create FASTA
	print('[scReadSim] Generating Synthetic Read FASTA files...')
	fasta_read2_cmd = "%s/bedtools getfasta -s -nameOnly -fi %s -bed %s/%s.read.bed -fo %s/%s.read2.strand.bed2fa.fa" % (bedtools_directory, referenceGenome_file, outdirectory, BED_filename_combined, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(fasta_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
			print(error.decode())
	# remove (-) or (+)
	org_fasta_read2_cmd = "sed '/^>/s/.\{3\}$//' %s/%s.read2.strand.bed2fa.fa > %s/%s.read2.bed2fa.fa" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(org_fasta_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to remove strand infomormation from synthetic read fasta file:\n', error.decode())
	# Create Read 1 from Read 2
	fasta_read1_cmd = "awk 'NR%%2==0 {print substr(p,2,26);} NR%%2 {p=$0;print p;}' %s/%s.read2.bed2fa.fa > %s/%s.read1.bed2fa.fa" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(fasta_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
			print(error.decode())
	# FASTA to FASTQ
	print('[scReadSim] Generating Synthetic Read FASTQ files...')
	fastq_read1_cmd = "%s/seqtk seq -F 'F' %s/%s.read1.bed2fa.fa > %s/%s.read1.bed2fa.fq" % (seqtk_directory, outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(fastq_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
			print('[ERROR] Fail to convert read1 synthetic fasta file to fastq file:', error.decode())
	fastq_read2_cmd = "%s/seqtk seq -F 'F' %s/%s.read2.bed2fa.fa > %s/%s.read2.bed2fa.fq" % (seqtk_directory, outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(fastq_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
			print('[ERROR] Fail to convert read2 synthetic fasta file to fastq file:', error.decode())
	print('[scReadSim] Sorting FASTQ files...')
	sort_fastq_read1_cmd = "cat %s/%s.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.read1.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(sort_fastq_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
			print('[ERROR] Fail to sort read1 synthetic fastq file:', error.decode())
	sort_fastq_read2_cmd = "cat %s/%s.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.read2.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(sort_fastq_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
			print('[ERROR] Fail to sort read2 synthetic fastq file:', error.decode())
	print("\n[scReadSim] Created:")
	print("[scReadSim] Read 1 FASTQ File: %s/%s.read1.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename))
	print("[scReadSim] Read 2 FASTQ File: %s/%s.read2.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename))
	print("[scReadSim] Done.")


def AlignSyntheticBam_Single(bowtie2_directory, samtools_directory, outdirectory, referenceGenome_name, referenceGenome_dir, synthetic_fastq_prename, output_BAM_pre):
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
		Base name of the synthetic FASTQ file output by function `scRNA_BED2FASTQ`.
	output_BAM_pre: `str`
		Specify the base name of the output BAM file.
	"""
	print('[scReadSim] Aligning FASTQ files onto Reference Genome Files with Bowtie2...')
	alignment_cmd = "%s/bowtie2 -x %s/%s -U %s/%s.read2.bed2fa.sorted.fq | %s/samtools view -bS - > %s/%s.synthetic.noCB.bam" % (bowtie2_directory, referenceGenome_dir, referenceGenome_name,  outdirectory, synthetic_fastq_prename, samtools_directory, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(alignment_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	print('[Bowtie2] Aligning:\n', error.decode())
	print('[scReadSim] Alignment Done.')
	print('[scReadSim] Generating Cell Barcode and UMI Barcode Tags...')
	addBC2BAM_header_cmd = "%s/samtools view %s/%s.synthetic.noCB.bam -H > %s/%s.synthetic.noCB.header.sam" % (samtools_directory, outdirectory, output_BAM_pre, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(addBC2BAM_header_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	addBC2BAM_cmd = "cat <( cat %s/%s.synthetic.noCB.header.sam ) <( paste <(%s/samtools view %s/%s.synthetic.noCB.bam ) <(%s/samtools view %s/%s.synthetic.noCB.bam | cut -f1 | cut -d':' -f1 | awk '{s=substr($1,1,16)}{g=substr($1,17,length($1))}{printf \"CB:Z:%%s\tUB:Z:%%s\\n\",s,g;}')) | %s/samtools view -bS - > %s/%s.synthetic.bam" % (outdirectory, output_BAM_pre, samtools_directory, outdirectory, output_BAM_pre, samtools_directory, outdirectory, output_BAM_pre, samtools_directory, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(addBC2BAM_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to add CB and UB tags to synthetic BAM file:\n', error.decode())
	print('[scReadSim] Sorting and Indexing BAM file...')
	sortBAMcmd = "%s/samtools sort %s/%s.synthetic.bam > %s/%s.synthetic.sorted.bam" % (samtools_directory, outdirectory, output_BAM_pre, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(sortBAMcmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to sort synthetic BAM file:\n', error.decode())
	indexBAMcmd = "%s/samtools index %s/%s.synthetic.sorted.bam" % (samtools_directory, outdirectory, output_BAM_pre)
	output, error = subprocess.Popen(indexBAMcmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[ERROR] Fail to index synthetic BAM file:\n', error.decode())
	print("\n[scReadSim] Created:")
	print("[scReadSim] Synthetic Read BAM File: %s/%s.synthetic.sorted.bam" % (outdirectory, output_BAM_pre))
	print("[scReadSim] Done.")


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
	print('[scReadSim] Substitution Error Calculating...')
	combine_read1_cmd = "java -jar %s ErrorRateByReadPosition -i %s -r %s -o %s/Real --collapse false" % (fgbio_jarfile, INPUT_bamfile, referenceGenome_file, outdirectory)
	output, error = subprocess.Popen(combine_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
		print('[Messages] Running fgbio on real bam file:\n', error.decode())
	# Generate Errors into fastq files
	real_error_rate_file = outdirectory + "/" + "Real.error_rate_by_read_position.txt"
	SubstiError(real_error_rate_file, outdirectory, synthetic_fastq_prename)
	# Combine FASTQs
	print('[scReadSim] Sorting FASTQ files...')
	sort_fastq_read1_cmd = "cat %s/%s.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.ErrorIncluded.read1.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(sort_fastq_read1_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
			print('[ERROR] Fail to sort read1 synthetic fastq file:', error.decode())
	sort_fastq_read2_cmd = "cat %s/%s.ErrorIncluded.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > %s/%s.ErrorIncluded.read2.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename, outdirectory, synthetic_fastq_prename)
	output, error = subprocess.Popen(sort_fastq_read2_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
			print('[ERROR] Fail to sort read2 synthetic fastq file:', error.decode())
	print("\n[scReadSim] Created:")
	print("[scReadSim] Read 1 FASTQ File with Substitution Error: %s/%s.ErrorIncluded.read1.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename))
	print("[scReadSim] Read 2 FASTQ File with Substitution Error: %s/%s.ErrorIncluded.read2.bed2fa.sorted.fq" % (outdirectory, synthetic_fastq_prename))
	print("[scReadSim] Done.")



