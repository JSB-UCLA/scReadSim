import numpy as np
import rpy2.robjects as robjects
import time
import pkg_resources
import os

def scATAC_GenerateSyntheticCount(count_mat_filename, directory, outdirectory, n_cell_new=None, total_count_new=None, celllabel_file=None):
	"""Simulate synthetic count matrix.

	Parameters
	----------
	count_mat_filename: `str`
		Base name of the count matrix output by function bam2countmat().
	directory: `str`
		Path to the count matrix.
	outdirectory: `str`
		Output directory of coordinate files.
	n_cell_new: `int` (default: None)
		Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
	total_count_new: `int` (default: None)
		Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
	celllabel_file: `str` (default: None)
		Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.
	"""
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	scATAC_runSyntheticCount = robjects.globalenv['scATAC_runSyntheticCount']
	if n_cell_new == None:
		n_cell_new = "default"
	if total_count_new == None:
		total_count_new = "default"
	if celllabel_file == None:
		celllabel_file = "default"
	scATAC_runSyntheticCount(count_mat_filename, directory, outdirectory, n_cell_new, total_count_new, celllabel_file)
	# if cluster_prestep == True:
	# 	scATAC_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 1)
	# else:
	# 	scATAC_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 0)


def scRNA_GenerateSyntheticCount(count_mat_filename, directory, outdirectory, UMI_modeling=False, UMI_count_mat_filename="UMI_countmat", n_cell_new=None, total_count_new=None, celllabel_file=None):
	"""Simulate synthetic count matrix.

	Parameters
	----------
	count_mat_filename: `str`
		Base name of the count matrix output by function scRNA_bam2countmat().
	directory: `str`
		Path to the count matrix.
	outdirectory: `str`
		Output directory of coordinate files.
	UMI_modeling: `bool` (default: False)
		Specify whether scReadSim should model UMI count of the input BAM file.
	UMI_count_mat_filename: `str` (default: 'UMI_countmat')
		Base name of the UMI count matrix output by function `scRNA_bam2countmat` with option UMI_modeling setting to Ture.
	n_cell_new: `int` (default: None)
		Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
	total_count_new: `int` (default: None)
		Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
	celllabel_file: `str` (default: None)
		Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file. If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.
	"""
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	scRNA_runSyntheticCount = robjects.globalenv['scRNA_runSyntheticCount']
	if n_cell_new == None:
		n_cell_new = "default"
	if total_count_new == None:
		total_count_new = "default"
	if celllabel_file == None:
		celllabel_file = "default"
	scRNA_runSyntheticCount(count_mat_filename, directory, outdirectory, UMI_modeling, UMI_count_mat_filename, n_cell_new, total_count_new, celllabel_file)
	# if cluster_prestep == True:
	# 	scRNA_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 1)
	# else:
	# 	scRNA_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 0)
