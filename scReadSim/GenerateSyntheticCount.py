import numpy as np
import rpy2.robjects as robjects
import time
import pkg_resources
import os

def scATAC_GenerateSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep=True):
	"""Simulate synthetic count matrix.

	Parameters
	----------
	count_mat_filename: `str`
		Base name of the count matrix output by function bam2countmat().
	directory: `str`
		Path of the count matrix.
	outdirectory: `str`
		Output directory of coordinate files.
	cluster_prestep: `bool`
		Set `cluster_prestep=True` to perform a Louvain clustering before implementing scDesign2.
	"""
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	scATAC_runSyntheticCount = robjects.globalenv['scATAC_runSyntheticCount']
	if cluster_prestep == True:
		scATAC_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 1)
	else:
		scATAC_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 0)

def scRNA_GenerateSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep=True):
	"""Simulate synthetic count matrix.

	Parameters
	----------
	count_mat_filename: `str`
		Base name of the count matrix output by function bam2countmat().
	directory: `str`
		Path of the count matrix.
	outdirectory: `str`
		Output directory of coordinate files.
	cluster_prestep: `bool`
		Set `cluster_prestep=True` to perform a Louvain clustering before implementing scDesign2.
	"""
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	scRNA_runSyntheticCount = robjects.globalenv['scRNA_runSyntheticCount']
	if cluster_prestep == True:
		scRNA_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 1)
	else:
		scRNA_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 0)
