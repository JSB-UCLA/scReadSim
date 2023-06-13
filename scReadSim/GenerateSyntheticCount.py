import numpy as np
import rpy2.robjects as robjects
import time
import pkg_resources
import os
from rpy2.robjects.vectors import StrVector
import rpy2.robjects.packages as rpackages

# import R's utility package
utils = rpackages.importr('utils')

# select a mirror for R packages
utils.chooseCRANmirror(ind=1) # select the first mirror in the list

# R package names
# Packages Rsubread, ROGUE no hosted on CRAN are installed in Rscript
packnames = ('pscl', 'tidyverse', 'Seurat')

# Selectively install what needs to be install.
# We are fancy, just because we can.
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))


def scATAC_GenerateSyntheticCount(count_mat_filename, directory, outdirectory, doub_classification_label_file=None, n_cell_new=None, total_count_new=None, celllabel_file=None, n_cluster=None):
	"""Simulate synthetic count matrix.

	Parameters
	----------
	count_mat_filename: `str`
		Base name of the count matrix output by function bam2countmat().
	directory: `str`
		Path to the count matrix.
	outdirectory: `str`
		Output directory of coordinate files.
	doub_classification_label_file: `str`
		Specify the absolute path to the doublet classification result `doublet_classification.Rdata`.
	n_cell_new: `int` (default: None)
		Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
	total_count_new: `int` (default: None)
		Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
	celllabel_file: `str` (default: None)
		Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file (and the columns of real count matrix). If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.
	"""
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	scATAC_runSyntheticCount = robjects.globalenv['scATAC_runSyntheticCount']
	if doub_classification_label_file == None:
		doub_classification_label_file = "default"
	if n_cell_new == None:
		n_cell_new = "default"
	if total_count_new == None:
		total_count_new = "default"
	if celllabel_file == None:
		celllabel_file = "default"
	if n_cluster == None:
		n_cluster = "default"
	scATAC_runSyntheticCount(count_mat_filename, directory, outdirectory, doub_classification_label_file, n_cell_new, total_count_new, celllabel_file, n_cluster)
	print("[scReadSim] Created:")
	print("[scReadSim] Synthetic count matrix: %s.scDesign2Simulated.txt" % count_mat_filename)
	print("[scReadSim] Synthetic cell label file: %s.scDesign2Simulated.CellTypeLabel.txt" % count_mat_filename)
	# if cluster_prestep == True:
	# 	scATAC_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 1)
	# else:
	# 	scATAC_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 0)


def scRNA_GenerateSyntheticCount(count_mat_filename, directory, outdirectory, doub_classification_label_file=None, n_cell_new=None, total_count_new=None, celllabel_file=None, n_cluster=None):
	"""Simulate synthetic count matrix.

	Parameters
	----------
	count_mat_filename: `str`
		Base name of the count matrix output by function scRNA_bam2countmat().
	directory: `str`
		Path to the count matrix.
	outdirectory: `str`
		Output directory of coordinate files.
	doub_classification_label_file: `str`
		Specify the absolute path to the doublet classification result `doublet_classification.Rdata`.
	n_cell_new: `int` (default: None)
		Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
	total_count_new: `int` (default: None)
		Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
	celllabel_file: `str` (default: None)
		Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file (and the columns of real count matrix). If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2.
	"""
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	scRNA_runSyntheticCount = robjects.globalenv['scRNA_runSyntheticCount']
	if doub_classification_label_file == None:
		doub_classification_label_file = "default"
	if n_cell_new == None:
		n_cell_new = "default"
	if total_count_new == None:
		total_count_new = "default"
	if celllabel_file == None:
		celllabel_file = "default"
	if n_cluster == None:
		n_cluster = "default"
	scRNA_runSyntheticCount(count_mat_filename, directory, outdirectory, doub_classification_label_file, n_cell_new, total_count_new, celllabel_file, n_cluster)
	print("[scReadSim] Created:")
	print("[scReadSim] Synthetic count matrix: %s.scDesign2Simulated.txt" % count_mat_filename)
	print("[scReadSim] Synthetic cell label file: %s.scDesign2Simulated.CellTypeLabel.txt" % count_mat_filename)

	# if cluster_prestep == True:
	# 	scRNA_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 1)
	# else:
	# 	scRNA_runSyntheticCount(count_mat_filename, directory, outdirectory, cluster_prestep = 0)
