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



def scATAC_GenerateSyntheticCount(count_mat_filename, directory, outdirectory, doub_classification_label_file=None, n_cell_new=None, total_count_new=None, celllabel_file=None, n_cluster=None, n_cores=1):
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
	n_cores: `int` (default: 1)
		Number of cores for parallel computing.
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
	scATAC_runSyntheticCount(count_mat_filename, directory, outdirectory, doub_classification_label_file, n_cell_new, total_count_new, celllabel_file, n_cluster, n_cores)
	print("\n[scReadSim] Created:")
	print("[scReadSim] Synthetic count matrix: %s.scDesign2Simulated.txt" % count_mat_filename)
	print("[scReadSim] Synthetic cell label file: %s.scDesign2Simulated.CellTypeLabel.txt" % count_mat_filename)



def scRNA_GenerateSyntheticCount(count_mat_filename, directory, outdirectory, doub_classification_label_file=None, n_cell_new=None, total_count_new=None, celllabel_file=None, n_cluster=None, n_cores=1):
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
	n_cores: `int` (default: 1)
		Number of cores for parallel computing.
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
	scRNA_runSyntheticCount(count_mat_filename, directory, outdirectory, doub_classification_label_file, n_cell_new, total_count_new, celllabel_file, n_cluster, n_cores)
	print("\n[scReadSim] Created:")
	print("[scReadSim] Synthetic count matrix: %s.scDesign2Simulated.txt" % count_mat_filename)
	print("[scReadSim] Synthetic cell label file: %s.scDesign2Simulated.CellTypeLabel.txt" % count_mat_filename)




def scATAC_GenerateSyntheticCount_MultiSample(INPUT_bamfile, outdirectory, n_cell_new=None, total_count_new=None, n_cluster=None, n_cores=1):
	"""Multi-sample/replicate implement of scReadSim for simulating scATAC-seq synthetic count matrix.

	Parameters
	----------
	INPUT_bamfile: `str`
		List of input BAM files (use absolute paths to the BAM files).
	outdirectory: `str`
		Specify the working directory of scReadSim for generating intermediate and final output files.
	n_cell_new: `int` (default: None)
		Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
	total_count_new: `int` (default: None)
		Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
	n_cores: `int` (default: 1)
		Number of cores for parallel computing.
	"""
	# Generate synthetic count matrix
	for rep_id in range(len(INPUT_bamfile)):
		print("\n[scReadSim] Simulating synthetic count matrices for sample %s..." % str(rep_id+1))
		sample_output_d = outdirectory + "/" + "Rep" + str(rep_id+1)
		# Specify the output count matrices' prenames
		count_mat_peak_filename = "Rep%s.peak.countmatrix" % str(rep_id+1)
		count_mat_nonpeak_filename = "Rep%s.nonpeak.countmatrix" % str(rep_id+1)
		# Generate synthetic count matrix for peak-by-cell count matrix
		print("\n[scReadSim] Generating synthetic counts for peaks...")
		scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_peak_filename, directory=sample_output_d, outdirectory=sample_output_d, n_cluster=n_cluster, n_cores=n_cores)
		# Generate synthetic count matrix for nonpeak-by-cell count matrix
		print("\n[scReadSim] Generating synthetic counts for non-peaks...")
		celllabel_file = sample_output_d + "/" + count_mat_peak_filename + ".LouvainClusterResults.txt"
		scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_nonpeak_filename, directory=sample_output_d, outdirectory=sample_output_d, celllabel_file=celllabel_file, n_cell_new=n_cell_new, total_count_new=total_count_new, n_cores=n_cores)


def scRNA_GenerateSyntheticCount_MultiSample(INPUT_bamfile, outdirectory, n_cell_new=None, total_count_new=None, n_cluster=None, n_cores=1):
	"""Multi-sample/replicate implement of scReadSim for simulating scRNA-seq synthetic count matrix.

	Parameters
	----------
	INPUT_bamfile: `str`
		List of input BAM files (use absolute paths to the BAM files).
	outdirectory: `str`
		Specify the working directory of scReadSim for generating intermediate and final output files.
	n_cell_new: `int` (default: None)
		Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
	total_count_new: `int` (default: None)
		Number of (expected) sequencing depth. If not specified, scReadSim uses the real sequencing depth.
	n_cores: `int` (default: 1)
		Number of cores for parallel computing.
	"""
	# Generate synthetic count matrix
	for rep_id in range(len(INPUT_bamfile)):
		print("\n[scReadSim] Simulating synthetic count matrices for sample %s..." % str(rep_id+1))
		sample_output_d = outdirectory + "/" + "Rep" + str(rep_id+1)
		# Specify the output count matrices' prenames
		UMI_gene_count_mat_filename = "Rep%s.gene.countmatrix" % str(rep_id+1)
		UMI_intergene_count_mat_filename = "Rep%s.intergene.countmatrix" % str(rep_id+1)
		# Generate synthetic count matrix for peak-by-cell count matrix
		print("\n[scReadSim] Generating synthetic counts for genes...")
		scRNA_GenerateSyntheticCount(count_mat_filename=UMI_gene_count_mat_filename, directory=sample_output_d, outdirectory=sample_output_d, n_cluster=n_cluster, n_cores=n_cores)
		# Generate synthetic count matrix for nonpeak-by-cell count matrix
		print("\n[scReadSim] Generating synthetic counts for inter-genes...")
		celllabel_file = sample_output_d + "/" + UMI_gene_count_mat_filename + ".LouvainClusterResults.txt"
		scRNA_GenerateSyntheticCount(count_mat_filename=UMI_intergene_count_mat_filename, directory=sample_output_d, outdirectory=sample_output_d, celllabel_file=celllabel_file, n_cell_new=n_cell_new, total_count_new=total_count_new, n_cores=n_cores)
