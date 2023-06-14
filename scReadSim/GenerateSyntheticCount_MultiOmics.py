import numpy as np
import rpy2.robjects as robjects
import os


def scMultiOmics_GenerateSyntheticCount(RNA_count_mat_filename, ATAC_count_mat_filename, directory, outdirectory, n_cell_new=None, celllabel_file=None, n_cluster=None, n_cores=1):
	"""Simulate synthetic multiomic count matrices.

	Parameters
	----------
	RNA_count_mat_filename: `str`
		Base name of the count matrix output by function scRNA_bam2countmat_paral().
	ATAC_count_mat_filename: `str`
		Base name of the count matrix output by function scATAC_bam2countmat_paral().
	directory: `str`
		Path to the count matrix.
	outdirectory: `str`
		Output directory of coordinate files.
	n_cell_new: `int` (default: None)
		Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
	celllabel_file: `str` (default: None)
		Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file (and the columns of real count matrix). If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign3.
	n_cores: `int` (default: 1)
		Number of cores for parallel computing.
	"""
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions_MultiOmics.R')
	print(rscript_dir)
	r['source'](rscript_dir)
	scMultiOmics_runSyntheticCount = robjects.globalenv['scMultiOmics_runSyntheticCount']
	if n_cell_new == None:
		n_cell_new = "default"
	if celllabel_file == None:
		celllabel_file = "default"
	if n_cluster == None:
		n_cluster = "default"
	scMultiOmics_runSyntheticCount(RNA_count_mat_filename, ATAC_count_mat_filename, directory, outdirectory, n_cell_new, celllabel_file, n_cluster, n_cores)
	print("[scReadSim] Created:")
	print("[scReadSim] Synthetic scMultiOmics (ATAC modality) count matrix: %s.scMultiOmics.scDesign3Simulated.ATAC.txt" % RNA_count_mat_filename)
	print("[scReadSim] Synthetic scMultiOmics (RNA modality) count matrix: %s.scMultiOmics.scDesign3Simulated.RNA.txt" % ATAC_count_mat_filename)
	print("[scReadSim] Synthetic cell label file: scMultiOmics.scDesign3Simulated.CellTypeLabel.txt")
