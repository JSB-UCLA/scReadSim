import numpy as np
import rpy2.robjects as robjects
import os

def ConditionEffect_GenerateSyntheticCount(count_mat_filename, cellcondition_file, directory, outdirectory, n_cell_new=None, celllabel_file=None, n_cluster=None, n_cores=1):
	"""Simulate synthetic count matrix.

	Parameters
	----------
	count_mat_filename: `str`
		Base name of the count matrix.
	cellcondition_file: `str`
		Specify the one-column text file containing the predefined cell conditions. Make sure that the order of cell conditions correspond to the cell barcode file (and the columns of real count matrix).
	directory: `str`
		Path to the count matrix.
	outdirectory: `str`
		Output directory of coordinate files.
    n_cell_new: `int` (default: None)
		Number of synthetic cells. If not specified, scReadSim uses the number of real cells.
	celllabel_file: `str` (default: None)
		Specify the one-column text file containing the predefined cell labels. Make sure that the order of cell labels correspond to the cell barcode file (and the columns of real count matrix). If no cell labels are specified, scReadSim performs a Louvain clustering before implementing scDesign2/3.
    n_cores: `int` (default: 1)
		Number of cores for parallel computing.
	"""
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions_ConditionEffect.R')
	print(rscript_dir)
	r['source'](rscript_dir)
	ConditionEffect_runSyntheticCount = robjects.globalenv['ConditionEffect_runSyntheticCount']
	if n_cell_new == None:
		n_cell_new = "default"
	if celllabel_file == None:
		celllabel_file = "default"
	if n_cluster == None:
		n_cluster = "default"
		
	ConditionEffect_runSyntheticCount(count_mat_filename, cellcondition_file, directory, outdirectory, n_cell_new, celllabel_file, n_cluster, n_cores)
	print("[scReadSim] Created:")
	print("[scReadSim] Synthetic count matrix: %s.ConditionEffect.scDesign3Simulated.txt" % count_mat_filename)
	print("[scReadSim] Synthetic cell label file: ConditionEffect.scDesign3Simulated.CellTypeLabel.txt")
	print("[scReadSim] Synthetic cell condition file: ConditionEffect.scDesign3Simulated.ConditionLabel.txt")
