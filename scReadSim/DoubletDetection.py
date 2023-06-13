import numpy as np
import rpy2.robjects as robjects
import os

def detectDoublet(count_mat_filename, directory, outdirectory, omic_choice):
	"""Simulate synthetic count matrix.

	Parameters
	----------
	count_mat_filename: `str`
		Base name of the count matrix output by function scRNA_bam2countmat().
	directory: `str`
		Path to the count matrix.
	outdirectory: `str`
		Output directory of coordinate files.
	omic_choice: `str`
		Specify the omic choice for doublet detection procedure: "ATAC" or "RNA".	
	"""
	r = robjects.r
	# TODO
    # Need to edit the Rscript path
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/DoubletDetection_Rscript.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	get_doublet_id = robjects.globalenv['get_doublet_id']
	get_doublet_id(count_mat_filename, directory, outdirectory, omic_choice)
	print("[scReadSim] Created:")
	print("[scReadSim] Doublet classification result: %s/doublet_classification.Rdata" % outdirectory)
