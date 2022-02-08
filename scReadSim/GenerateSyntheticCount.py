## This script takes use of Rscripts to generate synthetic count matrix.
import numpy as np
import rpy2.robjects as robjects
# import rpy2.robjects.numpy2ri as np2r
# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri
import time
import pkg_resources
import os

def scATAC_GenerateSyntheticCount(samplename, sample_format, directory, outdirectory, cluster_prestep=True):
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	scATAC_runSyntheticCount = robjects.globalenv['scATAC_runSyntheticCount']
	if cluster_prestep == True:
		scATAC_runSyntheticCount(samplename, sample_format, directory, outdirectory, cluster_prestep = 1)
	else:
		scATAC_runSyntheticCount(samplename, sample_format, directory, outdirectory, cluster_prestep = 0)

def scRNA_GenerateSyntheticCount(samplename, sample_format, directory, outdirectory, cluster_prestep=True):
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/SyntheticCountFunctions.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	scRNA_runSyntheticCount = robjects.globalenv['scRNA_runSyntheticCount']
	if cluster_prestep == True:
		scRNA_runSyntheticCount(samplename, sample_format, directory, outdirectory, cluster_prestep = 1)
	else:
		scRNA_runSyntheticCount(samplename, sample_format, directory, outdirectory, cluster_prestep = 0)
