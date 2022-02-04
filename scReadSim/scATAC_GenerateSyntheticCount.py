## This script takes use of Rscripts to generate synthetic count matrix.
import numpy as np
import rpy2.robjects as robjects
# import rpy2.robjects.numpy2ri as np2r
# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri
import time
import pkg_resources
import os

def GenerateSyntheticCount(samplename, sample_format, directory, outdirectory, cluster_prestep=True):
	r = robjects.r
	rscript_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Rscript/scATAC_SyntheticCountFunctions.R')
	print(rscript_dir)
	# rscript_dir = pkg_resources.resource_stream(__name__, 'Rscript/scATAC_SyntheticCountFunctions.R').read().decode()
	r['source'](rscript_dir)
	scReadSim_runSyntheticCount = robjects.globalenv['scReadSim_runSyntheticCount']
	if cluster_prestep == True:
		scReadSim_runSyntheticCount(samplename, sample_format, directory, outdirectory, cluster_prestep = 1)
	else:
		scReadSim_runSyntheticCount(samplename, sample_format, directory, outdirectory, cluster_prestep = 0)
