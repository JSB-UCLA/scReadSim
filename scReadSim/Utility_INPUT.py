import csv
import numpy as np
import pysam
import pandas as pd
from collections import defaultdict
from time import process_time
import sys
import subprocess
from tqdm import tqdm
import os

# Test
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/package_development/package_results/20220311_10X_scATACseq_INPUT"
# input_peakfile = "10X_ATAC_chr1_4194444_4399104.INPUT.peaks.bed"
# genome_size_file = "/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/mm10.chrom.sizes.removed"
# feature_file = outdirectory + input_peakfile
# comple_faeture_peakfile = input_comple_peakfile

def ComplementFeature(feature_file, comple_faeture_peakfile, genome_size_file, outdirectory, bedtools_directory):
	genome_size_df = pd.read_csv(genome_size_file, header=None, delimiter="\t")
	input_peak_df = pd.read_csv(feature_file, header=None, delimiter="\t")
	chromosomes_coverd = input_peak_df[0].unique()
	# Select chromosize according to input bed file
	search_string_chr = '|'.join(chromosomes_coverd)
	cmd = "cat %s | grep -Ew '%s' > %s/genome_size_selected.txt" % (genome_size_file, search_string_chr, outdirectory)
	output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	    print('[ERROR] Fail to extract corresponding chromosomes from genome size file:\n', error.decode())
	complement_cmd = "%s/bedtools complement -i %s -g %s/genome_size_selected.txt > %s/%s" % (bedtools_directory, feature_file, outdirectory, outdirectory, comple_faeture_peakfile)
	output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	if error:
	    print('[ERROR] Fail to create complementary feature set:\n', error.decode())
	print('Done!')   

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def fragment_length(open_peak):
    output = np.asarray([int(x[2]) - int(x[1]) for x in open_peak])
    return output


def match_peak(true_peakfile, ref_peakfile, outdirectory, assignment_file):
	with open(true_peakfile) as true_peak_file:
	    reader = csv.reader(true_peak_file, delimiter="\t")
	    true_peak_set = np.asarray(list(reader))
	with open(ref_peakfile) as ref_peak_file:
	    reader = csv.reader(ref_peak_file, delimiter="\t")
	    ref_peak_set = np.asarray(list(reader))
	ref_peak_fraglen = fragment_length(ref_peak_set)
    # ref_peak_fraglen.view('i8,i8,i8').sort(order=['f0'], axis=0)
    with open("%s/%s" % (outdirectory, assignment_file), 'w') as outsfile:
		for true_peak in true_peak_set:
		    true_length = int(true_peak[2]) - int(true_peak[1])
		    idx = find_nearest(ref_peak_fraglen, true_length)
            print("\t".join(true_peak[0:3]) + '\t' + "\t".join(ref_peak_set[idx][0:3]), file=outsfile)
	print('Done!')   

def bam2countmat(cells_barcode_file, bed_directory, bed_file, sam_filename, outdirectory, count_mat_filename):
    cells = pd.read_csv(cells_barcode_file, sep="\t")
    cells = cells.values.tolist()
    cells_barcode = [item[0] for item in cells]
    with open("%s/%s" % (outdirectory, count_mat_filename), 'w') as outsfile:
        samfile = pysam.AlignmentFile(sam_filename, "rb")
        with open("%s/%s" % (bed_directory, bed_file)) as open_peak:
            reader = csv.reader(open_peak, delimiter="\t")
            open_peak = np.asarray(list(reader))
        k = 0
        cellsdic = defaultdict(lambda: [None])
        for cell in cells_barcode:
            cellsdic[cell] = k
            k += 1
        k = 0
        peaksdic = defaultdict(lambda: [None])
        for rec in open_peak:
            rec_name = '_'.join(rec)
            peaksdic[rec_name] = k
            k += 1
        cells_n = len(cells_barcode)
        peaks_n = len(open_peak)
        # marginal_count_vec = [0] * len(open_peak)
        print("Converting count matrix...\n")
        # for rec in open_peak:
        for rec_id in tqdm(range(len(open_peak))):
            rec = open_peak[rec_id]
            rec_name = '_'.join(rec)
            currcounts = [0]*cells_n
            reads = samfile.fetch(rec[3], int(rec[4]), int(rec[5]))
            for read in reads:
                cell = read.qname.split(":")[0].upper()
                if cell in cells_barcode:
                    try:
                        currcounts[cellsdic[cell]] += 1
                    except KeyError:
                        pass
            # marginal_count_vec[rec_id] = sum(currcounts)
            # if sum(currcounts) > 0:
            print(rec_name + "\t" + "\t".join([str(x) for x in currcounts]),file = outsfile)




