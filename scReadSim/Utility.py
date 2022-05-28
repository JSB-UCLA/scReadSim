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

def CallPeak(macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre):
    """Perform peak calling using MACS3 

    Parameters
    ----------
    macs3_directory: `str`
        Path to software MACS3.
    INPUT_bamfile: `str`
        Input BAM file for anlaysis.
    outdirectory: `str`
        Output directory of peak calling.
    MACS3_peakname_pre: `str`
        Base name of peak calling results for MACS3.

    """
    macs_cmd = "%s/macs3 callpeak -f BAMPE -t %s -g mm -n %s/%s -B -q 0.01 --outdir %s" % (macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre, outdirectory)
    output, error = subprocess.Popen(macs_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to call peaks:\n', error.decode())


def ExtractBAMCoverage(INPUT_bamfile, samtools_directory, outdirectory):
    """Examine the covered chromosome names for the input bam file.

    Parameters
    ----------
    INPUT_bamfile: `str`
        Directory of input BAM file.
    samtools_directory: `str`
        Directory of software samtools.
    outdirectory: `str`
        Output directory.
    
    Return
    ------
    chromosomes_coverd: `list`
        List of chromosome names that the input bam files covers.
    """
    cmd = "%s/samtools idxstats %s > %s/bam.stats.txt" % (samtools_directory, INPUT_bamfile, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to rerturn the index BAM information:\n', error.decode())    
    bamstats = pd.read_csv("%s/bam.stats.txt" % outdirectory, header=None, delimiter="\t").to_numpy()
    chromosomes_coverd = bamstats[np.nonzero(bamstats[:,2])[0],0].tolist()
    return chromosomes_coverd


def scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_size_file, ref_peakfile, ref_comple_peakfile, MACS3_peakname_pre):
    """Create the foreground and background feature set for the input scATAC-seq bam file.

    Parameters
    ----------
    INPUT_bamfile: `str`
        Directory of input BAM file.
    samtools_directory: `str`
        Directory of software samtools.
    bedtools_directory: `str`
        Directory of software bedtools.
    outdirectory: `str`
        Output directory.
    genome_size_file: `str`
        Directory of Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicates the size.
    ref_peakfile: `str`
        Specify the base name of output foreground feature bed file.
    ref_comple_peakfile: `str`
        Specify the base name of output background feature bed file.    
    MACS3_peakname_pre: `str`
        Base name of peak calling results for MACS3.
    """
    cmd = "%s/bedtools sort -i %s/%s_peaks.narrowPeak | %s/bedtools merge  > %s/%s" % (bedtools_directory, outdirectory, MACS3_peakname_pre, bedtools_directory, outdirectory, ref_peakfile)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to create feature set:\n', error.decode())
    chromosomes_coverd = ExtractBAMCoverage(INPUT_bamfile, samtools_directory, outdirectory)
    search_string_chr = '|'.join(chromosomes_coverd)
    cmd = "cat %s | grep -Ew '%s' > %s/genome_size_selected.txt" % (genome_size_file, search_string_chr, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to extract gene regions from genome annotation file:\n', error.decode())
    complement_cmd = "%s/bedtools complement -i %s/%s -g %s/genome_size_selected.txt > %s/%s" % (bedtools_directory, outdirectory, ref_peakfile, outdirectory, outdirectory, ref_comple_peakfile)
    output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to create complementary feature set:\n', error.decode())
    print('Done!\n')   


def scRNA_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_annotation, genome_size_file, ref_peakfile, ref_comple_peakfile):
    """Create the foreground and background feature set for the input scRNA-seq bam file.

    Parameters
    ----------
    INPUT_bamfile: `str`
        Input BAM file.
    samtools_directory: `str`
        Path to software `samtools`.
    bedtools_directory: `str`
        Path to software `bedtools`.
    outdirectory: `str`
        Specify the output directory of the features files.
    genome_annotation: `str`
        Genome annotation file for the reference genome that the input BAM aligned on or the synthetic BAM should align on.
    genome_size_file: `str`
        Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicates the size.
    ref_peakfile: `str`
        Specify the name of output foreground feature bed file.
    ref_comple_peakfile: `str`
        Specify the name of output background feature bed file.    
    """
    genome_size_df = pd.read_csv(genome_size_file, header=None, delimiter="\t")
    chromosomes_coverd = ExtractBAMCoverage(INPUT_bamfile, samtools_directory, outdirectory)
    search_string_chr = '|'.join(chromosomes_coverd)
    cmd = "cat %s | grep -Ew '%s' > %s/genome_size_selected.txt" % (genome_size_file, search_string_chr, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to extract corresponding chromosomes from genome size file:\n', error.decode())
    cmd = """awk -F"\t" '$3=="gene"' %s | cut -f1,4,5 > %s/gene_region.bed""" % (genome_annotation, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to extract gene regions from genome annotation file:\n', error.decode())
    cmd = "%s/bedtools sort -i %s/gene_region.bed | %s/bedtools merge | grep -Ew '%s' > %s/%s" % (bedtools_directory, outdirectory, bedtools_directory, search_string_chr, outdirectory, ref_peakfile)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to create feature set:\n', error.decode())
    os.system("rm %s/gene_region.bed" % outdirectory)
    complement_cmd = "%s/bedtools complement -i %s/%s -g %s/genome_size_selected.txt > %s/%s" % (bedtools_directory, outdirectory, ref_peakfile, outdirectory, outdirectory, ref_comple_peakfile)
    output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to create complementary feature set:\n', error.decode())
    print('Done!\n')   


def bam2countmat(cells_barcode_file, bed_file, INPUT_bamfile, outdirectory, count_mat_filename):
    """Construct count matrix.

    Parameters
    ----------
    cells_barcode_file: `str`
        Cell barcode file corresponding to the input BAM file.
    bed_file: `str`
        Features bed file to generate the count matrix.
    INPUT_bamfile: `str`
        Input BAM file for anlaysis.
    outdirectory: `str`
        Specify the output directory of the count matrix file.
    count_mat_filename: `str`
        Specify the base name of output count matrix.
    """
    cells = pd.read_csv(cells_barcode_file, sep="\t")
    cells = cells.values.tolist()
    cells_barcode = [item[0] for item in cells]
    with open("%s/%s.txt" % (outdirectory, count_mat_filename), 'w') as outsfile:
        samfile = pysam.AlignmentFile(INPUT_bamfile, "rb")
        with open(bed_file) as open_peak:
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
            reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
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


def bam2countmat_INPUT(cells_barcode_file, assignment_file, INPUT_bamfile, outdirectory, count_mat_filename):
    """Construct count matrix for task with user input features set. 

    Parameters
    ----------
    cells_barcode_file: `str`
        Cell barcode file corresponding to the input BAM file.
    assignment_file: `str`
        Features assignment file output by function `match_peak`.
    INPUT_bamfile: `str`
        Input BAM file for anlaysis.
    outdirectory: `str`
        Specify the output directory of the count matrix file.
    count_mat_filename: `str`
        Specify the base name of output count matrix.
    """
    cells = pd.read_csv(cells_barcode_file, sep="\t")
    cells = cells.values.tolist()
    cells_barcode = [item[0] for item in cells]
    with open("%s/%s.txt" % (outdirectory, count_mat_filename), 'w') as outsfile:
        samfile = pysam.AlignmentFile(INPUT_bamfile, "rb")
        with open(assignment_file) as open_peak:
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


def find_nearest(array, value):
    """Find the index of peak from `array` with a length closest to `value`.

    Parameters
    ----------
    array: `numpy.array`
        Two-column array of peaks indicating the starting and ending positions.
    value: `int`
        Integer indicating the target peak length.

    Returns
    -------
    idx: `int`
        Index of the peak with a length closest to `value`.
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def match_peak(true_peakfile, ref_peakfile, outdirectory, assignment_file):
    """Find the reference features for the given features `true_peakfile` according to the similarity of peak length. The assignment would be stored as `assignment_file` within `outdirectory`.

    Parameters
    ----------
    true_peakfile: `str`
        Directory of the target features bed file.
    ref_peakfile: `str`
        Directory of the reference features bed file.
    outdirectory: `str`
        Output directory.
    assignment_file: `str`
        Specify the name of peak assignment file.
    """
    with open(true_peakfile) as true_peak_file:
        reader = csv.reader(true_peak_file, delimiter="\t")
        true_peak_set = np.asarray(list(reader))
    with open(ref_peakfile) as ref_peak_file:
        reader = csv.reader(ref_peak_file, delimiter="\t")
        ref_peak_set = np.asarray(list(reader))
    ref_peak_fraglen = np.asarray([int(x[2]) - int(x[1]) for x in ref_peak_set])
    with open("%s/%s" % (outdirectory, assignment_file), 'w') as outsfile:
        for true_peak in true_peak_set:
            true_length = int(true_peak[2]) - int(true_peak[1])
            idx = find_nearest(ref_peak_fraglen, true_length)
            print("\t".join(true_peak[0:3]) + '\t' + "\t".join(ref_peak_set[idx][0:3]), file=outsfile)
    print('Done!')


def ComplementFeature(feature_file, comple_feature_peakfile, genome_size_file, outdirectory, bedtools_directory):
    """Create background feature set given the user input foreground feature set `feature_file`.

    Parameters
    ----------
    feature_file: `str`
        Input feature set
    comple_feature_peakfile: `str`
        Specify the base name of output background feature bed file.    
    genome_size_file: `str`
        Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicating the size.
    outdirectory: `str`
        Output directory.
    bedtools_directory: `str`
        Path to software bedtools.
    """
    genome_size_df = pd.read_csv(genome_size_file, header=None, delimiter="\t")
    input_peak_df = pd.read_csv(feature_file, header=None, delimiter="\t")
    chromosomes_coverd = input_peak_df[0].unique()
    # Select chromosize according to input bed file
    search_string_chr = '|'.join(chromosomes_coverd)
    cmd = "cat %s | grep -Ew '%s' > %s/genome_size_selected.txt" % (genome_size_file, search_string_chr, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to extract corresponding chromosomes from genome size file:\n', error.decode())
    complement_cmd = "%s/bedtools complement -i %s -g %s/genome_size_selected.txt > %s/%s" % (bedtools_directory, feature_file, outdirectory, outdirectory, comple_feature_peakfile)
    output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to create complementary feature set:\n', error.decode())
    print('Done!')   
















