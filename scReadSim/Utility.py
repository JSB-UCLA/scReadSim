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
from joblib import Parallel, delayed
from collections import Counter


def CallPeak(macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre, qval=0.05):
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
    macs_cmd = "%s/macs3 callpeak -f BAMPE -t %s -g mm -n %s/%s -B -q %s --outdir %s" % (macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre, qval, outdirectory)
    output, error = subprocess.Popen(macs_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    print('[MACS3] Call peaks:\n', error.decode())


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


def scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_size_file, macs3_directory, INPUT_peakfile=None, INPUT_nonpeakfile=None, OUTPUT_peakfile=None):
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
    macs3_directory: `str`
        Path to software MACS3.
    INPUT_peakfile: `str` (default: None)
        Directory of user-specified input peak file.
    INPUT_nonpeakfile: `str` (default: None)
        Directory of user-specified input non-peak file.
    OUTPUT_peakfile: `str` (default: None)
        Directory of user-specified output peak file. Synthetic scATAC-seq reads will be generated taking `OUTPUT_peakfile` as ground truth peaks. Note that `OUTPUT_peakfile` does not name the generated feature files by function `scATAC_CreateFeatureSets`.
    """
    chromosomes_coverd = ExtractBAMCoverage(INPUT_bamfile, samtools_directory, outdirectory)
    search_string_chr = '|'.join(chromosomes_coverd)
    cmd = "cat %s | grep -Ew '%s' > %s/genome_size_selected.txt" % (genome_size_file, search_string_chr, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to extract gene regions from genome annotation file:\n', error.decode())
    # Input peaks and non-peaks
    if INPUT_peakfile is None or INPUT_nonpeakfile is None:
        # Define peaks and non-peaks using MACS3
        peakfile = outdirectory + "/" + "scReadSim.MACS3.peak.bed"
        nonpeakfile = outdirectory + "/" + "scReadSim.MACS3.nonpeak.bed"
        # Call peaks
        print("[scReadSim] No Input Peaks and Non-Peaks.")
        print("[scReadSim] Generate Peaks and Non-Peaks using MACS3 Instead.")
        print("[scReadSim] Generating Peaks using MACS3...")
        CallPeak(macs3_directory, INPUT_bamfile, outdirectory, "scReadSim_MACS3_Stringent", qval=0.01)
        cmd = "%s/bedtools sort -i %s/scReadSim_MACS3_Stringent_peaks.narrowPeak | %s/bedtools merge  > %s/scReadSim.MACS3.peak.bed" % (bedtools_directory, outdirectory, bedtools_directory, outdirectory)
        output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if error:
            print('[ERROR] Fail to create feature set:\n', error.decode())
        print("[scReadSim] Peaks Generated: %s/scReadSim.MACS3.peak.bed" % outdirectory)
        # Call non-peaks
        print("[scReadSim] Generating Non-Peaks using MACS3...")
        CallPeak(macs3_directory, INPUT_bamfile, outdirectory, "scReadSim_MACS3_LessStringent", qval=0.1)
        cmd = "%s/bedtools sort -i %s/scReadSim_MACS3_LessStringent_peaks.narrowPeak | %s/bedtools merge  > %s/scReadSim_MACS3.nonpeak.tmp.bed" % (bedtools_directory, outdirectory, bedtools_directory, outdirectory)
        output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if error:
            print('[ERROR] Fail to create feature set:\n', error.decode())
        # Calculate the inter peaks as non-peaks
        complement_cmd = "%s/bedtools complement -i %s/scReadSim_MACS3.nonpeak.tmp.bed -g %s/genome_size_selected.txt > %s/scReadSim.MACS3.nonpeak.bed" % (bedtools_directory, outdirectory, outdirectory, outdirectory)
        output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if error:
            print('[ERROR] Fail to create complementary feature set:\n', error.decode())
        print("[scReadSim] Non-Peaks Generated: %s/scReadSim.MACS3.nonpeak.bed" % outdirectory)
    else:
        peakfile = "scReadSim.UserInput.peak.bed"
        nonpeakfile = "scReadSim.UserInput.nonpeak.bed"
        # Merge and sort peak file
        print("[scReadSim] Merging Input Peak File: %s" % INPUT_peakfile)
        cmd = "%s/bedtools sort -i %s | %s/bedtools merge  > %s/%s" % (bedtools_directory, INPUT_peakfile, bedtools_directory, outdirectory, peakfile)
        output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if error:
            print('[ERROR] Fail to sort input peak set:\n', error.decode())
        print("[scReadSim] Merging Input Non-Peak File: %s" % INPUT_nonpeakfile)
        print("[scReadSim] Peaks Generated: %s/%s" % (outdirectory, peakfile))
        # Merge and sort non-peak file
        cmd = "%s/bedtools sort -i %s | %s/bedtools merge  > %s/%s" % (bedtools_directory, INPUT_nonpeakfile, bedtools_directory, outdirectory, nonpeakfile)
        output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if error:
            print('[ERROR] Fail to merge input non-peak feature set:\n', error.decode())
        print("[scReadSim] Non-Peaks Generated: %s/%s" % (outdirectory, nonpeakfile))
    # Union peaks and non-peaks
    print("[scReadSim] Generating Gray Areas...")
    cmd = "cat %s/%s %s/%s | bedtools sort | bedtools merge > %s/Peak_NonPeaks_Union.bed" % (outdirectory, peakfile, outdirectory, nonpeakfile, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    complement_cmd = "%s/bedtools complement -i %s/Peak_NonPeaks_Union.bed -g %s/genome_size_selected.txt > %s/scReadSim.grayareas.bed" % (bedtools_directory, outdirectory, outdirectory, outdirectory)
    output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to create gray area feature set:\n', error.decode())
    print("[scReadSim] Gray Areas Generated: %s/scReadSim.grayareas.bed" % outdirectory)
    print('\n[scReadSim] Created:')
    print('[scReadSim] Peak File: %s/%s' % (outdirectory, peakfile))
    print('[scReadSim] Non-Peak File: %s/%s' % (outdirectory, nonpeakfile))
    print('[scReadSim] Gray Area File: %s/scReadSim.grayareas.bed' % (outdirectory))
    # Output Peak Mode
    if OUTPUT_peakfile is not None:
        output_peakfile = "scReadSim.output.peak.bed"
        output_nonpeakfile = "scReadSim.output.nonpeak.bed"
        print("\n[scReadSim] User-specified Output Peaks Detected: %s" % OUTPUT_peakfile)
        # Merge and sort output peak file
        print("[scReadSim] Merging Output Peak File: %s" % OUTPUT_peakfile)
        cmd = "%s/bedtools sort -i %s | %s/bedtools merge  > %s/%s" % (bedtools_directory, OUTPUT_peakfile, bedtools_directory, outdirectory, output_peakfile)
        output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if error:
            print('[ERROR] Fail to sort output peak set:\n', error.decode())
        print("Generating Output Non-Peaks...")
        complement_cmd = "%s/bedtools complement -i %s -g %s/genome_size_selected.txt > %s/%s" % (bedtools_directory, OUTPUT_peakfile, outdirectory, outdirectory, output_nonpeakfile)
        output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if error:
            print('[ERROR] Fail to create output non-peak feature set:\n', error.decode())
        print('\n[scReadSim] Created:')
        print('[scReadSim] Output Peak File: %s/%s' % (outdirectory, output_peakfile))
        print('[scReadSim] Output Non-Peak File: %s/%s' % (outdirectory, output_nonpeakfile))
    print('[scReadSim] Done!')


def scRNA_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_annotation, genome_size_file):
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
    cmd = "%s/bedtools sort -i %s/gene_region.bed | %s/bedtools merge | grep -Ew '%s' > %s/scReadSim.Gene.bed" % (bedtools_directory, outdirectory, bedtools_directory, search_string_chr, outdirectory)
    print("[scReadSim] Generating Bed File for Genes...")
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to create feature set:\n', error.decode())
    os.system("rm %s/gene_region.bed" % outdirectory)
    print("[scReadSim] Generating Bed File for InterGenes...")
    complement_cmd = "%s/bedtools complement -i %s/scReadSim.Gene.bed -g %s/genome_size_selected.txt > %s/scReadSim.InterGene.bed" % (bedtools_directory, outdirectory, outdirectory, outdirectory)
    output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to create complementary feature set:\n', error.decode())
    print('\n[scReadSim] Created:')
    print('[scReadSim] Gene Bed File: %s/scReadSim.Gene.bed' % (outdirectory))
    print('[scReadSim] InterGene Bed File: %s/scReadSim.InterGene.bed' % (outdirectory))
    print('[scReadSim] Done!')

def countmat_mainloop(rec_id):
    """Construct count vector for each scATAC-seq feature.

    """
    count_array = np.zeros(cells_n, dtype=int) # initialize
    rec = open_peak[rec_id]
    rec_name = '_'.join((rec[0], str(rec[1]), str(rec[2])))
    samfile = pysam.AlignmentFile(INPUT_bamfile_glb, "rb")
    reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))  # question: what about feature intersection or half overlap?
    cell_iter = []
    cell_idx_ls = []
    for read in reads:
        cell_iter.append(read.qname.split(":")[0].upper())
    for cell in cell_iter:
        if cell in cells_barcode:
            cell_idx_ls.append(cells_barcode.index(cell))
    counter = Counter(cell_idx_ls)
    keys = list(counter.keys())
    values = list(counter.values())
    count_array[keys] = values
    count_array_withPeak = np.insert(count_array.astype(str), 0, rec_name)
    return count_array_withPeak

def scATAC_bam2countmat_paral(cells_barcode_file, bed_file, INPUT_bamfile, outdirectory, count_mat_filename, n_cores=1):
    """Construct count matrix for scATAC-seq BAM file.

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
    n_cores: `int` (default: 1)
        Specify the number of cores for parallel computing when generating count matrix.
    """
    cells = pd.read_csv(cells_barcode_file, sep="\t")
    cells = cells.values.tolist()
    # Specify global vars
    global open_peak, cells_n, cells_barcode, INPUT_bamfile_glb
    INPUT_bamfile_glb = INPUT_bamfile
    cells_barcode = [item[0] for item in cells]
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
    print("[scReadSim] Generating read count matrix...\n")
    mat_array = Parallel(n_jobs=n_cores, backend='multiprocessing')(delayed(countmat_mainloop)(rec_id) for rec_id in (range(len(open_peak))))
    para_countmat = np.array(mat_array)
    with open("%s/%s.txt" % (outdirectory, count_mat_filename), 'w') as outsfile:
        for rec_id in tqdm(range(len(open_peak))):
            print("\t".join([str(x) for x in para_countmat[rec_id,:]]),file = outsfile)
    print('[scReadSim] Created:')
    print('[scReadSim] Read Count Matrix: %s/%s.txt' % (outdirectory, count_mat_filename))
    print('[scReadSim] Done!')


def scRNA_UMIcountmat_mainloop(rec_id):
    """Construct count vector for each scRNA-seq feature.

    """
    UMI_currlist  = [["empty UMI"] for _ in range(cells_n)] 
    rec = open_peak[rec_id]
    rec_name = '_'.join((rec[0], str(rec[1]), str(rec[2])))
    samfile = pysam.AlignmentFile(INPUT_bamfile_glb, "rb")
    reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))  
    UMI_iter = []
    cell_idx_ls = []
    for read in reads:
        cell = read.qname.split(":")[0].upper()
        if cell in cells_barcode:
            try:
                if read.has_tag(UMI_tag_glb):
                    UMI = read.get_tag(UMI_tag_glb)
                    UMI_currlist[cellsdic[cell]].append(UMI)
            except KeyError:
                pass
    UMI_count_array = [len(set(UMIs_percell))-1 for UMIs_percell in UMI_currlist]
    UMI_count_array.insert(0,rec_name)
    return UMI_count_array


def scRNA_bam2countmat_paral(cells_barcode_file, bed_file, INPUT_bamfile, outdirectory, count_mat_filename, UMI_modeling=True, UMI_tag = "UB:Z", n_cores=1):
    """Construct read (or UMI) count matrix for scRNA-seq BAM file.

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
        Specify the base name of output read (or UMI) count matrix.
    UMI_modeling: `bool` (default: True)
        Specify whether scReadSim should model UMI count of the input BAM file.
    UMI_tag: `str` (default: 'UB:Z')
        If UMI_modeling is set to True, specify the UMI tag of input BAM file, default value 'UB:Z' is the UMI tag for 10x scRNA-seq.
    n_cores: `int` (default: 1)
        Specify the number of cores for parallel computing when generating count matrix.
    """
    cells = pd.read_csv(cells_barcode_file, sep="\t")
    cells = cells.values.tolist()
    # Specify global vars
    global open_peak, cells_n, cells_barcode, INPUT_bamfile_glb, UMI_tag_glb, cellsdic
    UMI_tag_glb = UMI_tag
    INPUT_bamfile_glb = INPUT_bamfile
    cells_barcode = [item[0] for item in cells]
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
    if UMI_modeling == True:
        print("[scReadSim] UMI Mode Detected.")
        print("[scReadSim] Generating UMI Count Matrix...")
        UMI_countmat_array = Parallel(n_jobs=n_cores, backend='multiprocessing')(delayed(scRNA_UMIcountmat_mainloop)(rec_id) for rec_id in (range(len(open_peak))))
        print("[scReadSim] Generated UMI Count Matrix.")
        UMI_countmat = np.array(UMI_countmat_array)
        print("[scReadSim] Writing UMI Count Matrix TXT File...")
        with open("%s/%s.txt" % (outdirectory, count_mat_filename), 'w') as outsfile:
            for rec_id in tqdm(range(len(open_peak))):
                print("\t".join([str(x) for x in UMI_countmat[rec_id,:]]),file = outsfile)
        print("[scReadSim] Created:")
        print("[scReadSim] UMI Count Matrix %s.txt" % count_mat_filename)
    else:
        print("[scReadSim] Detected that UMI Mode Is Off.")
        print("[scReadSim] Generating Read Count Matrix...")
        read_countmat_array = Parallel(n_jobs=n_cores, backend='multiprocessing')(delayed(countmat_mainloop)(rec_id) for rec_id in (range(len(open_peak))))
        read_countmat = np.array(read_countmat_array)
        print("[scReadSim] Writing Read Count Matrix TXT File...")
        with open("%s/%s.txt" % (outdirectory, count_mat_filename), 'w') as outsfile:
            for rec_id in tqdm(range(len(open_peak))):
                print("\t".join([str(x) for x in read_countmat[rec_id,:]]),file = outsfile)
        print("[scReadSim] Created:")
        print("[scReadSim] Read count matrix %s.txt" % count_mat_filename)
    print("[scReadSim] Done.\n")


def scATAC_bam2countmat_OutputPeak(cells_barcode_file, assignment_file, INPUT_bamfile, outdirectory, count_mat_filename):
    """Construct count matrix for task with user input features set. 

    Parameters
    ----------
    cells_barcode_file: `str`
        Cell barcode file corresponding to the input BAM file.
    assignment_file: `str`
        Features mapping file output by function `Utility.FeatureMapping`.
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
        print("[scReadSim] Generating count matrix...")
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


def find_nearest_peak(array, value, k, ref_read_density):
    """Find the index of largest-read-density input peak from `array` with a length similar to `value`.

    """
    array = np.asarray(array)
    # idx = (np.abs(array - value)).argmin()   
    idx = (np.abs(array - value)).argpartition(k)
    idx_ReadDensity = ref_read_density[idx[:k]].argmax()
    final_idx = idx[idx_ReadDensity]
    return final_idx


def find_nearest_nonpeak(array, value, k, ref_read_density):
    """Find the index of largest-read-density input non-peak from `array` with a length similar to `value`.

    """
    array = np.asarray(array)
    # idx = (np.abs(array - value)).argmin()   
    idx = (np.abs(array - value)).argpartition(k)
    idx_ReadDensity = ref_read_density[idx[:k]].argmin()
    final_idx = idx[idx_ReadDensity]
    return final_idx


def fragment_length(open_peak):
	output = np.asarray([int(x[2]) - int(x[1]) for x in open_peak])
	return output


def match_peak(MarginalCountList, output_peaks, input_peaks, outdirectory, assignment_file, n_top):
    """Map from output peaks to input peaks according to the similarity of peak length. 

    """
    with open("%s" % (output_peaks)) as true_peak_file:
        reader = csv.reader(true_peak_file, delimiter="\t")
        true_peak_set = np.asarray(list(reader))
    with open("%s" % (input_peaks)) as ref_peak_file:
        reader = csv.reader(ref_peak_file, delimiter="\t")
        ref_peak_set = np.asarray(list(reader))
    ref_marginal_count = np.asarray(MarginalCountList)
    ref_peak_fraglen = fragment_length(ref_peak_set)
    ref_read_density = ref_marginal_count.astype(int) / ref_peak_fraglen
    # ref_peak_fraglen.view('i8,i8,i8').sort(order=['f0'], axis=0)
    n_top_select = min(len(ref_marginal_count)-1, n_top)
    with open("%s/%s" % (outdirectory, assignment_file), 'w') as outsfile:
        for true_peak_id in tqdm(range(len(true_peak_set))):
            true_peak = true_peak_set[true_peak_id]
            true_length = int(true_peak[2]) - int(true_peak[1])
            idx = find_nearest_peak(ref_peak_fraglen, true_length, n_top_select, ref_read_density)
            print("\t".join(true_peak[0:3]) + '\t' + "\t".join(ref_peak_set[idx][0:3]), file=outsfile)


def match_nonpeak(MarginalCountList ,output_nonpeaks, input_nonpeaks, outdirectory, assignment_file, n_top):
    """Map from output non-peaks to input non-peaks according to the similarity of peak length. 

    """
    with open("%s" % (output_nonpeaks)) as true_peak_file:
        reader = csv.reader(true_peak_file, delimiter="\t")
        true_peak_set = np.asarray(list(reader))
    with open("%s" % (input_nonpeaks)) as ref_peak_file:
        reader = csv.reader(ref_peak_file, delimiter="\t")
        ref_peak_set = np.asarray(list(reader))
    ref_marginal_count = np.asarray(MarginalCountList)
    ref_peak_fraglen = fragment_length(ref_peak_set)
    ref_read_density = ref_marginal_count.astype(int) / ref_peak_fraglen
    n_top_select = min(len(ref_marginal_count)-1, n_top)
    # ref_peak_fraglen.view('i8,i8,i8').sort(order=['f0'], axis=0)
    with open("%s/%s" % (outdirectory, assignment_file), 'w') as outsfile:
        for true_peak_id in tqdm(range(len(true_peak_set))):
            true_peak = true_peak_set[true_peak_id]
            true_length = int(true_peak[2]) - int(true_peak[1])
            idx = find_nearest_nonpeak(ref_peak_fraglen, true_length, n_top_select, ref_read_density)
            print("\t".join(true_peak[0:3]) + '\t' + "\t".join(ref_peak_set[idx][0:3]), file=outsfile)


def bam2MarginalCount(bed_file, sam_filename):
    """Count read overlapped for each feature from input bed file. 

    """
    samfile = pysam.AlignmentFile(sam_filename, "rb")
    with open("%s" % (bed_file)) as open_peak:
        reader = csv.reader(open_peak, delimiter="\t")
        open_peak = np.asarray(list(reader))
    k = 0
    peaksdic = defaultdict(lambda: [None])
    for rec in open_peak:
        rec_name = '_'.join(rec)
        peaksdic[rec_name] = k
        k += 1
    peaks_n = len(open_peak)
    MarginalCountList = np.empty((peaks_n), dtype="int")
    print("[scReadSim] Converting marginal count vector...")
    # for rec in open_peak:
    for rec_id in tqdm(range(len(open_peak))):
        rec = open_peak[rec_id]
        rec_name = '_'.join(rec)
        currcounts =  0
        reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
        for read in reads:
            currcounts += 1
        # if sum(currcounts) > 0:
        MarginalCountList[rec_id] = int(currcounts)
        # print(rec_name + "\t" + str(currcounts),file = outsfile)
    return MarginalCountList


def FeatureMapping(INPUT_bamfile, input_peaks, input_nonpeaks, output_peaks, output_nonpeaks, outdirectory, assignment_peak_file, assignment_nonpeak_file, n_top=50):
    """Obtain mappings between input and output peaks, input and output non-peaks. The mappings are output as `assignment_peak_file` and `assignment_nonpeak_file` within `outdirectory`.  

    Parameters
    ----------
    INPUT_bamfile: `str`
        Input BAM file for anlaysis.
    input_peaks: `str`
        BED file of user specified (or generated by scReadSim+MACS3) input peaks.
    input_nonpeaks: `str`
        BED file of user specified (or generated by scReadSim+MACS3) input nonpeaks.
    output_peaks: `str`
        BED file of user specified output peaks.
    output_nonpeaks: `str`
        BED file of output non-peaks, generated by function 'scATAC_CreateFeatureSets'.
    outdirectory: `str`
        Output directory.
    assignment_peak_file: `str`
        Specify the name of peak mapping file.
    assignment_nonpeak_file: `str`
        Specify the name of nonpeak mapping file.
    n_top: 'int'(default: 50)
        Specify the number of input peaks (or non-peaks) with the most similar length as the candidate mapped input peaks (or non-peaks) for each the output peak (or non-peak). From the candidate input peaks (or non-peaks), scReadSim further selects the one with largest read density for peak mapping (smallest read density for non-peak mapping).
    """
    peak_MarginalCountList = bam2MarginalCount(input_peaks, INPUT_bamfile)
    print("[scReadSim] Mapping Input Peaks and Output Peaks...")
    match_peak(peak_MarginalCountList, output_peaks, input_peaks, outdirectory, assignment_peak_file, n_top)
    nonpeak_MarginalCountList = bam2MarginalCount(input_nonpeaks, INPUT_bamfile)
    print("[scReadSim] Mapping Input Non-Peaks and Output Non-Peaks...")
    match_nonpeak(nonpeak_MarginalCountList, output_nonpeaks, input_nonpeaks, outdirectory, assignment_nonpeak_file, n_top)
    print('[scReadSim] Created:')
    print('[scReadSim] Mapping File between Input and Output Peaks: %s/%s.txt' % (outdirectory, assignment_peak_file))
    print('[scReadSim] Mapping File between Input and Output Peaks: %s/%s.txt' % (outdirectory, assignment_nonpeak_file))
    print('[scReadSim] Done!')








