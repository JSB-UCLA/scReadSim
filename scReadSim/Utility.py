import csv
import numpy as np
import pysam
import pandas as pd
from collections import defaultdict
from time import process_time
import sys
import subprocess
from tqdm import tqdm

def CallPeak(macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre):
    macs_cmd = "%s/macs3 callpeak -f BAMPE -t %s -g mm -n %s/%s -B -q 0.01 --outdir %s" % (macs3_directory, INPUT_bamfile, outdirectory, MACS3_peakname_pre, outdirectory)
    output, error = subprocess.Popen(macs_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to call peaks:\n', error.decode())


def scATAC_CreateFeatureSets(bedtools_directory, outdirectory, genome_file, ref_peakfile, ref_comple_peakfile, MACS3_peakname_pre):
    cmd = "%s/bedtools sort -i %s/%s_peaks.narrowPeak | %s/bedtools merge  > %s/%s" % (bedtools_directory, outdirectory, MACS3_peakname_pre, bedtools_directory, outdirectory, ref_peakfile)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to create feature set:\n', error.decode())
    complement_cmd = "%s/bedtools complement -i %s/%s -g %s > %s/%s" % (bedtools_directory, outdirectory, ref_peakfile, genome_file, outdirectory, ref_comple_peakfile)
    output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to create complementary feature set:\n', error.decode())
    print('Done!\n')   

def scRNA_CreateFeatureSets(bedtools_directory, outdirectory, genome_annotation, genome_size, ref_peakfile, ref_comple_peakfile):
    cmd = "awk -F\"\t\" '$3==\"gene\"' %s | cut -f1,4,5 > %s/gene_region.bed" % (genome_annotation, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to extract gene regions from genome annotation file:\n', error.decode())
    cmd = "%s/bedtools sort -i %s/gene_region.bed | %s/bedtools merge  > %s/%s" % (bedtools_directory, outdirectory, bedtools_directory, outdirectory, ref_peakfile)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to create feature set:\n', error.decode())
    os.system("rm %s/gene_region.bed", outdirectory)
    complement_cmd = "%s/bedtools complement -i %s/%s -g %s > %s/%s" % (bedtools_directory, outdirectory, ref_peakfile, genome_size, outdirectory, ref_comple_peakfile)
    output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to create complementary feature set:\n', error.decode())
    print('Done!\n')   

# ## Construct expression matrix for single paired read
# def read_pair_generator(bam, contig_str=None,start_int = None,stop_int=None):
#     """
#     Generate read pairs in a BAM file or within a region string.
#     Reads are added to read_dict until a pair is found.
#     """
#     read_dict = defaultdict(lambda: [None, None])
#     for read in bam.fetch(contig=contig_str, start=start_int, stop=stop_int):
#         if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_duplicate:
#             continue
#         qname = read.query_name
#         if qname not in read_dict:
#             if read.is_read1:
#                 read_dict[qname][0] = read
#             else:
#                 read_dict[qname][1] = read
#         else:
#             if read.is_read1:
#                 yield read, read_dict[qname][1]
#             else:
#                 yield read_dict[qname][0], read
#             del read_dict[qname]

# def read_pair_generator_supp(bam, contig_str=None,start_int = None,stop_int=None):
#     """
#     Generate read pairs in a BAM file or within a region string.
#     Reads are added to read_dict until a pair is found.
#     """
#     read_dict = defaultdict(lambda: [None, None])
#     for read in bam.fetch(contig=contig_str, start=start_int, stop=stop_int):
#         if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_duplicate:
#             continue
#         qname = read.query_name
#         if qname not in read_dict:
#             if read.is_read1:
#                 read_dict[qname][0] = read
#             else:
#                 read_dict[qname][1] = read
#         else:
#             del read_dict[qname]
#     single_reads = []
#     for qname in read_dict:
#         if read_dict[qname][0] == None:
#             single_reads.append(read_dict[qname][1])
#         else:
#             single_reads.append(read_dict[qname][0])
#     return single_reads

# def bam2countmat_countPairsandSingle(cells_barcode_file, bed_directory, bed_file, sam_filename, outdirectory, count_mat_filename):
#     cells = pd.read_csv(cells_barcode_file, sep="\t")
#     cells = cells.values.tolist()
#     cells_barcode = [item[0] for item in cells]
#     with open("%s/%s" % (outdirectory, count_mat_filename), 'w') as outsfile:
#         samfile = pysam.AlignmentFile(sam_filename, "rb")
#         with open("%s/%s" % (bed_directory, bed_file)) as open_peak:
#             reader = csv.reader(open_peak, delimiter="\t")
#             open_peak = np.asarray(list(reader))
#         k = 0
#         cellsdic = defaultdict(lambda: [None])
#         for cell in cells_barcode:
#             cellsdic[cell] = k
#             k += 1
#         k = 0
#         peaksdic = defaultdict(lambda: [None])
#         print("Converting count matrix...\n")
#         for rec in open_peak:
#             rec_name = '_'.join(rec)
#             peaksdic[rec_name] = k
#             k += 1
#         cells_n = len(cells_barcode)
#         peaks_n = len(open_peak)
#         # for rec in open_peak:
#         for rec_id in tqdm(range(len(open_peak))):
#             rec = open_peak[rec_id]
#             rec_name = '_'.join(rec)
#             currcounts = [0]*cells_n
#             for read1, read2 in read_pair_generator(samfile, rec[0], int(rec[1]), int(rec[2])):
#                 read1_str = str(read1).split("\t")
#                 read2_str = str(read2).split("\t")
#                 if read1.has_tag('CB:Z:'):
#                     cell = read1.get_tag('CB:Z:')
#                     # cell = read1_str[0].split(':')[0]
#                     if cell in cells_barcode:
#                         try:
#                             currcounts[cellsdic[cell]] += 2
#                         except KeyError:
#                             pass
#             for read in read_pair_generator_supp(samfile, rec[0], int(rec[1]), int(rec[2])):
#                 read_str = str(read).split("\t")
#                 if read.has_tag('CB:Z:'):
#                     cell = read.get_tag('CB:Z:')
#                     # cell = read_str[0].split(':')[0]
#                     if cell in cells_barcode:
#                         try:
#                             currcounts[cellsdic[cell]] += 1
#                         except KeyError:
#                             pass
#             # print(sum(currcounts))
#             # if sum(currcounts) > 0:
#             print(rec_name + "\t" + "\t".join([str(x) for x in currcounts]),file = outsfile)


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
            # if sum(currcounts) > 0:
            print(rec_name + "\t" + "\t".join([str(x) for x in currcounts]),file = outsfile)

# # Test
# cells_barcode_file = INPUT_cells_barcode_file
# bed_directory = outdirectory
# bed_file = ref_peakfile
# sam_filename = INPUT_bamfile
# outdirectory = outdirectory
# count_mat_filename = "%s.%s" % (count_mat_filename_new, count_mat_format)

# def bam2countmat_new(cells_barcode_file, bed_directory, bed_file, sam_filename, outdirectory, count_mat_filename):
#     cells = pd.read_csv(cells_barcode_file, sep="\t")
#     cells = cells.values.tolist()
#     cells_barcode = [item[0] for item in cells]
#     nonzero_peak_list = list()
#     nonzero_bed_file = "nonzero_%s" % bed_file
#     nonzero_peak_list = list()
#     with open("%s/%s" % (outdirectory, count_mat_filename), 'w') as outsfile:
#         samfile = pysam.AlignmentFile(sam_filename, "rb")
#         with open("%s/%s" % (bed_directory, bed_file)) as open_peak:
#             reader = csv.reader(open_peak, delimiter="\t")
#             open_peak = np.asarray(list(reader))
#         k = 0
#         cellsdic = defaultdict(lambda: [None])
#         for cell in cells_barcode:
#             cellsdic[cell] = k
#             k += 1
#         k = 0
#         peaksdic = defaultdict(lambda: [None])
#         for rec in open_peak:
#             rec_name = '_'.join(rec)
#             peaksdic[rec_name] = k
#             k += 1
#         cells_n = len(cells_barcode)
#         peaks_n = len(open_peak)
#         print("Converting count matrix...\n")
#         # for rec in open_peak:
#         for rec_id in tqdm(range(len(open_peak))):
#             rec = open_peak[rec_id]
#             rec_name = '_'.join(rec)
#             currcounts = [0]*cells_n
#             reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
#             for read in reads:
#                 cell = read.qname.split(":")[0].upper()
#                 if cell in cells_barcode:
#                     try:
#                         currcounts[cellsdic[cell]] += 1
#                     except KeyError:
#                         pass
#             if sum(currcounts) > 0:
#                 nonzero_peak_list.append(rec)
#                 print(rec_name + "\t" + "\t".join([str(x) for x in currcounts]),file = outsfile)
#     with open("%s/%s" % (outdirectory, nonzero_bed_file), 'w') as outsfile2:
#         for rec in nonzero_peak_list:
#             print("\t".join(rec),file = outsfile2)



















