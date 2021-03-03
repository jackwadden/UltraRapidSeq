#!/usr/bin/python

import sys
import os
import operator
import pysam
import vcfpy
import datetime
from Bio import SeqIO


####################################################
####################################################
####################################################

def parseVCF(vcf_fn):

    reader = vcfpy.Reader.from_path(vcf_fn)

    # extract variant roi
    for record in reader:
        return record

    return

####################################################
# expects 1) fastq file with timestamps, and 2) bam file as input, 3) VCF file describing variant

if len(sys.argv) != 4:
    print("usage: python3 compute_vaf_over_time.py <fastq> <bamfile> <vcffile>")
    sys.exit(1)

combined_fastq_fn = sys.argv[1]
combined_fastq_bam_fn = sys.argv[2]
vcf_fn = sys.argv[3]

####################################################
# extract roi from VCF file
vcf_record = parseVCF(vcf_fn)


####################################################
# for each read, identify the basecall over the target
bamfile = pysam.AlignmentFile(combined_fastq_bam_fn, 'r')

pileup_iter = bamfile.pileup(vcf_record.CHROM, vcf_record.POS - 1, vcf_record.POS,
                             truncate=True,
                             stepper="samtools",
                             min_base_quality=0,
                             ignore_overlaps=False,
                             ignore_orphans=False)

readid_to_base_map = dict()

for column in pileup_iter:

    if(column.reference_pos == vcf_record.POS - 1):
        for read in column.pileups:

            if read.is_del:
                readid_to_base_map[read.alignment.query_name] = 'del'
            elif read.is_refskip:
                readid_to_base_map[read.alignment.query_name] = 'ins'
            else:
                readid_to_base_map[read.alignment.query_name] = read.alignment.query_sequence[read.query_position]
                # query position is None if is_del or is_refskip is set.
                #print ('\tbase in read %s = %s' % (read.alignment.query_name, read.alignment.query_sequence[read.query_position]))

####################################################
# sort fastqs based on time

record_time_list = list()

print("# Reading FASTQ records...")
for record in SeqIO.parse(combined_fastq_fn, "fastq"):

        # record.description holds the entire first line of the read
        # split by space it gives name, runid, sampleid, read, channel, start_time, barcode
        read_info = record.description.split()
        start_time = read_info[5].split('=')[1]

        # 2019-07-10T20:53:45Z
        parsed_time = datetime.datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%SZ")
        time_record_entry = (parsed_time, record)
        record_time_list.append(time_record_entry)
        #print(record_time_list[len(record_time_list)-1])

print("   Found {} reads.".format(len(record_time_list)))
    
def get_first(elem):
    return elem[0]

print("# Sorting based on read time stamps...")
sorted_record_time_list = sorted(record_time_list, key=get_first )
print("   Done.")

#### emit sorted fastq?
sorted_records = list()
for time, record in sorted_record_time_list:
    sorted_records.append(record)
    
with open("time_sorted_reads.fastq", "w") as output_handle:
    SeqIO.write(sorted_records, output_handle, "fastq")


####################################################
# print mut/(mut+wt) over time

start_time = sorted_record_time_list[0][0]
mut_total = 0
wt_total = 0
vaf = 0.0

# extract ref/alt values
REF = vcf_record.REF.upper()
for alt in vcf_record.ALT:
    ALT = alt.value.upper()


for read in sorted_record_time_list:

    #print(read)
    time = read[0]
    delta = time - start_time
    seconds_elapsed = delta.total_seconds()
    record = read[1]
    read_id = record.id
    basecall = 'X'
    
    if read_id in readid_to_base_map.keys():
        basecall = readid_to_base_map[read_id]
        if basecall == ALT:
            mut_total = mut_total + 1
        elif basecall == REF:
            wt_total = wt_total + 1
        vaf = float(mut_total)/(float(mut_total)+float(wt_total))
    
    print(seconds_elapsed, read_id, basecall, mut_total, wt_total, mut_total + wt_total, vaf)


