import sys, os, argparse, collections

import pysam
import vcfpy
##########
def rc(base):

    if base == 'A':
        return 'T'
    elif base =='T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    else:
        return ''
    

##########
def is_in_target_region(vcf_record, entry):

    target_width = 1000

    
    #print(entry.reference_name, entry.reference_start, entry.reference_end, vcf_record.POS, vcf_record.CHROM)
    
    if entry.reference_name == vcf_record.CHROM:
        if entry.reference_start > (vcf_record.POS - target_width) and entry.reference_start < (vcf_record.POS + target_width) or entry.reference_end < (vcf_record.POS + target_width) and entry.reference_end > (vcf_record.POS - target_width):
            return True

    return False

def covers_target_locus(vcf_record, entry):
    if entry.reference_name == vcf_record.CHROM:
        if entry.reference_start <= vcf_record.POS and entry.reference_end > vcf_record.POS:
            return True

    return False
    

##########

if len(sys.argv) != 3:
    print("usage: python3 bam_classifier.py <bamfile> <vcffile>")
    exit(1)
    
bam_fn = sys.argv[1]
vcf_fn = sys.argv[2]

# parse samfile
bamfile = pysam.AlignmentFile(bam_fn, 'rb')

# parse target vcf file
vcf_reader = vcfpy.Reader.from_path(vcf_fn)
vcf_record = None
for record in vcf_reader:
    vcf_record = record

vcf_chrom = vcf_record.CHROM
vcf_locus = vcf_record.POS
wt_base = vcf_record.REF
mut_base = vcf_record.ALT[0].value

# for each entry, categorize as "target, fragment, background, unmapped"
read_classification_list = list()
read_classification_counters = {'target' : 0, 'fragment': 0, 'background' : 0, 'unmapped' : 0}

mut_counter = 0
wt_counter = 0
C_counter = 0
G_counter = 0

# collect map of read_ids to base according to ref position
pileup_iter = bamfile.pileup(vcf_record.CHROM, vcf_record.POS - 1, vcf_record.POS,
                             truncate=True,
                             stepper="samtools",
                             min_base_quality=0,
                             ignore_overlaps=False,
                             ignore_orphans=False)

read_id_to_base_map = dict()

for column in pileup_iter:
    if(column.reference_pos == vcf_record.POS - 1):
        for read in column.pileups:
            if read.is_del:
                read_id_to_base_map[read.alignment.query_name] = 'del'
            elif read.is_refskip:
                read_id_to_base_map[read.alignment.query_name] = 'ins'
            else:
                read_id_to_base_map[read.alignment.query_name] = read.alignment.query_sequence[read.query_position]
bamfile.close()

# parse bamfile again
bamfile = pysam.AlignmentFile(bam_fn, 'rb')

# attach read id, classification ,and base into one tuple
for entry in bamfile:

    classification = ''
    call_at_locus = ''

    if entry.is_unmapped:
        classification = 'unmapped'
    else:

        if is_in_target_region(vcf_record, entry):

            if covers_target_locus(vcf_record, entry):
                classification = 'target'
            else:
                classification = 'fragment'

        else:
            classification = 'background'

    #print(entry.query_name, classification, read_id_to_base_map[entry.query_name])
    if entry.query_name not in read_id_to_base_map.keys():
        base = 'X'
    else:
        base = read_id_to_base_map[entry.query_name]
        
    if base == mut_base:
        mut_counter = mut_counter + 1
    elif base == wt_base:
        wt_counter = wt_counter + 1
    
    read_classification_list.append((entry.query_name, classification,base))
    read_classification_counters[classification] += 1
    
for name, classification,base in read_classification_list:
    print("{},{},{}".format(name, classification,base))

print(mut_counter, wt_counter)

# write classifications to file
with open('bam_classifications.csv', 'w') as filehandle:
    filehandle.writelines("{},{},{}\n".format(name, classification, base) for name,classification,base in read_classification_list)
        

# print summary of classifications
#print(read_classification_counters)
