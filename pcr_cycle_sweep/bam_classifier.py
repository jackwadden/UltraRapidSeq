import sys, os, argparse, collections

import pysam
import vcfpy


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

print(vcf_chrom, vcf_locus)

# for each entry, categorize as "target, fragment, background, unmapped"
read_classification_list = list()
read_classification_counters = {'target' : 0, 'fragment': 0, 'background' : 0, 'unmapped' : 0}

for entry in bamfile:
    #print(entry.query_name, entry.reference_name)

    classification = ''

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

    read_classification_list.append((entry.query_name, classification))
    read_classification_counters[classification] += 1
    
#for name, classification in read_classification_list:
#    print("{},{}".format(name, classification))

# write classifications to file
with open('bam_classifications.csv', 'w') as filehandle:
    filehandle.writelines("{},{}".format(name, classification) for name,classification in read_classification_list)
        

# print summary of classifications
print(read_classification_counters)
