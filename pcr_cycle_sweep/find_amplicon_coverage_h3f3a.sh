#!/bin/bash
BAMFILE=$1

# NOTE: MUST PRE FILTER FOR SUPPLEMENTARY AND SECONDARY ALIGNMENTS

samtools mpileup -Q 0 -q 0 -A -d 200000 ${BAMFILE} -o pileup.txt

# get max depth
max_depth="$(cut -d$'\t' -f 4 pileup.txt | sort -n | uniq | tail -n 1)"

python count_all_amplicons.py > amplicons.bed

total_alignments="$(samtools flagstat ${BAMFILE} | head -n 1 | cut -d' ' -f 1)"
secondary_alignments="$(samtools flagstat ${BAMFILE} | head -n 2 | tail -n 1 | cut -d' ' -f 1)"

supp_alignments="$(samtools flagstat ${BAMFILE} | head -n 3 | tail -n 1 | cut -d' ' -f 1)"
total_reads="$((total_alignments - secondary_alignments - supp_alignments))"
unmapped_reads="$(samtools view -f4 -c ${BAMFILE})"

echo "Total Reads:"
echo $total_reads
echo "Mapped Reads:"
echo "$((total_reads - unmapped_reads))"
echo "Total Amplicons:"
samtools view -q 0 -L amplicons.bed ${BAMFILE} | wc -l
echo "Target Coverage:"
samtools view -q 0 -L h3f3a_coverage.bed ${BAMFILE} | wc -l
echo "Max Target Depth:"
echo $max_depth
echo "SNP Depth:"
samtools view -q 0 -L h3f3a_depth.bed ${BAMFILE} | wc -l
