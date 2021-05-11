# UltraRapidSeq
This repository is meant to accompany the paper "Ultra-rapid Sequencing-based Molecular Diagnostics via Real-Time Threshold Sequencing." It is a collection of scripts and data used to perform certain analysis from the paper, and shared for posterity and reproducibility.

## VAF Over Time
Usage: ```python3 <guppy_generated_fastq> <bamfile> <VCF record>```

Will sort fastq file based on Guppy timestamps and compute VAF support over time based on supplied VCF record.

## PCR Cycle Sweep
Scripts to classify alignments into background reads, amplicon fragments, and on target reads.

Usage: ```find_amplicon_coverage_hist1h3b.sh <bamfile>```
Usage: ```find_amplicon_coverage_h3f3a.sh <bamfile>```

### bam_classifier.py 
Re-implementation in pure python that should produce identical results to find_amplicon_coverage_xxx.sh.
Usage: ```python3 bam_classifier.py <bamfile> <vcffile>```

## Calc Sampling Rate
Scripts to compute MinION Sampling Rate per channel, and the number of active/participating pores over time.

Usage: ```python3 calc_sampling_rate.py <guppy_generated_fastq> <analysis start time (s)> <analysis finish time (s)>```

## Common

### basecall_fast5.sh
Simple guppy wrapper. Script was used to basecall SQK-RAD004 runs. These are the end-to-end experiments.

### basecall_barcoded_fast5.sh
Simple guppy_wrapper. Script used to basecall SQK-RBK004 runs. These are the time-sweep/cycle-sweep experiments.

### align_fastq.sh
Simple minimap2 wrapper.

### sam_to_bam.sh
SAM to BAM conversion and secondary/supplementary alignment filtering.
