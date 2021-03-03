#!/bin/bash

guppy_basecaller \
    -i fast5 \
    -s fastq_barcoded \
    --trim_barcodes \
    --kit SQK-RBK004 \
    --barcode_kits "SQK-RBK004" \
    --flowcell FLO-MIN106 \
    --num_callers 4 \
    -x "cuda:0"
