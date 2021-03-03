#!/bin/bash

guppy_basecaller \
    -i fast5 \
    -s fastq \
    --kit SQK-RAD004 \
    --flowcell FLO-MIN106 \
    --num_callers 1 \
    -x "cuda:0"
