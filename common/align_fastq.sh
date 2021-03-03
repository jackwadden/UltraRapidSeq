#!/bin/bash

REF_IDX=${1}
FASTQ=${2}
KEY=${FASTQ%.fastq}
echo ${KEY}

minimap2 -a \
         -x map-ont \
	 -t 1 \
	 --eqx \
	 -R '@RG\tID:foo\tSM:bar' \
	 -o ${KEY}.sam \
	 ${REF_IDX} \
	 ${FASTQ}
