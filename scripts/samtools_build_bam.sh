#!/bin/bash

CONTIGS=$1
OUTPUT_FILE=$2
N_THREADS=$3

samtools faidx $CONTIGS
samtools view -bt $CONTIGS.fai $OUTPUT_FILE.sam > $OUTPUT_FILE.bam
samtools sort -o $OUTPUT_FILE.bam $OUTPUT_FILE.bam
samtools index $OUTPUT_FILE.bam