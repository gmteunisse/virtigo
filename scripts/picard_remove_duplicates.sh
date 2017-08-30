#!/bin/bash

BIN_DIR=$1
BAM_FILE=$2
OUTPUT_FILE=$3

java -jar $BIN_DIR/picard.jar MarkDuplicates \
    INPUT=$BAM_FILE.bam \
    OUTPUT=$BAM_FILE.bam \
    METRICS_FILE=$OUTPUT_FILE.metrics \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=FALSE

samtools sort -o $BAM_FILE.bam $BAM_FILE.bam
samtools index $BAM_FILE.bam
samtools flagstat $BAM_FILE.bam > $OUTPUT_FILE.flagstat