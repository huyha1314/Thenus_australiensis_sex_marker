#!/bin/bash
source config.sh

echo "Starting 1. Quailty control of raw reads with fastp (fastp 0.23.4)"
pixi run fastp \
        -i "$i1" -I "$i2" \
        -o result/fastp/trim.${sample}_1.fq.gz \
        -O result/fastp/trim.${sample}_2.fq.gz \
        --trim_front1 8 --trim_front2 8 \
        --length_required 50 \
        --qualified_quality_phred 25 \
        --thread 16 \
        --html result/multiqc/${sample}.report.html \
        --json result/multiqc/${sample}.report.json" >> fastp.txt
