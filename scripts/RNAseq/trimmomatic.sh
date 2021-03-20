#!/bin/sh

sample=$1
describer=$(echo ${sample} | sed 's/_R1.fastq.gz//')

java -jar trimmomatic-0.35.jar PE ${describer}_R1.fastq.gz ${describer}_R2.fastq.gz ${describer}_forward_paired.fq.gz \
${describer}_forward_unpaired.fq.gz ${describer}_reverse_paired.fq.gz ${describer}_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

seqtk mergepe ${describer}_forward_paired.fq.gz ${describer}_reverse_paired.fq.gz > ${describer}.pe.fq.gz
cat ${describer}_forward_unpaired.fq.gz ${describer}_reverse_unpaired.fq.gz > ${describer}.se.fq.gz