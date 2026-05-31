#!/bin/bash

module load samtools/1.22.1

samtools view -Sb sample.sam -o sample.bam;
samtools sort sample.bam > sample.sorted.bam;
samtools index sample.sorted.bam;
