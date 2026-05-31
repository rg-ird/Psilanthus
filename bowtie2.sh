#!/bin/bash

module load bowtie2/2.3.4.1 

bowtie2 --no-unal --very-sensitive -q -x reference.fa -U sample.fq.gz -S sample.sam
