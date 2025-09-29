#!/bin/bash                           
module load SAMtools/1.19.2-GCC-13.2.0  
samtools view -c $1 > $1.counts.txt     
module purge                            