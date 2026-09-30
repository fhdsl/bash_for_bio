#!/bin/bash                           
module load SAMtools/1.19.2-GCC-13.2.0  
samtools view -c CALU1_combined_final.sam > CALU1_combined_final.sam.counts.txt
module purge                            
