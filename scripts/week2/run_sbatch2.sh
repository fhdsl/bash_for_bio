#!/bin/bash

file_array=(../../data/*.fastq)   
ind=$((SLURM_ARRAY_TASK_ID-1))
echo "$ind"
current_file=${file_array[$ind]}
echo "$current_file"
# execute the `run_bwa.sh` script on $current_file
./run_bwa.sh $current_file
