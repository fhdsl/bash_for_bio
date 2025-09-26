#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=1-3
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=00:10:00
file_array=(../../data/*.fastq)
ind=$((SLURM_ARRAY_TASK_ID-1))
echo "$ind"
current_file=${file_array[$ind]}
echo "$current_file"
./run_bwa.sh $current_file
