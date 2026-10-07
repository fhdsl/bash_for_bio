file_array=(../../data/*.fastq)   
num_files=${#file_array[@]}
max_index=$((num_files - 1))
sbatch --array=0-${max_index} run_sbatch2.sh
