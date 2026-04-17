#!/bin/bash
#SBATCH --array=1-3
#SBATCH --nodes=1
sleep 120
echo "${SLURM_ARRAY_TASK_ID} job"