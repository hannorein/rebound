#!/bin/bash
#SBATCH --array=0-239 
#SBATCH --partition=stars
#SBATCH --partition=stars
#SBATCH --job-name=1jac
#SBATCH --mem=10M    
#SBATCH --output=slurm.log
#
echo `date`": starting taskid="$SLURM_ARRAY_TASK_ID "on host" `hostname`
./rebound $SLURM_ARRAY_TASK_ID 
echo `date`": finished taskid="$SLURM_ARRAY_TASK_ID
