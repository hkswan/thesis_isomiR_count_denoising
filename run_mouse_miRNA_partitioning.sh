#!/bin/bash
#SBATCH --partition=debug --time=01:00:00  --output=pairwise_align.log
#SBATCH --cpus-per-task=2
##SBATCH --mem-per-cpu 16GB
#SBATCH -o ./logs/log_%j.txt 
#SBATCH -e ./errors/err_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=hannah_swan@URMC.Rochester.edu 

#SBATCH -a 1-5

module load r/4.2.1
##echo This is miRNA $i
##echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_I."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)
Rscript /scratch/hswan/thesis_isomiR_count_denoising/mouse_isomiR_partitioning.R $SLURM_ARRAY_TASK_ID

