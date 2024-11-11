#!/bin/bash
#SBATCH --partition=debug  --time=01:00:00  
#SBATCH -o ./logs/log_debug_correct_isomiR_%j.txt
#SBATCH -e ./errors/err_debug_correct_isomiR_%j.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=hannah_swan@URMC.Rochester.edu 

module load r/4.2.1
Rscript /scratch/hswan/thesis_isomiR_count_denoising/debug_correct_isomiR_count_step.R

