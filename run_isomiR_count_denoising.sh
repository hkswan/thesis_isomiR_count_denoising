#!/bin/bash
#SBATCH --partition=standard --time=01:00:00  
#SBATCH --cpus-per-task=2
##SBATCH --mem-per-cpu 16GB
#SBATCH -o ./logs/denoise_isomiR_counts_log.txt 
#SBATCH -e ./errors/denoise_isomiR_counts_err.txt
#SBATCH --mail-type=all
#SBATCH --mail-user=hannah_swan@URMC.Rochester.edu 



module load r/4.2.1
Rscript /scratch/hswan/thesis_isomiR_count_denoising/isomiR_count_denoising.R 

