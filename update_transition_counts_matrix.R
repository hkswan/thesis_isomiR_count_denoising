rm(list=ls())

library(tidyverse)
library(Biostrings)
source("/scratch/hswan/thesis_isomiR_count_denoising/correct_technical_length_variant_functions.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/function_dev/get_miRNA_level_transition_counts_FUNCTION.R")

f1 = list.files("/scratch/hswan/thesis_isomiR_count_denoising/rmv_all_err_seqs")
f1 = f1[grepl(".Rds", f1)]


args = commandArgs(trailingOnly = TRUE)
idx = as.numeric(args[1])

cat("File name:", f1[idx], "\n")

filepath = paste0("/scratch/hswan/thesis_isomiR_count_denoising/rmv_all_err_seqs/", 
                  f1[idx], collapse="")

partition_df =  readRDS(filepath)$partition_df

trns_counts = get_miRNA_level_transition_counts(partition_df)


miRNA = partition_df$miRNA_name[1]


save_file_name = paste0(miRNA, "_update_trns_counts.Rds")
save_file_path=paste0("/scratch/hswan/thesis_isomiR_count_denoising/transition_count_matrices/update/", save_file_name)
cat("saving to", save_file_path, "\n")
saveRDS(trns_counts, save_file_path)
