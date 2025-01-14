#write a script that denoises simulated length variant isomiR counts and returns the number of false positives, true negatives, etc.
#saves partition objects output and results to separate .Rds files 

rm(list = ls())

source("/scratch/hswan/thesis_isomiR_count_denoising/denoise_isomiR_counts_WORKING_FUNCTION.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/load_mouse_miRNA_data_function.R")

library(tidyverse)
library(Biostrings)

transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

sim_data_filepath = "/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/length_variants/"
sim_data_files = list.files(sim_data_filepath)

args = commandArgs(trailingOnly=TRUE)
datafile_idx = as.numeric(args[1])

cat("Loading datasets from file", sim_data_files[datafile_idx], "\n")

datasets = readRDS(paste0(sim_data_filepath, sim_data_files[datafile_idx], collapse=""))

partition_objs = list()
for(i in 1:length(datasets)){
  rowdata_sim = datasets[[i]]$rowdata
  countdf_sim = datasets[[i]]$counts
  partition_objs[[i]] = denoise_isomiR_counts(rowdata_sim, countdf_sim, transition_probs, rowdata_sim$miRNA_name[1], 0.05, 10, "BH")
}

cat("Calculating false positive rate \n")

false_positive_counts = vector(length = length(partition_objs))
true_negative_counts = vector(length = length(false_positive_counts))
true_positive_counts = vector(length = length(false_positive_counts))

for(i in 1:length(false_positive_counts)){
  df = partition_objs[[i]]$partition_df
  false_positive_counts[i] = filter(df, center == 1) %>% filter(., uniqueSequence != datasets[[i]]$center_seq) %>% nrow()
  true_negative_counts[i] = filter(df, center == 0) %>% filter(., uniqueSequence != datasets[[i]]$center_seq) %>% nrow()
  true_positive_counts[i] = filter(df, center == 1) %>% filter(., uniqueSequence == datasets[[i]]$center_seq) %>% nrow()
    
}

negative_counts = false_positive_counts+true_negative_counts

false_positive_rates = false_positive_counts/negative_counts 

results = list(false_positive_rates=false_positive_rates, avg_fp_rate=mean(false_positive_rates), false_positive_counts=false_positive_counts,
               true_negative_counts=true_negative_counts, true_positive_counts=true_positive_counts)


save_path = '/scratch/hswan/thesis_isomiR_count_denoising/sims/01_14_2025/'

results_save_file = paste0(miRNA,"_simulated_length_variant_results.Rds")
cat("Saving results to", paste0(save_path, "results/", results_save_file, collapse=""), "\n")
saveRDS(results, paste0(save_path, "results/", results_save_file, collapse=""))

partition_objs_save_file = paste0(miRNA, "_simulated_length_variant_partition_objs.Rds")
cat("Saving results to", paste0(save_path, "partition_objs/", partition_objs_save_file, collapse=""), "\n")
saveRDS(partition_objs, paste0(save_path, "partition_objs/", partition_objs_save_file, collapse=""))
