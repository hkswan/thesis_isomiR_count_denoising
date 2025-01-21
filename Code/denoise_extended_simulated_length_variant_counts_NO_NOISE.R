rm(list = ls())

source("/scratch/hswan/thesis_isomiR_count_denoising/denoise_isomiR_counts_WORKING_FUNCTION.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/load_mouse_miRNA_data_function.R")

transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

length_variant_data_path = "/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/extended_length_variants_NO_NOISE/01_14_2025/"
length_variant_data_files = list.files(length_variant_data_path)

args = commandArgs(trailingOnly = TRUE)
data_file_idx = as.numeric(args[1])

cat("Loading datasets from", length_variant_data_files[data_file_idx], "\n")

datasets = paste0(length_variant_data_path, length_variant_data_files[data_file_idx], collapse="") %>% readRDS()

partition_objs = list()

for(i in 1:length(datasets)){
  rowdata_sim = datasets[[i]]$sim_rowdata
  countdata_sim = datasets[[i]]$sim_counts 
  cat("Denoising isomiR counts from dataset", i, "\n")
  partition_objs[[i]] = denoise_isomiR_counts(rowdata_sim, countdata_sim, transition_probs, rowdata_sim$miRNA_name[1], 0.05, 10, "BH")
}

num_partitions_created = lapply(partition_objs, function(x) return(x[['partition_df']]$partition %>% max())) %>% unlist()
num_seqs = lapply(partition_objs, function(x) return(x[['partition_df']] %>% nrow())) %>% unlist()
num_false_positives = num_partitions_created-1
num_true_negatives = lapply(partition_objs, function(x) return(x[['partition_df']] %>% filter(., center == 0 & partition ==1) %>% nrow())) %>% unlist()
num_true_positives = lapply(partition_objs, function(x) return(x[['partition_df']] %>% filter(., center == 1 & partition == 1) %>% nrow())) %>%
  unlist()
num_false_negatives = num_seqs-num_false_positives-num_true_negatives-num_true_positives

miRNA = strsplit(length_variant_data_files[data_file_idx], "_100")[[1]][1]

save_file_path = "/scratch/hswan/thesis_isomiR_count_denoising/sims/01_14_2025/extended_length_variants_NO_NOISE/"

save_file_name_partition = paste0(save_file_path, "partition_objs/", miRNA, "_simulated_length_variant_partitions.Rds", collapse="")
save_file_name_results = paste0(save_file_path, "results/", miRNA, "_simulated_length_variant_results.Rds", collapse="")

results = list(num_partitions_created=num_partitions_created, num_seqs=num_seqs, num_false_positives=num_false_positives, 
               num_true_negatives=num_true_negatives, num_true_positives=num_true_positives, num_false_negatives=num_false_negatives)

cat("Saving partition_objs object\n")
saveRDS(partition_objs, save_file_name_partition)
cat("Saving results")
saveRDS(results, save_file_name_results)
