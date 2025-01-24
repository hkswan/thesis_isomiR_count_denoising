
#load functions that we've written
source("/scratch/hswan/thesis_isomiR_count_denoising/Code/simulate_length_variant_isomiRs_multi_partition_FUNCTION.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/Code/denoise_isomiR_counts_WORKING_FUNCTION.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/Code/load_mouse_miRNA_data_function.R")

#load mirna data
mousedata = load_mouse_miRNA_data()
rowdata = mousedata$rowdata
countdata = mousedata$countdata %>% rowSums()

#load transition probabilities matrix
transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

true_miRNAs = get_true_mouse_miRNAs(rowdata)

args = commandArgs(trailingOnly = TRUE)
miRNA_idx = as.numeric(args[1])
num_datasets = as.numeric(args[2])
num_true_partitions = as.numeric(args[3])
date_run = args[4]

miRNA = true_miRNAs[miRNA_idx]
cat("miRNA:", miRNA, "\n")
cat("Simulating", num_datasets,  "datasets each with", num_true_partitions,  "true isomiR sequences \n")

datasets_list = list()
for(n in 1:num_datasets){
  cat("Now generating dataset number", n, "\n")
  datasets_list[[n]] = simulate_length_variant_isomiRs(rowdata,  countdata, miRNA,  FALSE, num_true_partitions, 100, 1000, 
                                  transition_probs, n)
}

dataset_save_filepath = paste0("/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/length_variants_multi_partition/", 
                               date_run, "/", collapse="")
filename = paste0(miRNA,"_simulated_length_variants_", num_datasets, "_datasets_", num_true_partitions, 
                  "_true_partitions.Rds")
cat("Saving generated datasets to", paste0(dataset_save_filepath, filename, collapse=""), "\n")

saveRDS(datasets_list, paste0(dataset_save_filepath, filename))
