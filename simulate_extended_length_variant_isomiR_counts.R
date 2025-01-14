
rm(list = ls())

args = commandArgs(trailingOnly = TRUE)
mirna_idx = as.numeric(args[1])
num_datasets = as.numeric(args[2])

source("/scratch/hswan/thesis_isomiR_count_denoising/simulate_length_variant_isomiRs_EXTENDED_FUNCTION.R")
transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

mousedata = load_mouse_miRNA_data()
rowdata = mousedata$rowdata
countdata = mousedata$countdata

count_df = data.frame(rowdata, count = rowSums(countdata))

true_miRNAs = get_true_mouse_miRNAs(rowdata)

miRNA = true_miRNAs[mirna_idx]
cat("Generating simulated length variant isomiRs for miRNA", miRNA, "\n")

datasets = list()
for(i in 1:num_datasets){
  datasets[[i]] = simulate_length_variant_isomiRs_extened(rowdata, count_df, miRNA, i, 1000, 100, transition_probs)
}

save_file_path = "/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/extended_length_variants_NO_NOISE/"
save_file_name = paste0(miRNA, "_", num_datasets, "_simulated_length_variant_isomiRs_datasets.Rds")
cat("Saving datasets to", paste0(save_file_path, save_file_name, collapse=""), "\n")

saveRDS(datasets, paste0(save_file_path, save_file_name, collapse=""))