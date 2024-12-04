
source("/scratch/hswan/thesis_isomiR_count_denoising/denoise_isomiR_counts_WORKING_FUNCTION.R")
args = commandArgs(trailingOnly=TRUE)
idx = as.numeric(args[1])
cat("miRNA_idx:", idx, "\n")

#load mouse data
mouse_data = load_mouse_miRNA_data()

#get rowdata, countdata objects
rowdata = mouse_data$rowdata
countdata = mouse_data$countdata

#get rowsums to create count_df object (complete pooling of samples)
count_df = rowSums(countdata)

#get unique miRNAs, id miRNA whose isomiR counts we are denoising
unique_miRNAs = unique(rowdata$miRNA_name)
x = unique_miRNAs[idx]

#load transition probability matrix
transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

obj = denoise_isomiR_counts(rowdata, count_df, transition_probs, x, 0.05, 100, "BH")

filename = paste0(x, "_test_run_isomiR_count_denoising.Rds", collapse="")
filepath = "/scratch/hswan/thesis_isomiR_count_denoising/test_run_isomiR_count_denoising/"
saveRDS(obj, paste0(filepath, filename, collapse=""))