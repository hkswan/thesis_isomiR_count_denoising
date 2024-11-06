rm(list = ls())

#load necessary R packages
library(Biostrings)
library(DESeq2)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
miRNA_idx = as.numeric(args)[1]

#load functions
source("/scratch/hswan/thesis_isomiR_count_denoising/correct_isomiR_counts_v2_FUNCTION.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/correct_technical_length_variant_functions.R")


#load ernesto's benchmark dataset and get unique miRNAs out of it:
load("/scratch/mmccall2_lab/miRNA_reference/Rdata/summarized_experiment_isomiR.RData")
isomiR_se_object = se_object
rm(se_object)

#make rowdata
rowdata = rowData(isomiR_se_object) %>% data.frame()
#make sure rowdata has appropriate column names that correct_isomiR_counts_v2 function is expecting 
colnames(rowdata) = c("uniqueSequence", "miRNA_name")

#make countdata
countdata = isomiR_se_object@assays@data$counts %>% data.frame()
colnames(countdata) = str_remove_all(colnames(countdata), "X")

exp_miRNAs = unique(rowdata$miRNA_name)

my_miRNA = exp_miRNAs[miRNA_idx]

cat("Removing technical isomiRs from sample 1L for miRNA", my_miRNA, "\n")

file_name = paste0(my_miRNA, "_rmv_all_isomiRs_sample_1L.Rds")
file_path = paste0("/scratch/hswan/thesis_isomiR_count_denoising/isomiR_partitioning_results/")

cat("Partitioning results will be saved to", paste0(file_path, file_name, collapse=""), "if it doesn't already exist.\n")

if(!(file_name %in% list.files(file_path))){
  tst_remove_all = correct_isomiR_counts_v2(my_miRNA, rowdata=rowdata, countdata=countdata, sample_idx=2, max_iterations = 20, OMEGA_A=0.05)
  saveRDS(tst_remove_all, paste0(file_path, file_name, collapse=""))
} else{
  cat(file_name, "already exists in file path. Do you want to load it instead?\n")
}


