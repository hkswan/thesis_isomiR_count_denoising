#R script that analyzes simulated data created in the generate_simulated_data.R script and exports 
#resulting partition_df, power, and type I error rate to sims/results folder


source("/scratch/hswan/thesis_isomiR_count_denoising/denoise_isomiR_counts_WORKING_FUNCTION.R")

args = commandArgs(trailingOnly = TRUE)

dataset_idx = as.numeric(args[1])
file_idx = as.numeric(args[2])

transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

dataset_files = list.files("/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/updated_simulated_data/")
filepath = "/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/updated_simulated_data/"
datasets = readRDS(paste0(filepath, dataset_files[file_idx], collapse=""))

data = datasets[[dataset_idx]]

simulated_rowdata=data$simulated_rowdata
simulated_count_df=data$simulated_counts 

miRNA = simulated_rowdata$miRNA_name[1]
obj = denoise_isomiR_counts(simulated_rowdata, simulated_count_df, transition_probs, miRNA, 0.05, 10, "BH")

partition_df = obj$partition_df

save_file_name =  paste0(miRNA, "_", length(datasets)-1, "_simulated_datasets_partition_obj_", dataset_idx,
".Rds", collapse="")

save_file_path = paste0("/scratch/hswan/thesis_isomiR_count_denoising/sims/partition_objs/", miRNA, "/", 
                        collapse="")

if(save_file_name %in% list.files(save_file_path)){
  #load 
  partitions = readRDS(paste0(save_file_path, save_file_name, collapse=""))
  #add to list
  partitions[[length(partitions)+1]] = partition_df
  #save again
  saveRDS(partitions, paste0(save_file_path, save_file_name, collapse=""))
} else{
  dir.create(save_file_path)
  partitions = list(partition_df)
  saveRDS(partitions, paste0(save_file_path, save_file_name, collapse=""))
  
}
