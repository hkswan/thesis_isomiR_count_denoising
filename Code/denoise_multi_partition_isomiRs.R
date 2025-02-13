
source("/scratch/hswan/thesis_isomiR_count_denoising/Code/denoise_isomiR_counts_update_partitions_at_end.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/Code/compute_performance_metrics_FUNCTION.R")
transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

#set working directory just in case we forgot
setwd("/scratch/hswan/thesis_isomiR_count_denoising")

file_path = "./data/simulated_data/multi_partition_isomiRs/"

simulated_isomiRs = list.files("./data/simulated_data/multi_partition_isomiRs")

#load a datalist:
idx = commandArgs(trailingOnly = TRUE) %>% as.numeric()
data_list = readRDS(paste0(file_path, simulated_isomiRs[idx], collapse=""))

#pick out first dataset:
rowdata = data_list[[1]]$sim_rowdata
countdata = data_list[[1]]$sim_counts


cat("Rowdata dimensions:", dim(rowdata), "\n")
cat("countdata dimensions:", length(countdata), "\n")

head(rowdata)
head(countdata)
table(countdata)

miRNA = strsplit(simulated_isomiRs[[idx]], "_25_simulated")[[1]][1]
cat("miRNA:", miRNA, "\n")

denoising_list = list()

start=Sys.time()
for(i in 1:length(data_list)){
  cat(i, "\n")
  rowdata = data_list[[i]]$sim_rowdata
  countdata = data_list[[i]]$sim_counts
  denoising_list[[i]] = denoise_isomiR_counts(rowdata, countdata, transition_probs, miRNA, 0.05, 
                                              10, "BH")
}
end=Sys.time() - start
cat("\n Time to denoise all 25 simulated datasets:\n")
print(end)

partitions = lapply(denoising_list, function(x) return(x[['partition_df']]))

lapply(partitions, function(x) return(max(x$partition))) %>% unlist() %>% table()

performance_metrics_list = list()
for(i in 1:length(data_list)){
  performance_metrics_list[[i]]=compute_performance_metrics(partitions[[i]], data_list[[i]]$center_seqs)
}

save_path = "./sims/02_03_2025/"

results_file = paste0(save_path, "results/", miRNA, "_results.Rds")
partition_file = paste0(save_path, "partition_objs/", miRNA, "partition_objs.Rds")

cat("Saving results to", results_file, "\n")
saveRDS(performance_metrics_list, results_file)
cat("Saving partitions to", partition_file, "\n")
saveRDS(denoising_list, partition_file)


#fps = lapply(performance_metrics_list, function(x) return(x[['fp_rate']])) %>% unlist()
#tps  = lapply(performance_metrics_list, function(x) return(x[['tp_rate']])) %>% unlist()
