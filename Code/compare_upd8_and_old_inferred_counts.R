rm(list = ls())

source("/scratch/hswan/thesis_isomiR_count_denoising/Code/denoise_isomiR_counts_WORKING_FUNCTION.R")

transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

partition_obj_path = "/scratch/hswan/thesis_isomiR_count_denoising/sims/01_24_2025/multi_partitions/partition_objs/"
data_path = "/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/length_variants_multi_partition/01_24_2025/"

partition_obj_files = list.files(partition_obj_path)
data_list_files = list.files(data_path)

args = commandArgs(trailingOnly = TRUE)
dataset_idx = args[1] %>% as.numeric()
miRNA_idx = args[2] %>% as.numeric()

partition_objs = readRDS(paste0(partition_obj_path, partition_obj_files[miRNA_idx]))
miRNA = strsplit(partition_obj_files[miRNA_idx], "_partition")[[1]][1]

data_file = data_list_files[grepl(miRNA, data_list_files)]
data_list = readRDS(paste0(data_path, data_file, collapse=""))



partitioning = partition_objs[[dataset_idx]]$partition_df

center_seqs = data_list[[dataset_idx]]$center_seqs
isomiR_seqs = data_list[[dataset_idx]]$isomiR_seqs

#create true_partition column
partitioning$true_partition = rep(0, nrow(partitioning))
#get isomiR sequences from associated data object:
true_isomiR_seqs = data_list[[dataset_idx]]$isomiR_seqs
center_seqs = names(true_isomiR_seqs)
for(j in center_seqs){
  p = filter(partitioning, uniqueSequence == j) %>% select(., partition) %>% unlist() %>% unname()
  partitioning$true_partition[partitioning$uniqueSequence %in% true_isomiR_seqs[[j]]] = p
}
master_alignments = list()
master_transitions = list()
master_lambdas = list()
inferred_isomiR_seqs = filter(partitioning, center == 0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
for(j in center_seqs){
  master_alignments[[j]] = lapply(inferred_isomiR_seqs, Biostrings::pairwiseAlignment, pattern = j)
  names(master_alignments[[j]]) = inferred_isomiR_seqs
}
for(j in names(master_alignments)){
  master_transitions[[j]] = lapply(master_alignments[[j]], get_transitions)
}
for(j in names(master_transitions)){
  master_lambdas[[j]] = lapply(master_transitions[[j]], compute_lambda, transition_probs = transition_probs)
}

partitioning$updated_partition = rep(0, nrow(partitioning))
#get updated partition memberships 
for(i in inferred_isomiR_seqs){
  l = lapply(master_lambdas, function(x) return(x[[i]])) %>% unlist()
  upd8_p = which(l == max(l))
  partitioning$updated_partition[partitioning$uniqueSequence==i] = upd8_p
}

for(i in 1:3){
  partitioning$updated_partition[partitioning$uniqueSequence == center_seqs[i]] = i
  partitioning$true_partition[partitioning$uniqueSequence==center_seqs[i]] = i
}

true_counts = sapply(1:3, function(x) filter(partitioning, true_partition==x) %>% select(., count) %>% unlist() %>% sum())
upd8_counts = sapply(1:3, function(x) filter(partitioning, updated_partition==x) %>% select(., count) %>% unlist() %>% sum())
old_infer_counts = sapply(1:3, function(x) filter(partitioning, partition ==x) %>% select(., count) %>% unlist() %>% sum())

infer_counts_mat = rbind(true_counts, upd8_counts, old_infer_counts)

save_file_path = "/gpfs/fs2/scratch/hswan/thesis_isomiR_count_denoising/sims/01_24_2025/multi_partitions/upd8_count_estimation/"
if(!(miRNA) %in% list.files(save_file_path)){
  #setwd(save_file_path)
  dir.create(file.path(paste0(save_file_path, miRNA, collapse="")))
}

filename = paste0("infer_counts_mat_dataset_", dataset_idx, ".Rds", collapse="")
saveRDS(infer_counts_mat, paste0(save_file_path, miRNA, "/", filename))
