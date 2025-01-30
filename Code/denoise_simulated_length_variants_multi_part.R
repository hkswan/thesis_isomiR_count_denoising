source("/scratch/hswan/thesis_isomiR_count_denoising/Code/denoise_isomiR_counts_update_partitions_at_end.R")
transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

data_file_path = "/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/length_variants_multi_partition/01_24_2025/"
multi_partition_data_files = list.files(data_file_path)

#load datafile 

data_set_idx = commandArgs(trailingOnly=TRUE) %>% as.numeric()

data_list = readRDS(paste0(data_file_path, multi_partition_data_files[data_set_idx]))

#always the same 3 sequences 
#i think they're so closely related that it's possible for some of the isomiRs we made from the 1st center sequence to also be from
#the 2nd center sequence and so on - so after 1st iteration they don't have a significant p-value so they aren't moving 


miRNA = strsplit(multi_partition_data_files[[data_set_idx]], "_simulated")[[1]][1]

start = Sys.time()
denoise_objs = list()
for(i in 1:length(data_list)){
  dataset = data_list[[i]]
  cat("Denoising simulated data from dataset", i, "\n")
  denoise_objs[[i]] = denoise_isomiR_counts(dataset$sim_rowdata, dataset$counts, transition_probs, dataset$sim_rowdata$miRNA_name[1],
                                            0.05, 10, "BH")
}
end = Sys.time()
elapsed = end - start 
cat("Time it took to denoise all 50 datasets:\n")
print(elapsed)
cat("\n")

calculate_results = function(denoised_obj, true_center_seqs){
  inferred_center_seqs = filter(denoised_obj, center == 1) %>% select(., uniqueSequence) %>% unlist()
  inferred_tech_isomiRs = filter(denoised_obj, center == 0) %>% select(., uniqueSequence) %>% unlist()
  num_true_positives = sum(true_center_seqs %in% inferred_center_seqs)
  num_false_positives = sum(!(inferred_center_seqs %in% true_center_seqs))
  num_true_negatives = sum(!(inferred_tech_isomiRs %in% true_center_seqs))
  num_false_negatives = length(inferred_tech_isomiRs) - num_true_negatives
  
  fp_rate = num_false_positives / (num_false_positives + num_true_negatives)
  tp_rate = num_true_positives / (num_true_positives + num_false_negatives)
  fn_rate = 1 - tp_rate
  tn_rate = 1 - fp_rate 
  
  return(list(num_true_positives=num_true_positives, num_false_positives=num_false_positives, num_true_negatives=num_true_negatives,
              num_false_negatives=num_false_negatives, fp_rate=fp_rate, tp_rate=tp_rate, fn_rate=fn_rate, tn_rate=tn_rate))

}

generate_true_partitions = function(data_list_isomiR_seqs, partition_uniqueSequences){
  num_partitions = length(data_list_isomiR_seqs)
  true_partitions = rep(0, length(partition_uniqueSequences))
  names(true_partitions) = partition_uniqueSequences
  for(j in 1:num_partitions){
    true_partitions[names(true_partitions) %in% data_list_isomiR_seqs[[j]] | names(true_partitions) == names(data_list_isomiR_seqs[j])] = j
  }
  
  return(true_partitions)
}

calculate_correctly_labeled_seqs = function(true_partition, inferred_partition){
  num_correct = sum(true_partition == inferred_partition)
  return(num_correct)
}

results_list = list()
for(i in 1:length(data_list)){
  results_list[[i]] = calculate_results(denoise_objs[[i]]$partition_df, data_list[[i]]$center_seqs)
}

save_file_path = "/scratch/hswan/thesis_isomiR_count_denoising/sims/01_24_2025/multi_partitions_MOVE_ALL/"

partition_path = paste0(save_file_path, "partition_objs/", collapse="")
results_path = paste0(save_file_path, "results/", collapse="")

partition_file = paste0(miRNA, "_partition_objs.Rds")
results_file = paste0(miRNA, "_results.Rds")

saveRDS(results_list, paste0(results_path, results_file))
saveRDS(denoise_objs, paste0(partition_path, partition_file))

# fp_rates = lapply(results_list, function(x) return(x[['fp_rate']])) %>% unlist() %>% round(., digits=3)
# table(fp_rates)
# 
# #power = 1-false negative rate 
# 
# fn_rates = lapply(results_list, function(x) return(x[["fn_rate"]])) %>% unlist() %>% round(., digits=3)
# power = 1 - fn_rates
# table(power)

# fp_df = data.frame(x = 1:length(fp_rates), y = fp_rates)
# ggplot(fp_df, aes(x=x, y=y)) + geom_point() + geom_abline(slope = 0, intercept = 0.05, col='red') + 
#   scale_y_continuous(limits = c(0.0, 0.05)) + xlab("Dataset idx") + ylab("False positive rate") + 
#   ggtitle(paste0("Observed false positive rate vs. dataset idx miRNA ", denoise_objs[[1]]$partition_df$miRNA_name[1]))
# 
# power_df = data.frame(x= 1:length(fn_rates), y = power)
# ggplot(power_df, aes(x=x, y=y)) + geom_point() + geom_abline(slope = 0, intercept = 0.80, col='darkgreen') + 
#   scale_y_continuous(limits=c(0.75, 1.04)) + xlab("Dataset idx") + ylab("Power") + 
#   ggtitle(paste0("Observed power vs. dataset idx miRNA ", denoise_objs[[1]]$partition_df$miRNA_name[1]))


# data1 = data_list[[1]]
# center_seqs = data1$center_seqs
# 
# center = center_seqs[1]
# center2 = center_seqs[2]
# center3 = center_seqs[3]
# 
# p1 = denoise_objs[[1]]$partition_df
# 
# isomiR_seqs2 = data1$isomiR_seqs[[center2]]
# mislabeled_isomiRs = filter(p1, uniqueSequence %in% isomiR_seqs2 & partition != 2) %>% filter(., count != 0) %>% 
#   select(., uniqueSequence) %>% unlist()
# 
# isomiR_seqs3 = data1$isomiR_seqs[[center3]]
# mislabeled_isomiRs3 = filter(p1, uniqueSequence %in% isomiR_seqs3 & partition != 3) %>% filter(., count != 0) %>%
#   select(., uniqueSequence) %>% unlist()
# lambda_1s = vector()
# lambda_2s = vector()
# 
# lambda_1_3s = vector()
# lambda_3s = vector()
# for(i in 1:length(mislabeled_isomiRs)){
#   a1 = Biostrings::pairwiseAlignment(pattern = center, mislabeled_isomiRs[i])
#   a2 = Biostrings::pairwiseAlignment(pattern = center2, mislabeled_isomiRs[i])
#   
# 
#   t1 = get_transitions(a1)
#   t2 = get_transitions(a2)
# 
#   l1 = compute_lambda(t1, transition_probs)
#   lambda_1s = c(lambda_1s, l1)
#   l2 = compute_lambda(t2, transition_probs)
#   lambda_2s = c(lambda_2s, l2)

  # nj1 = filter(p1, center == 1 & partition == 1) %>% select(., count) %>% unlist()
  # nj2 = filter(p1, center == 1 & partition == 2) %>% select(., count) %>% unlist()
  # 
  # mu1 = l1*nj1
  # mu2 = l2*nj2
  # 
  # y = filter(p1, uniqueSequence == mislabeled_isomiRs[i]) %>% select(., count) %>% unlist()

  # ppois(y-1, mu1, lower.tail=F)/(1-dpois(0, mu1))
  # ppois(y-1, mu2, lower.tail=F)/(1-dpois(0, mu2))
}

# for(i in 1:length(mislabeled_isomiRs3)){
#   a1 = Biostrings::pairwiseAlignment(pattern = center, mislabeled_isomiRs3[i])
#   a3 = Biostrings::pairwiseAlignment(pattern = center3, mislabeled_isomiRs3[i])
#   
#   t1 = get_transitions(a1)
#   t3 = get_transitions(a3)
#   
#   l1 = compute_lambda(t1, transition_probs)
#   lambda_1_3s = c(lambda_1_3s, l1)
#   l3 = compute_lambda(t3, transition_probs)
#   lambda_3s = c(lambda_3s, l3)
# }
# lambda_1s
# lambda_2s
# 
# sum(lambda_2s > lambda_1s)
# 
# 
# lambda_1_3s
# lambda_3s
