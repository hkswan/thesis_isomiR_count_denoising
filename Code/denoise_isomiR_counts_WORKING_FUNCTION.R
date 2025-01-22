

library(Biostrings)
library(tidyverse)
library(DESeq2)

source("/scratch/hswan/thesis_isomiR_count_denoising/Code/load_mouse_miRNA_data_function.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/Code/correct_technical_length_variant_functions.R")


denoise_isomiR_counts = function(rowdata, count_df, transition_probability_matrix,
                                 miRNA, omega_A=NULL, max_iter = NULL, adjust_method=c("BH", "Bonferroni")){
  #initialize partition_df object 
  cat("Creating initial partition_df object for miRNA", miRNA, "\n")
  partition_df = cbind(rowdata, count=count_df) %>% data.frame() %>% filter(., miRNA_name == miRNA)
  partition_df$partition = rep(1, nrow(partition_df))
  partition_df$center = rep(0, nrow(partition_df))
  
  #id initial center sequence 
  initial_center_seq = filter(partition_df, count==max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()
  
  partition_df$center[partition_df$uniqueSequence == initial_center_seq] = 1
  
  #with initial partition, complete 1st iteration outside of a while loop
  
  #STEP1 - get pairwise alignments between center sequence and all other sequences mapping to user-specified miRNA
  isomiR_seqs = partition_df$uniqueSequence[partition_df$center==0]
  cat("Getting initial alignments\n")
  alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern=initial_center_seq)
  names(alignments) = isomiR_seqs
  cat("Getting initial transitions \n")
  transitions = lapply(alignments, get_transitions)
  cat("Calculating initial lambdas\n")
  lambdas = lapply(transitions, compute_lambda, transition_probs=transition_probability_matrix)
  master_lambdas = list(lambdas)
  names(master_lambdas) = 1 
  
  niter = 1
  no_change = 0 
  while(niter <= max_iter & no_change < 1){
    cat("Beginning iteration", niter, "\n")
    cat("Computing abundance p-values \n")
    unique_partitions = unique(partition_df$partition) %>% as.numeric() %>% sort()
    
    isomiR_seqs = partition_df$uniqueSequence[partition_df$center == 0]
    raw_p = vector(length = length(isomiR_seqs))
    names(raw_p) = isomiR_seqs
    
    for(j in unique_partitions){
      nj = filter(partition_df, partition == j & center == 1) %>% select(., count) %>% unlist() %>% unname()
      cat("Partition", j, "center sequence has observed read count", nj, "\n")
      partition_isomiRs = filter(partition_df, partition == j & center == 0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
      cat("Partition", j, "contains", length(partition_isomiRs), "non center sequences. \n")
      for(i in partition_isomiRs){
        #cat("seq:", i, "\n")
        ni = filter(partition_df, uniqueSequence == i) %>% select(., count) %>% unlist()
        #cat("lambda:",lambda, "\n")
        lambda = master_lambdas[[j]][[i]]
        num = ppois(ni-1, nj*lambda, lower.tail=FALSE)
        denom = 1-dpois(0, nj*lambda)
        #cat("Raw p-value:", num/denom, "\n")
        raw_p[i] = num/denom
      }
    }
    
    cat("Adjusting for multiple testing and getting hypothesis testing results \n")
    
    if(adjust_method == "BH"){
      cat("Adjusting for multiple testing using Benjamini-Hochberg procedure\n")
      adjusted_p = p.adjust(raw_p, method="BH")
      results = rep(0, length(adjusted_p))
      results[adjusted_p < omega_A] = 1
      results_df = cbind(raw_p=raw_p, adj_p=adjusted_p, sig = results) %>% data.frame()
    } else if(adjust_method == "Bonferroni"){
      cat("Adjusting for multiple testing using Bonferroni\n")
      results = rep(0, length(raw_p))
      results[raw_p < omega_A/length(results)] = 1
      results_df = cbind(adj_p=raw_p, sig = results) %>% data.frame()
    }
    
    #now we need to use the results to potentially ID a new center sequence for each existing partition
    
    #first - create copy of partition_df for updating 
    update_df = partition_df 
    
    for(j in unique_partitions){
      #for each partition, get list of partition isomiRs 
      partition_isomiRs = filter(partition_df, partition == j  & center == 0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
      #create subset of hypothesis testing results for partition isomiRs 
      results_subset_df = results_df[partition_isomiRs,]
      # filter subset for significant p-values only if we have more than 0 rows in results_subset_df
      if(nrow(results_subset_df) > 0){
        new_center_seq = filter(results_subset_df, sig == 1)
        if(nrow(new_center_seq) > 0){
          new_center_seq = filter(new_center_seq, adj_p==min(adj_p)) %>% row.names()
          new_partition = max(update_df$partition) + 1 
          cat("Creating partition", new_partition, "from partition", j, ". Checking for new center sequence to  create partition. \n")
          if(length(new_center_seq) > 1){
            cat("Multiple candidates for new center sequence. Picking one at random.\n")
            new_center_seq = sample(new_center_seq, 1)
          }
          update_df$partition[update_df$uniqueSequence == new_center_seq] = new_partition
          update_df$center[update_df$uniqueSequence == new_center_seq] = 1
        }
      }
    }
    
    #after IDing new center sequences / creating new partitions, we need to update our list of unique partitions and center sequences
    center_seqs = update_df$uniqueSequence[update_df$center == 1]
    
    #check to see if we've created any new partitions: 
    if(setequal(unique_partitions, unique(update_df$partition))==TRUE){
      cat("No new partitions created.\n")
      cat("Ending at iteration", niter, "\n")
      return(list(partition_df=partition_df, initial_center_seq=initial_center_seq, alignments=alignments, transitions=transitions, results_subset_df=results_subset_df, raw_p=raw_p, results_df=results_df, master_lambdas=master_lambdas, niter=niter, no_change=no_change, update_df=update_df))
    }
    
    unique_partitions = unique(update_df$partition) %>% as.numeric() %>% sort()
    
    #get list of sequences with significant p-values so we can consider those sequences for new group membership
    
    significant_seqs = row.names(results_df[results_df$sig==1,])
    significant_seqs = significant_seqs[!(significant_seqs %in% center_seqs)]
    
    for(j in unique_partitions){
      if(!(j %in% names(master_lambdas))){
        center_seq = filter(update_df, partition == j & center == 1) %>% select(., uniqueSequence) %>% unlist() %>% unname()
        new_alignments = lapply(significant_seqs, Biostrings::pairwiseAlignment, pattern = center_seq)
        names(new_alignments) = significant_seqs
        new_transitions = lapply(new_alignments, get_transitions)
        new_lambdas = lapply(new_transitions, compute_lambda, transition_probs = transition_probability_matrix)
        master_lambdas[[length(master_lambdas)+1]] = new_lambdas
      }
    }
    names(master_lambdas) = unique_partitions
    
    #now we need to determine which partition significant sequences are joining 
    
    for(i in significant_seqs){
      l = lapply(master_lambdas, function(x) return(x[[i]])) %>% unlist()
      new_p = names(l[which(l==max(l))]) %>% as.numeric()
      if(length(new_p) > 1){
        new_p = sample(new_p, 1)
      }
      update_df$partition[update_df$uniqueSequence == i] = new_p
    }
    
    cat("Finishing iteration", niter, "\n")
    niter = niter+1
    partition_df=update_df
  }
  #return everything we can rn because we're in the development stages (ugh lol)
  return(list(partition_df=partition_df, initial_center_seq=initial_center_seq, alignments=alignments, 
              transitions=transitions, results_subset_df=results_subset_df, raw_p=raw_p, results_df=results_df, 
              master_lambdas=master_lambdas, niter=niter, no_change=no_change, update_df=update_df, new_center_seq=new_center_seq))
}
