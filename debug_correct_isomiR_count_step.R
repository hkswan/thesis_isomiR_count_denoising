
#clear environment
rm(list=ls())

#load necessary libraries
library(Biostrings)
library(tidyverse)

#load benchmark dataset function
source("/scratch/hswan/thesis_isomiR_count_denoising/load_mouse_miRNA_data_function.R")
#load other necessary functions
source("/scratch/hswan/thesis_isomiR_count_denoising/correct_technical_length_variant_functions.R")

# #load data
data = load_mouse_miRNA_data()

isomiR_se_object=data$isomiR_se_object
countdata=data$countdata
rowdata=data$rowdata
print(head(rowdata))

max_iter = 2
OMEGA_A = 0.05
test_miRNA = unique(rowdata$miRNA_name)[2]
cat("miRNA for testing is", test_miRNA, "\n")

initial_partition_df = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_partition_df.Rds")
print(head(initial_partition_df))

transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/transition_probs.Rds")
print(transition_probs)

correct_isomiR_counts_step = function(test_miRNA, initial_partition_df, transition_probs, max_iter, OMEGA_A){
  iter =  0
  no_change = 0
  partition_df = initial_partition_df
  while(iter < max_iter & no_change == 0){
    unique_partitions = unique(partition_df$partition)
    print(length(unique_partitions))
    if(length(unique_partitions) == 1){
      center_seq = filter(partition_df, center == 1) %>% select(., uniqueSequence)  %>% unlist() %>% unname()
      cat("Center seq is", center_seq, "\n")
      isomiR_seqs = filter(partition_df, center == 0 & count != 0) %>% select(., uniqueSequence) %>% unlist()  %>% unname()
      #alignments
      cat("Getting alignments \n")
      alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = center_seq)
      names(alignments)  = isomiR_seqs
      #transitions
      cat("Getting transitions\n")
      transitions = lapply(alignments, get_transitions)
      #lambdas
      cat("Calculating lambdas\n")
      lambdas = lapply(transitions, compute_lambda, transition_probs=transition_probs)
      master_lambdas = list()
      master_lambdas[[center_seq]] = lambdas
      cat("Calculating abundance p_values\n")
      #abundance p values:
      n_j = filter(partition_df, uniqueSequence  == center_seq) %>% select(., count) %>% unlist() %>% unname()
      raw_p_values = vector()
      for(i in isomiR_seqs){
        lambda = master_lambdas[[center_seq]][[i]]
        n_i = filter(partition_df, uniqueSequence  == i) %>% select(., count) %>% unlist() %>% unname()
        num = ppois(n_i, n_j*lambda, lower.tail=FALSE)
        denom  = 1-dpois(0, n_j*lambda)
        p = c(num/denom)
        names(p) = i
        raw_p_values = c(raw_p_values, p)
      }
      cat("Getting hypothesis testing results \n")
      results_df = data.frame(raw_p_values)
      results_df$adjusted_p_values = p.adjust(raw_p_values, method="BH")
      results_df$significant = rep(0, nrow(results_df))
      results_df$significant[results_df$adjusted_p_values < OMEGA_A] = 1

      cat("Determining if we need to create a new partition\n")
      significant_seqs = row.names(results_df[results_df$significant==1,])
      update_df = partition_df
      #if at least one seq in initial partition has significant p-value, start new partition
      if(length(significant_seqs) != 0){
        new_center_seq = filter(partition_df, uniqueSequence %in% significant_seqs) %>% filter(., count == max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()
        new_partition = max(update_df$partition) + 1
        update_df$center[update_df$uniqueSequence == new_center_seq] = 1
        update_df$partition[update_df$uniqueSequence == new_center_seq] = new_partition
        #update center sequences list
        center_seqs = filter(update_df, center == 1) %>% select(., uniqueSequence) %>% unlist() %>% unname()
        isomiR_seqs = filter(update_df, center == 0 & count !=0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
        for(seq in center_seqs){
          if(!(seq %in% names(master_lambdas))){
            new_alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = seq)
            names(new_alignments) = isomiR_seqs
            new_transitions = lapply(new_alignments, get_transitions)
            new_lambdas = lapply(new_transitions, compute_lambda, transition_probs=transition_probs)
            master_lambdas[[seq]] = new_lambdas
          }
        }
        cat("Updating partition membership, letting some isomiR sequences move to newly created partition\n")
        for(seq in isomiR_seqs){
          l = lapply(master_lambdas, function(x) return(x[[seq]])) %>% unlist()
          l = which(l == max(l))
          update_df$partition[update_df$uniqueSequence == seq] = l

        }
        partition_df = update_df
        cat("Iter", iter, "done\n")
        iter = iter + 1
        cat("Beginning iter", iter,  "\n")
      } else{
        no_change = 1
      }

     
    } else{
      raw_p_values  = vector()
      for(j in unique_partitions){
        center_seq = filter(partition_df, partition ==  j & center == 1) %>% select(., uniqueSequence) %>% unlist() %>% unname()
        cat("Center sequence of partition", j, "is", center_seq, "\n")
        partition_isomiRs = filter(partition_df, partition == j & center == 0 & count != 0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
        n_j = filter(partition_df, uniqueSequence == center_seq) %>% select(., count) %>% unlist() %>% unname()
        cat("Observed read count of center sequence is",  n_j, "\n")
        for(i in partition_isomiRs){
          lambda = master_lambdas[[center_seq]][[i]]
          n_i = filter(partition_df, uniqueSequence ==i) %>% select(., count) %>% unlist() %>% unname()
          num = ppois(n_i, n_j*lambda, lower.tail=FALSE)
          denom = 1-dpois(0, n_j*lambda)
          p=c(num/denom)
          names(p) = i
          raw_p_values = c(raw_p_values, p)
        }
      }
      cat("Getting results of hypothesis testing \n")
      results_df = data.frame(raw_p_values)
      results_df$adjusted_p_values = p.adjust(raw_p_values)
      results_df$significant = rep(0, nrow(results_df))
      results_df$significant[results_df$adjusted_p_values < OMEGA_A] = 1

      cat("Determining which existing partitions, if any, have sequences that should start their own partition\n")

      significant_seqs = row.names(filter(results_df, significant == 1))
      #make copy of partition_df for updating
      update_df = partition_df

      if(length(significant_seqs) > 0){
        for(j in unique_partitions){
          partition_isomiRs = filter(partition_df, partition == j & center == 0 & count !=0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
          partition_significant_seqs = significant_seqs[significant_seqs %in% partition_isomiRs]
          if(length(partition_significant_seqs) > 0){
            #ID new center seq
            new_center_seq = filter(partition_df, uniqueSequence %in% partition_significant_seqs) %>% filter(., count == max(count)) %>% select(., uniqueSequence)  %>% unlist() %>% unname()
            new_center_seq = new_center_seq[1]
            #start new grp
            new_partition = max(update_df$partition) + 1
            update_df$partition[update_df$uniqueSequence == new_center_seq] = new_partition
            update_df$center[update_df$uniqueSequence == new_center_seq] = 1

          }
        }
        #update list of center sequences
        center_seqs = filter(update_df, center == 1) %>% select(., uniqueSequence) %>% unlist() %>% unname()
        #update list of isomiR sequences
        isomiR_seqs = filter(update_df, center == 0 & count != 0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
        #align isomiR sequences to newly identified center sequences, then calculate lambdas so we can update grp membership
        for(seq in center_seqs){
          if(!(seq %in% names(master_lambdas))){
            new_alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = seq)
            names(new_alignments) = isomiR_seqs
            new_transitions = lapply(new_alignments, get_transitions)
            new_lambdas = lapply(new_transitions, compute_lambda, transition_probs=transition_probs)
            master_lambdas[[seq]] = new_lambdas
          }

        }
        #update partition membership
        for(i in isomiR_seqs){
          l = lapply(master_lambdas, function(x) return(x[[i]])) %>% unlist()
          l = which(l==max(l))
          update_df$partition[update_df$uniqueSequence  == i] = l
        }
        partition_df=update_df
        cat("Finishing iteration", iter, "\n")
        iter = iter + 1
        cat("Beginning iteration", iter, "\n")
      } else{
        cat("No new partitions being formed. No change from last iteration \n")
        no_change = 1
      }
    }

  }
  return(list(partition_df=partition_df, master_lambdas=master_lambdas))
}

tst = correct_isomiR_counts_step(test_miRNA = test_miRNA, initial_partition_df = initial_partition_df, max_iter = 10,
                           OMEGA_A=0.05, transition_probs = transition_probs)



cat("Successfully loaded function in environment \n")