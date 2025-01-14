library(tidyverse)
library(Biostrings)


source("/scratch/hswan/thesis_isomiR_count_denoising/denoise_isomiR_counts_WORKING_FUNCTION.R")

simulate_length_variant_isomiRs_extened = function(rowdata, countdata, miRNA, seed, max_iter, N, transition_probs = NULL, verbose = FALSE){
  cat("Identifying center sequence for miRNA", miRNA, "\n")
  center_seq = filter(count_df, miRNA_name == miRNA) %>% filter(., count == max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()
  isomiR_seqs = vector()
  set.seed(seed)
  seq_length_diffs = -7:7
  ends = c("3p", "5p")
  nucleotides = c("A", "C", "G", "T")
  differences = vector()
  iter_counter = 0 
  max_iter = max_iter
  N = N
  
  cat("Generating length variant isomiR sequences \n")
  while(length(isomiR_seqs) < N & iter_counter < max_iter){
    diff = sample(seq_length_diffs, 1)
    
    seq = center_seq
    #if difference in length is negative that means there are additional nt in the isomiR sequence when compared to the center sequence 
    if(diff < 0){
      for(i in 1:-diff){
        #first, sample a nucleotide at random:
        nt = sample(nucleotides, 1)
        end = sample(ends, 1)
        if(verbose){
          cat("nt:", nt, "\n")
          cat("end:", end, "\n")
        }
        if(end == "3p"){
          seq = paste0(seq, nt, collapse="")
        } else if(end == "5p"){
          seq = paste0(nt, seq, collapse="")
        }
      }
      if(!(seq %in% isomiR_seqs)){
        #cat("Adding seq to isomiR_seqs")
        isomiR_seqs = c(isomiR_seqs, seq)
      }
      iter_counter = iter_counter + 1
    } else if(diff > 0){
      for(i in 1:diff){
        end = sample(ends, 1)
        if(end == "3p"){
          seq = substr(seq, start=1, stop = nchar(seq)-1)
        } else if(end == "5p"){
          seq = substr(seq, start = 2, stop = nchar(seq))
        }
      }
      if(!(seq %in% isomiR_seqs)){
        isomiR_seqs = c(isomiR_seqs, seq)
      }
      iter_counter = iter_counter + 1 
    }
  }
  
  cat("Getting alignments\n")
  alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = center_seq)
  names(alignments) = isomiR_seqs
  cat("Getting transitions from alignments \n")
  transitions = lapply(alignments, get_transitions)
  cat("Calculating lambdas \n")
  lambdas = lapply(transitions, compute_lambda, transition_probs=transition_probs)
  
  cat("Identifying observed read count of center sequence, n_j:\n")
  
  nj = filter(count_df, uniqueSequence == center_seq) %>% select(., count) %>% unlist() %>% unname()
  if(verbose){
    cat("nj:", nj, "\n")
  }
  
  mus = nj * unlist(lambdas)
  sim_counts = vector()
  cat("Simulating read counts from error distribution \n")
  for(i in 1:length(mus)){
    sim_counts = c(sim_counts, rpois(1, mus[i]))
  }
  
  cat("Assembling data \n")
  sim_rowdata = data.frame(c(center_seq,isomiR_seqs), rep(miRNA, length(isomiR_seqs)+1))
  colnames(sim_rowdata) = c("uniqueSequence", "miRNA_name")
  
  sim_counts = c(nj, sim_counts)
  return(list(sim_counts=sim_counts, sim_rowdata=sim_rowdata, mus=mus, seed=seed, N=N, iter=iter_counter))
}