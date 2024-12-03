
library(Biostrings)
library(tidyverse)
library(DESeq2)

args = commandArgs(trailingOnly = TRUE)
miRNA_idx = as.numeric(args[1])

source("/scratch/hswan/thesis_isomiR_count_denoising/correct_technical_length_variant_functions.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/load_mouse_miRNA_data_function.R")

data = load_mouse_miRNA_data()

rowdata = data$rowdata
countdata = data$countdata

uniq_miRNAs = unique(rowdata$miRNA_name)
miRNA = uniq_miRNAs[miRNA_idx]

count_df = cbind(rowdata, count = rowSums(countdata))

isomiR_subset = filter(count_df, miRNA_name == miRNA)

center_seq = filter(isomiR_subset, count == max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()
isomiR_seqs = filter(isomiR_subset, uniqueSequence != center_seq) %>% select(., uniqueSequence) %>% unlist() %>% unname()

alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = center_seq)

transitions = lapply(alignments, get_transitions)

transition_counts = initialize_transition_counts_matrix(5)

if(length(transitions) != 0){
  for(t in transitions){
    alignment_length = ncol(t)
    for(i in 1:alignment_length){
      trns = get_transition_at_idx(i, t)
      x=trns$x
      y=trns$y
    if(i == 1 & x != "-" & y != "-"){
      transition_counts[x,y]=transition_counts[x,y]+1
      transition_counts['-', '-'] = transition_counts['-', '-']+1
    }
    if(i == alignment_length & x != "-" & y !=  "-"){
      transition_counts['-', '-'] = transition_counts['-', '-'] + 1 
      transition_counts[x,y] = transition_counts[x,y]+1
    } else{
      transition_counts[x,y]=transition_counts[x,y]+1
    }
    
    }
  }
}


filename = paste0(miRNA, "_transition_counts.RDS")
saveRDS(transition_counts, paste0('/scratch/hswan/thesis_isomiR_count_denoising/transition_count_matrices/', filename))

