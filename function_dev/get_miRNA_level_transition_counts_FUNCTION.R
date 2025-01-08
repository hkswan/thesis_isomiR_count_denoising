
get_miRNA_level_transition_counts=function(partition_df){
  miRNA = partition_df$miRNA_name[1]
  cat("miRNA:", miRNA, "\n")
  
  transition_counts = initialize_transition_counts_matrix(5)
  
  unique_partitions = unique(partition_df$partition) %>% sort()
  
  for(j in unique_partitions){
    center_seq = filter(partition_df, center == 1  & partition == j) %>% select(., uniqueSequence) %>% unlist() %>% unname()
    partition_isomiRs = filter(partition_df, center == 0 & partition == j) %>% select(., uniqueSequence) %>% unlist()
    a = lapply(partition_isomiRs, Biostrings::pairwiseAlignment, pattern = center_seq)
    t = lapply(a, get_transitions)
    for(i in t){
      alignment_length = ncol(i)
      for(l in 1:alignment_length){
        trns = get_transition_at_idx(l, i)
        x=trns$x
        y=trns$y
        if(x != "-" & y !="-" & (l == 1 | l == alignment_length)){
          transition_counts["-","-"] = transition_counts["-", "-"]+1
        } else{
          transition_counts[x,y]=transition_counts[x,y]+1
        }
      }
    }
  }
  return(transition_counts)
}
