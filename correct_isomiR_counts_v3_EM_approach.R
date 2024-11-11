

correct_isomiR_counts_v3 = function(user_miRNA, rowdata, countdata, inner_loop_max_iterations=NULL, outer_loop_max_iterations, OMEGA=0.05, sample_idx=2){
  
  cat("Initializing partition_df object\n")
  partition_df = cbind(rowdata, count = countdata[,sample_idx]) %>% data.frame() %>% filter(., miRNA_name == user_miRNA)
  partition_df$partition = rep(1, nrow(partition_df))
  initial_center_seq = filter(partition_df, count == max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()
  partition_df$center = rep(0, nrow(partition_df))
  partition_df$center[partition_df$uniqueSequence == initial_center_seq] = 1
  
  isomiR_seqs = filter(partition_df, center == 0 & count != 0 & partition == 1) %>% select(., uniqueSequence) %>% unlist()
  
  
  outer_iter = 0 
  inner_iter = 0
  outer_loop_no_change = 0 
  inner_loop_no_change = 0
  
  #initialize master_lambdas list 
  master_lambdas = list()
  master_lambdas[[initial_center_seq]] = rep(0, length(isomiR_seqs))
  
  while(outer_loop_no_change == 0 & outer_iter < outer_loop_max_iterations){
    #get unique partitions 
    unique_partitions = unique(partition_df$partition)
    center_seqs = filter(partition_df, center == 1)  %>% select(., uniqueSequence) %>% unlist() %>% unname()
    isomiR_seqs = filter(partition_df, center == 0 & count !=0) %>% select(., uniqueSequence) %>% unlist() %>%
      unname()
    #get alignments between elements of each partition and that partition's center sequence:
    alignments = list()
    for(j in unique_partitions){
      center_seq = partition_df$uniqueSequence[partition_df$center == 1 & partition_df$partition == j]
      partition_isomiRs = partition_df$uniqueSequence[partition_df$center == 0 & partition_df$count !=0 & partition_df$partition == j]
      a = lapply(partition_isomiRs, Biostrings::pairwiseAlignment, pattern = center_seq)
      names(a) = partition_isomiRs
      alignments=c(alignments,a)
    }
    #get transitions 
    transitions = lapply(alignments, get_transitions)
    #initialize transition counts matrix 
    transition_counts = initialize_transition_counts_matrix(5)
    for(t in transitions){
      alignment_length = ncol(t)
      for(i in 1:alignment_length){
        trns = get_transition_at_idx(i,t)
        x = trns$x
        y = trns$y
        if(i == 1 & x!= "-" & y!= "-"){
          transition_counts["-", "-"] = transition_counts["-", "-"] + 1
        } else if(i == alignment_length & x!= "-" & y != "-"){
          transition_counts["-", "-"] = transition_counts["-", "-"] + 1
        } else{
          transition_counts[x,y] = transition_counts[x,y] + 1 
        }
      }
    }
    transition_probs = transition_counts
    for(i in 1:nrow(transition_probs)){
      transition_probs[i,] = transition_probs[i,]/sum(transition_probs[i,])
    }
    
    #for existing partitions, calculate lambdas:
    for(seq in center_seqs){
      lam = master_lambdas[[seq]]
      
    }
    outer_iter = outer_iter + 1 
  }
  
    
    
  
  return(list(partition_df=partition_df, transition_probs=transition_probs))
   
}

  