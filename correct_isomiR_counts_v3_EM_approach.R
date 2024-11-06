

correct_isomiR_counts_v3 = function(user_miRNA, rowdata, countdata, max_iterations, OMEGA=0.05, sample_idx=2){
  
  cat("Initializing partition_df object\n")
  partition_df = cbind(rowdata, count = countdata[,sample_idx]) %>% data.frame() %>% filter(., miRNA_name == user_miRNA)
  partition_df$partition = rep(1, nrow(partition_df))
  initial_center_seq = filter(partition_df, count == max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()
  partition_df$center = rep(0, nrow(partition_df))
  partition_df$center[partition_df$uniqueSequence == initial_center_seq] = 1
  
  cat("Calculating initial transition probabilities matrix using initial partitioning\n")
  transition_counts = initialize_transition_counts_matrix(5)
  isomiR_seqs = filter(partition_df, center == 0 & count != 0 & partition == 1) %>% select(., uniqueSequence) %>% unlist()
  
  cat("Getting initial alignments\n")
  initial_alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = initial_center_seq)
  names(initial_alignments) = isomiR_seqs
  cat("Getting transitions from those alignments\n")
  initial_transitions = lapply(initial_alignments, get_transitions)
  
  for(t in initial_transitions){
    alignment_length = ncol(t)
    for(i in 1:alignment_length){
      trns = get_transition_at_idx(i, t)
      x = trns$x
      y = trns$y
      if(i == 1 & x != "-" & y != "-"){
        transition_counts["-", "-"] + transition_counts["-", "-"] + 1 
      } else if(i == alignment_length & x != "-" & y != "-"){
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
  
  iter = 0 
  inner_loop_no_change = 0
  
  #get initial alignments 
  alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = initial_center_seq)
  names(alignments) = isomiR_seqs
  #get initial transitions
  transitions = lapply(alignments, get_transitions)
  #get initial lambdas 
  lambdas = lapply(transitions, compute_lambda, transition_probs=transition_probs)
  #use those initial lambdas we get by aligning all isomiR sequences to initialize master_lambdas list 
  master_lambdas = list(lambdas)
  names(master_lambdas) = initial_center_seq
  
  while(inner_loop_no_change == 0 & iter < max_iterations){
    #get unique partitions
    unique_partitions = unique(partition_df$partition)
    #calculate raw p values 
    raw_p_values = vector()
    for(j in unique_partitions){
      n_j = filter(partition_df, partition == j & center == 1) %>% select(., count) %>% unlist() %>% unname()
      cat("In partition", j, "observed read count of center sequence is", n_j, "\n") 
      partition_isomiRs = filter(partition_df, partition == j & center == 0 & count != 0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
      for(i in partition_isomiRs){
        lambda = master_lambdas[[j]][[i]]
        n_i = filter(partition_df, uniqueSequence == i) %>% select(., count) %>% unlist() %>% unname()
        num = ppois(n_i, n_j*lambda, lower.tail=FALSE)
        denom = 1-dpois(0, n_j*lambda)
        p = c(num/denom)
        names(p) = i
        raw_p_values = c(raw_p_values, p)
      }
    }
    #get hypothesis testing results 
    results_df = data.frame(raw_p_values)
    results_df$adjusted_p_values = p.adjust(raw_p_values, method = "BH")
    results_df$significant = rep(0, nrow(results_df))
    results_df$significant[results_df$adjusted_p_values < OMEGA_A] = 1 
    
    #get significant sequences 
    significant_isomiRs = row.names(results_df[results_df$significant==1,])
    
    #for each partition, determine if we need to create a new partition:
    
    #first create update_df, a copy of the current partition_df 
    update_df = partition_df 
    
    cat("Determining which partitions will spawn new partitions \n")
    for(j in unique_partitions){
      #get isomiRs in partition that have significant p-value
      significant_partition_isomiRs = filter(partition_df, partition == j & center == 0 & count != 0) %>% 
        filter(., uniqueSequence %in% significant_isomiRs) %>% select(., uniqueSequence) %>% unlist() %>% unname()
      #if there are any, then create new partition:
      if(length(significant_partition_isomiRs) > 0){
        #create new grp 
        new_partition = max(update_df$partition) + 1
        #ID new center sequence
        new_center = filter(partition_df, partition == j & uniqueSequence %in% significant_isomiRs) %>% filter(., count == max(count)) %>% 
          select(., uniqueSequence) %>% unlist() %>% unname()
        new_center = new_center[1]
        
        update_df$center[update_df$uniqueSequence == new_center] = 1
        update_df$partition[update_df$uniqueSequence == new_center] = new_partition
        
      }
      
    }
    
    #now update list of center sequences:
    cat("Updating center sequence vector \n")
    center_seqs = filter(update_df, center == 1) %>% select(., uniqueSequence) %>% unlist() %>% unname()
    cat("Updating isomiR sequence vector \n")
    isomiR_seqs = filter(update_df, center == 0 & count!= 0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
    
    cat("Aligning isomiR sequences to newly identifed center sequences and calculating lambdas \n")
    for(seq in center_seqs){
      if(!(seq %in% names(master_lambdas))){
        new_alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = seq)
        names(new_alignments) = isomiR_seqs
        new_transitions = lapply(new_alignments, get_transitions)
        new_lambdas = lapply(new_transitions, compute_lambda, transition_probs=transition_probs)
        master_lambdas[[seq]] = new_lambdas
      }
    }

    cat("Updating group assignment, assigning isomiRs to newly created partitions \n")
    for(seq in isomiR_seqs){
      l = lapply(master_lambdas, function(x) return(x[[seq]])) %>% unlist()
      l = which(l == max(l))
      update_df$partition[update_df$uniqueSequence == seq] = l
    }

    if(all.equal(update_df$partition, partition_df$partition) == TRUE){
      inner_loop_no_change = 1 
      cat("All technical isomiRs removed \n")
    } else{
      partition_df = update_df 
      cat("Finishing iteration", iter, "\n")
      iter = iter+1
      cat("Beginning iteration", iter, "\n")
    }
  }
    
    
  
  return(list(partition_df=partition_df, transition_probs=transition_probs, results_df=results_df, center_seqs=center_seqs, isomiR_seqs=isomiR_seqs,
              master_lambdas=master_lambdas, update_df=update_df))
   
}

  