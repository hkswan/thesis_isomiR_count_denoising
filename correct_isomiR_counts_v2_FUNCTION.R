#function to correct and remove all technical isomiRs from sequencing data 
#assumes read counts of technical isomiRs are generated from a true sequence by process following a Poisson distribution 
#takes following arguments: test_miRNA, a string giving the name of the miRNA whose technical isomiRs we'd like
#to remove from the data 
#rowdata, a dataframe containing mapping information from running FASTQ files through miRge/sRNAbench
#must have column titled miRNA_name indicating name of miRNA sequence was mapped to and uniqueSequence containing isomiR sequence
#countdata is a dataframe of read counts with a row for each unique sequence and a column for each sample
#sample_idx is an integer indicating which columns counts' we want to denoise and use in error modeling
#max_iterations is an integer, maximum number of iterations algorithm will run for without removing all technical isomiRs from count data
#OMEGA_A is a real number between 0 and 1, acts like alpha in traditional hypothesis testing framework
#maximum probability of an isomiR sequence read count as extreme or more extreme than the true sequence read count and still 
#consideer isomiR sequence its own true sequence 
#returns a dataframe called partition_df, 1 row for each unique sequence, partition column indicates which group
#each sequence belongs to, used to collapse counts 
#function uses some helper fuunctions defined in file in following path
#/scratch/hswan/thesis_lengthvariant_ID/correct_technical_length_variant_functions.R

correct_isomiR_counts_v2 = function(test_mirna, rowdata, countdata, sample_idx, max_iterations, OMEGA_A){
  #INITIALIZE PARTITION_DF W/FOLLOWING STEPS
  #get sequence, miRNA, and countdata from a single sample together in a dataframe  - currently we have multiple miRNAs in partition_df
  partition_df = cbind(rowdata, count = countdata[,sample_idx]) %>% data.frame()
  #filter partition_df so we have information for just the user-specified miRNA 
  partition_df = filter(partition_df, miRNA_name == test_mirna)
  #create the initial partition, group 1.  all sequences mapping to the same miRNA start in this initial partition 
  partition_df$partition = rep(1, nrow(partition_df))
  #create indicator variable called center, if 1 that indicates that sequence is center of corresponding partition 
  partition_df$center = rep(0, nrow(partition_df))
  #ID initial center sequence - most abundant sequence mapping to miRNA. is often but not always canonical sequence 
  initial_center_seq = filter(partition_df, count == max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()
  #update center indicator variable for our initial center seq to be 1 
  partition_df$center[partition_df$uniqueSequence == initial_center_seq] = 1 
  #adding a column for sequence lengths for now because i want to see the amount of degradation we have 
  partition_df = mutate(partition_df, seqlen = nchar(uniqueSequence))
  
  unique_isomirs = unique(partition_df$uniqueSequence)
  
  for(seq in unique_isomirs){
    if(filter(partition_df, uniqueSequence == seq) %>% nrow() > 1){
      row_idx_to_rmv = which(partition_df$uniqueSequence == seq & partition_df$count == 0)
      if(length(row_idx_to_rmv) == filter(partition_df, uniqueSequence == seq) %>% nrow()){
        row_idx_to_rmv = row_idx_to_rmv[-1]
      }
      partition_df = partition_df[-row_idx_to_rmv,]
    }
  }
  
  cat("Initial partition created. All isomiR sequences mapping to same miRNA belong to same initial partition \n")
  cat("Center sequence of initial partition is", initial_center_seq, "\n")
  
  no_change = 0 
  iter = 1 
  
  while(iter < max_iterations & no_change == 0){
    #STEP 1 - identify unique partitions
    unique_partitions = unique(partition_df$partition)
    cat("Unique partitions of isomiR sequences are", unique_partitions, "\n")
    #for each partition perform pairwise alignments between center sequence and all partition isomiRs 
    #if center sequence is only element of partition then can skip performing alignments
    cat("Performing alignments\n")
    alignments = list()
    for(j in unique_partitions){
      center_seq = filter(partition_df, partition == j & center == 1) %>% select(., uniqueSequence) %>% unlist() %>% unname()
      partition_isomiRs = filter(partition_df, partition == j & center == 0 & count != 0) %>% select(., uniqueSequence) %>% unlist() %>% unname()
      if(length(partition_isomiRs) > 0){
        a = lapply(partition_isomiRs, Biostrings::pairwiseAlignment, pattern = center_seq)
        names(a) = partition_isomiRs
        alignments = c(alignments, a)
      }
    }
    cat("Getting transitions from pairwise alignments\n")
    transitions = lapply(alignments, get_transitions)
    
    cat("Estimating transition probabilities from the data\n")
    transition_counts = initialize_transition_counts_matrix(5)
    cat("Initialized transition_counts matrix:\n")
    print(transition_counts)
    for(t in transitions){
      alignment_length = ncol(t)
      for(i in 1:alignment_length){
        trs = get_transition_at_idx(i, t)
        x=trs$x
        y=trs$y
        if(i == 1 & x != '-' & y != '-'){
          transition_counts['-','-'] = transition_counts['-','-']+1
        } else if(i == alignment_length & x != '-' & y != '-'){
          transition_counts['-','-'] = transition_counts['-','-']+1
        } else{
          transition_counts[x,y] = transition_counts[x,y]+1
        }
      }
    }
    cat("transition_counts matrix:\n")
    print(transition_counts)
    
    transition_probs = transition_counts
    for(i in 1:nrow(transition_probs)){
      transition_probs[i,] = transition_probs[i,]/sum(transition_probs[i,])
    }
    cat("transition_probs:\n")
    print(transition_probs)
    
    cat("Calculating lambdas\n")
    lambdas = lapply(transitions, compute_lambda, transition_probs=transition_probs)
    
    cat("Calculating abundance p-values\n")
    raw_p_values = vector()
    for(j in unique_partitions){
      n_j = filter(partition_df, center == 1 & partition == j) %>% select(., count) %>% unlist()
      cat("Calculating abundance p-values in partition", j, "\n")
      cat("n_j is", n_j, "\n")
      partition_isomiRs = filter(partition_df, center == 0 & count != 0 & partition == j) %>% select(., uniqueSequence) %>% unlist()
      for(i in partition_isomiRs){
        n_i = partition_df$count[partition_df$uniqueSequence == i]
        lambda = lambdas[[i]]
        num = ppois(n_i, n_j*lambda, lower.tail=FALSE)
        denom = 1-dpois(0, n_j*lambda)
        p = c(num/denom)
        names(p) = i
        raw_p_values = c(raw_p_values, p)
      }
    }
    
    cat("Getting hypothesis testing results\n")
    results_df = data.frame(raw_p_values)
    results_df$adjusted_p_values = p.adjust(raw_p_values, method='BH')
    results_df$seq = names(raw_p_values)
    row.names(results_df)=NULL
    results_df$significant = rep(0, nrow(results_df))
    results_df$significant[results_df$adjusted_p_values < OMEGA_A] = 1
    
    cat("Determining if new partitions need to be added\n")
    update_df=partition_df
    for(j in unique_partitions){
      partition_isomiRs = filter(partition_df, count !=0 & center == 0 & partition == j) %>% select(., uniqueSequence) %>% unlist()
      significant_isomiRs = filter(results_df, seq %in% partition_isomiRs) %>% filter(., significant == 1) %>% select(., seq) %>% unlist()
      if(length(significant_isomiRs) > 0){
        #ID new center sequence 
        new_center = filter(partition_df, uniqueSequence %in% significant_isomiRs) %>% filter(., count == max(count)) %>% select(., uniqueSequence) %>% unlist()
        new_center = new_center[1]
        #create new partition
        new_partition = max(update_df$partition) + 1
        update_df$partition[update_df$uniqueSequence == new_center] = new_partition
        update_df$center[update_df$uniqueSequence == new_center] = 1 
        
        #align isomiR sequences with new center
        new_alignments = lapply(partition_isomiRs, Biostrings::pairwiseAlignment, pattern=new_center)
        names(new_alignments) = partition_isomiRs
        new_transitions = lapply(new_alignments, get_transitions)
        new_lambdas = lapply(new_transitions, compute_lambda, transition_probs=transition_probs)
        for(i in partition_isomiRs){
          old_lambda = lambdas[[i]]
          new_lambda = new_lambdas[[i]]
          if(new_lambda > old_lambda){
            update_df$partition[update_df$uniqueSequence == i] = new_partition
          }
        }
      }
    }
    
    if(all.equal(partition_df$partition, update_df$partition)==TRUE){
      no_change = 1 
      cat("All technical isomiRs have been identified\n")
      partition_df = update_df
      return(partition_df)
    }
    
    cat("Iteration", iter, "completed\n")
    partition_df=update_df
    iter = iter + 1
  }
  return(list(partition_df=partition_df, iter=iter))}