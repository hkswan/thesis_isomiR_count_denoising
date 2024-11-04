
#takes output of pairwise alignment of partition center sequence and isomiR sequence as input
#iterates over each position of the alignment of the two sequences and returns the nucleotide at that position 
#in pattern sequence and subject sequence 
#output is a 2 x L dimensional matrix where L is the length of the alignment 
#output intended to be used as input in computing lamba for a given center sequence / isomiR sequence pair
get_transitions = function(alignment){
  pattern = Biostrings::alignedPattern(alignment)
  subject = Biostrings::alignedSubject(alignment)
  alignment_length = nchar(pattern)
  transitions = sapply(1:alignment_length, function(x) return(c(Biostrings::substr(pattern, x, x), Biostrings::substr(subject, x,x))))
  return(transitions)
}


#takes output of get_transitions function as input, along with matrix of transition probabilities
#iterate through the columns, element in first row gives row name of transition probability value we need
#element in 2nd row gives column name of transition probability value we need 
#then multiply all transition probabilities together to get lambda and return 

get_individual_transition_prob = function(idx, transition, transition_probs){
  x = transition[1,idx]
  y = transition[2,idx]
  pr = transition_probs[x,y]
  return(pr)
}


compute_lambda = function(transition, transition_probs){
  alignment_length = ncol(transition)
  pr = sapply(1:alignment_length, get_individual_transition_prob, transition=transition, transition_probs=transition_probs)
  lambda = prod(pr)
  return(lambda)
}

get_obs_read_count = function(seq, countdata, sample_idx, rowdata, miRNA){
  df = cbind(rowdata, countdata[,sample_idx])
  colnames(df) = c(colnames(rowdata), 'count')
  n = filter(df, uniqueSequence == seq & exact.miRNA == miRNA) %>% select(., count) %>% unlist()
  return(n)
}

compute_abundance_p_value = function(lambdas, partition_df, seq){
  #browser()
  lambda = lambdas[seq] %>% unlist()
  n_i = partition_df$count[partition_df$seq == seq]
  partition = partition_df$partition[partition_df$seq ==  seq]
  n_j = partition_df$count[partition_df$partition == partition & partition_df$center == 1]
  denom = 1-dpois(0, n_j*lambda)
  num = ppois(n_i, n_j*lambda, lower.tail=FALSE)
  if(n_i == 0){
    p  = 1 
  }
  # if(denom == 0){
  #   p = 1
  # }
  else{
    p = num / denom
  }
  return(p)
}

get_hypothesis_test_results = function(p_values, OMEGA_A, isomiR_seqs){
  corrected_p_values  = p.adjust(p_values)
  results_matrix = cbind(p_values, corrected_p_values)
  results_matrix = cbind(results_matrix, rep(0, nrow(results_matrix)))
  results_matrix = data.frame(results_matrix)
  results_matrix[results_matrix[,2] < OMEGA_A,3] = 1
  colnames(results_matrix) = c("raw_p_values", "adjusted_p_values", "significance_indicator")
  results_matrix$seq = names(p_values)
  #row.names(results_matrix) = isomiR_seqs
  return(results_matrix)
}

update_partition = function(results, partition_df){
  unique_partitions = unique(partition_df$partition)
  update_partition_df = partition_df
  significant_seqs = filter(results, significance_indicator == 1) %>% row.names()
  for(i in 1:length(unique_partitions)){
    seqs_to_move = filter(update_partition_df, seq %in% significant_seqs & partition == i) %>% select(., seq) %>% unlist() %>% unname()
    new_partition = max(unique(update_partition_df$partition)) + 1
    update_partition_df$partition[update_partition_df$seq %in% seqs_to_move] = new_partition
  }
  return(update_partition_df)
}

#draws an epsilon from a beta(alpha, beta) distribution and uses that to create a 5 x 5 matrix of nucleotide transition probabilities, including probability
#of a nucleotide to gap transition 
update_transition_probabilities_matrix = function(seed, alpha, beta){
  set.seed(seed)
  epsilon = rbeta(1, alpha, beta)
  transition_probs = matrix(epsilon/4, 5, 5)
  row.names(transition_probs) = c("A", "C", "G", "T", "-")
  colnames(transition_probs) = c("A", "C", "G", "T", "-")
  diag(transition_probs) = 1-epsilon
  return(transition_probs)
}

get_transition_at_idx = function(idx, transition){
  return(list(x=transition[1,idx], y=transition[2,idx]))
}

initialize_transition_counts_matrix = function(shape){
  mat = matrix(0, shape, shape)
  row.names(mat) = c("A", "C", "G", "T", "-")
  colnames(mat) = row.names(mat)
  return(mat)
}