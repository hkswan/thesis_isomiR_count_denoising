simulate_length_variant_counts = function(miRNA, count_df, seed, transition_probs = NULL, add_noise=FALSE, max_k =  1000){
  #first identify center sequence
  center_seq = filter(count_df, miRNA_name == miRNA) %>% filter(., count==max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()
  #initialize isomiR_seqs vector:
  isomiR_seqs = vector()
  set.seed(seed)
  #create vector of differences in length between center sequence and isomiR sequences to sample from:
  #(this is currently set aribitrarily but we can change this)
  seq_length_diffs = 1:7
  ends = c("3p", "5p")
  differences = vector()
  k=0
  cat("Generating length variant sequences \n")
  while(length(isomiR_seqs) < 50 & k <= max_k){
    diff = sample(seq_length_diffs, 1)
    #cat("Difference:", diff, "\n")
    
    #for each difference:
    seq = center_seq
    for(i in 1:diff){
      #cat("i:", i, "\n")
      #sample end of miRNA sequence:
      end = sample(ends, 1)
      #cat("end:", end, "\n")
      if(end == "3p"){
        seq = substr(seq, start=1, stop=(nchar(seq)-1))
      } else if(end == "5p"){
        seq = substr(seq, 2, stop=nchar(seq))
      }
      
    }
    if(!(seq %in% isomiR_seqs)){
      isomiR_seqs = c(isomiR_seqs, seq)
      differences = c(differences, diff)
    }
    k=k+1
  }
  
  cat("Generating counts for length variant sequences \n")
  nj = filter(count_df, uniqueSequence == center_seq) %>% select(., count) %>% unlist() %>% unname()
  
  cat("Alignments between center sequence and generated isomiR sequences \n")
  alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern=center_seq)
  names(alignments) = isomiR_seqs
  cat("Getting transitions\n")
  transitions = lapply(alignments, get_transitions)
  cat("Computing error distribution means\n")
  lambdas = lapply(transitions, compute_lambda, transition_probs=transition_probs) %>% unlist()
  error_dist_means = lambdas * nj 
  
  isomiR_sequence_counts = vector(length=length(isomiR_seqs))
  names(isomiR_sequence_counts) = isomiR_seqs
  cat("Drawing counts\n")
  set.seed(seed)
  for(seq in isomiR_seqs){
    isomiR_sequence_counts[seq] = rpois(1, error_dist_means[seq])
    if(add_noise){
      sd = sqrt(error_dist_means[seq])
      noise = sample(0:sd, 1)
      isomiR_sequence_counts[seq] = isomiR_sequence_counts[seq]+noise
    }
  }
  
  cat("Assembling rowdata, countdata objects in appropriate format\n")
  rowdata = data.frame(miRNA_name = rep(miRNA, length(isomiR_seqs)+1))
  rowdata$uniqueSequence = c(center_seq, isomiR_seqs)
  
  counts = c(nj, isomiR_sequence_counts)
  
  return(list(center_seq=center_seq, differences=differences, isomiR_seqs=isomiR_seqs, isomiR_sequence_counts=isomiR_sequence_counts, rowdata=rowdata, counts=counts))
  
}