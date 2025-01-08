count_number_differences = function(center_seq, error_seq){
  #step 1 - alignment
  alignment = Biostrings::pairwiseAlignment(pattern=center_seq, error_seq)
  #step 2 - transition
  transition = get_transitions(alignment)
  #step 3 
  num_diff = 0 
  for(i in 1:ncol(transition)){
    trns = get_transition_at_idx(i, transition)
    x=trns$x
    y=trns$y
    if(x != y){
      num_diff = num_diff+1
    }
  }
  return(num_diff)
}