#function that takes simulated dataset true positives and inferred partitioning and returns
#performance metrics

library(tidyverse)
library(Biostrings)

compute_performance_metrics = function(inferred_partition_df, true_center_seqs){
  #compute true positive rate - pull out inferred center sequences
  #any inferred center sequence that is also in true_center_seqs vector is a true positive
  inferred_center_seqs = filter(inferred_partition_df, center == 1) %>% select(., uniqueSequence) %>%
    unlist() %>% unname()
  num_tp = sum(inferred_center_seqs %in% true_center_seqs)
  #any inferred center sequence not in true center sequences vector is a false positive
  #can get this number by subtracting from the total number of positives the number of true positives
  num_fp = length(inferred_center_seqs) - num_tp 
  inferred_technical_isomiRs = inferred_partition_df$uniqueSequence[!(inferred_partition_df$uniqueSequence 
                                                                     %in% inferred_center_seqs)]
  #any inferred technical isomiR that appears in the true center seqs vector is a false negative 
  num_fn  = sum(inferred_technical_isomiRs %in% true_center_seqs)
  num_tn = length(inferred_technical_isomiRs) - num_fn
  
  fp_rate = num_fp/(num_fp + num_tn)
  tp_rate = num_tp/(num_tp + num_fn)
  tn_rate = 1-fp_rate
  fn_rate  = 1-tp_rate
  return(list(num_tp=num_tp, num_fp=num_fp, num_fn=num_fn, num_tn=num_tn, fp_rate=fp_rate,
         tp_rate=tp_rate, tn_rate=tn_rate, fn_rate=fn_rate))
}
