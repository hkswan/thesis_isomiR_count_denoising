rm(list = ls())
source("/scratch/hswan/thesis_isomiR_count_denoising/Code/denoise_isomiR_counts_WORKING_FUNCTION.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/Code/load_mouse_miRNA_data_function.R")


simulate_length_variant_isomiRs = function(rowdata, countdata, miRNA, add_noise = FALSE, num_true_partitions = 1, num_seqs_per_grp = 100, 
                                           max_iter = 1000, transition_probs, seed = 1989){
  count_df = data.frame(rowdata, count = countdata) 
  count_df = filter(count_df, miRNA_name == miRNA) %>% arrange(., desc(count))
  center_seqs = count_df$uniqueSequence[1:num_true_partitions]
  
  ends = c("3p", "5p")
  nucleotides = c("A", "C", "G", "T")
  
  isomiR_seqs = list()
  set.seed(seed)
  for(seq in center_seqs){
    seqs = vector()
    i = 0 
    cat("Generating isomiR sequences from center sequence", seq, "\n")
    while(i < max_iter & length(seqs) < num_seqs_per_grp){
      #sample difference
      x = seq 
      num_differences = sample(-7:7, 1)
      if(num_differences < 0){
        for(j in 1:-num_differences){
          #1st sample end:
          e = sample(ends, 1)
          #then sample nt: 
          nt = sample(nucleotides, 1)
          if(e == "3p"){
            x = paste0(x, nt, collapse = "")
          } else if(e == "5p"){
            x = paste0(nt, x, collapse="")
          }
        }
      } else if(num_differences > 0){
        for(j in 1:num_differences){
          e = sample(ends, 1)
          if(e == "3p"){
            x = substr(x, start = 1,  stop = nchar(x)-1)
          } else if(e == "5p"){
            x = substr(x, start = 2, stop = nchar(x))
          }
        }
      }
      if(!(x %in% seqs) & !(x %in% center_seqs)){
        seqs = c(seqs, x)
      }
      i = i + 1 
    }
    isomiR_seqs[[seq]] = seqs 
  }
  names(isomiR_seqs) = center_seqs
  
  #rmv any duplicate sequences
  iso = unique(unlist(isomiR_seqs))
  duplicate_seqs = vector()
  for(seq in iso){
    membership = lapply(isomiR_seqs, function(x) seq %in% x) %>% unlist()
    if(sum(membership) > 1){
      duplicate_seqs = c(duplicate_seqs, seq)
    }
  }
  
  unique_isomiR_seqs = lapply(isomiR_seqs, function(x) return(x[!(x %in% duplicate_seqs)]))
  
  cat("Getting alignments \n")
  alignments = list()
  for(seq in center_seqs){
    s = unique_isomiR_seqs[[seq]]
    a = lapply(s, Biostrings::pairwiseAlignment, pattern = seq)
    alignments[[seq]] = a 
  }
  cat("Getting transitions \n")
  
  transitions = list()
  for(seq in center_seqs){
    a = alignments [[seq]]
    t = lapply(a, get_transitions)
    transitions[[seq]] = t
  }
  
  cat("Calculating lambdas \n")
  lambdas = list()
  for(seq in center_seqs){
    t = transitions[[seq]]
    l = lapply(t, compute_lambda, transition_probs = transition_probs)
    lambdas[[seq]] = l
  }
  
  cat("Calculating means of error distributions \n")
  mus = list()
  for(seq in center_seqs){
    l = lambdas[[seq]] %>% unlist()
    nj = filter(count_df, uniqueSequence == seq) %>% select(., count) %>% unlist() %>% unname()
    mu = l*nj
    mus[[seq]] = mu
  }

  cat("Drawing counts \n")
  counts = vector()
  set.seed(seed)
  for(seq in center_seqs){
    mu = mus[[seq]]
    for(m in mu){
      counts = c(counts, rpois(1, m))
    }
  }
  
  cat("Packing simulated data for use with isomiR count denoising functions \n")
  sim_rowdata = data.frame(c(unlist(unique_isomiR_seqs), center_seqs))
  colnames(sim_rowdata) = "uniqueSequence"
  sim_rowdata$miRNA_name = rep(miRNA, nrow(sim_rowdata))
  
  for(seq in center_seqs){
    counts = c(counts, filter(count_df, uniqueSequence == seq) %>% select(., count) %>% unlist() %>% unname())
  }
  
  
  return(list(center_seqs=center_seqs, isomiR_seqs=isomiR_seqs, duplicate_seqs=duplicate_seqs, unique_isomiR_seqs = unique_isomiR_seqs, 
              mus = mus, counts=counts, sim_rowdata = sim_rowdata))
}

# tst_datasets = list()
# 
# for(k in 1:5){
#   cat("Generating dataset", k, "\n")
#   tst_datasets[[k]] = simulate_length_variant_isomiRs(rowdata, countdata, miRNA, FALSE, 3, 100, 1000, transition_probs, k)
# }

# sim_rowdatas = lapply(tst_datasets, function(x) return(x[['sim_rowdata']]))
# sim_counts = lapply(tst_datasets, function(x) return(x[['counts']]))
# 
# partition_objs = list()
# for(k in 1:5){
#   cat("Dataset:", k, "\n")
#   partition_objs[[k]] = denoise_isomiR_counts(tst_rowdatas[[k]], tst_counts[[k]], transition_probs, miRNA, 0.05, 10, "BH")
# }
