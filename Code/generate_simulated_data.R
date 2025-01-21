#R script that makes simulated data from sequences appearing in ernesto's benchmark dataset #
#simulated count data is for isomiRs mapping to a user-specified miRNA

rm(list=ls())

args = commandArgs(trailingOnly = TRUE)
ndatasets = 500
miRNA_idx = as.numeric(args[1])



source("/scratch/hswan/thesis_isomiR_count_denoising/denoise_isomiR_counts_WORKING_FUNCTION.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/function_dev/count_number_differences_FUNCTION.R")
#load data
mousedata = load_mouse_miRNA_data()

rowdata = mousedata$rowdata
countdata = mousedata$countdata

count_df = rowSums(countdata)

transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

unique_miRNAs = unique(rowdata$miRNA_name)
true_miRNAs = sapply(unique_miRNAs, function(x) strsplit(x,"Mmu") %>% 
                              unlist() %>% length()-1)
true_miRNAs = unique_miRNAs[true_miRNAs==1]

miRNA = true_miRNAs[miRNA_idx]

cat("Generating",  ndatasets, "null datasets from isomiR-level count data of miRNA", miRNA, "\n")

#make a psuedo partition_df object:
df = cbind(rowdata, count=count_df) %>% filter(., miRNA_name == miRNA)

#id a center sequence for user-specified miRNA
cat("Identifying center sequence for miRNA", miRNA,  "\n")
center_seq = filter(df, count == max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()

#all other sequences mapping to miRNA are error_seqs:
cat("All other sequences are error sequences under null hypothesis \n")
error_seqs = df$uniqueSequence[df$uniqueSequence != center_seq]

#get number of differences between each error sequence and the center sequence we identified
cat("There are", length(error_seqs), "error sequences \n")
cat("Getting number of differences between each error sequence and center sequence we identified\n")
start = Sys.time()
num_diffs = sapply(error_seqs, count_number_differences, center_seq=center_seq)
end = Sys.time()-start
cat("Getting number of differences between center sequence and error sequences took", end, "\n")
#names(num_diffs) = error_seqs 

#get simple isomiRs 
simple_isomiRs = names(num_diffs[num_diffs == 1])

#check that there is at least 1 simple isomiR in the dataset: 

num_simple_isomiRs = length(simple_isomiRs)

if(num_simple_isomiRs !=  0){
  cat("Aligning simple isomiRs to center sequence\n")
  alignments = lapply(simple_isomiRs, Biostrings::pairwiseAlignment, pattern=center_seq)
  names(alignments) = simple_isomiRs
  cat("Getting transition objects from alignments \n")
  transitions = lapply(alignments, get_transitions)
  cat("Calculating lambdas")
  lambdas  = lapply(transitions, compute_lambda, transition_probs=transition_probs) %>% unlist()
  nj = filter(df, uniqueSequence == center_seq) %>% select(., count) %>% unlist() %>% unname()
  cat("Observed read count of center sequence is", nj, "reads.\n")
  cat("Calculating means of error distributions by multiplying the lambda vector by  nj")
  error_dist_means = nj*lambdas
  
  #generate list of random seeds:
  seeds = 1:ndatasets
  
  datasets = list()
  #draw counts from error distributions:
  for(i in 1:ndatasets){
    set.seed(seeds[i])
    x = sapply(error_dist_means, function(y) rpois(1, y))
    sd = sapply(error_dist_means, sqrt)
    noise = sapply(sd, function(y) sample(-y:y, 1))
    simulated_counts = c(x+noise, nj)
    #round to integer
    simulated_counts=round(simulated_counts, 0)
    simulated_rowdata = data.frame(uniqueSequence = names(error_dist_means), miRNA_name = rep(miRNA, length(error_dist_means)))
    simulated_rowdata = rbind(simulated_rowdata, c(center_seq, miRNA))
    datasets[[i]] = list(simulated_counts=simulated_counts, simulated_rowdata=simulated_rowdata)
  }
  
  #datasets[[length(datasets)+1]] = seeds
  
  filepath = "/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/updated_simulated_data/"
  filename = paste0(miRNA, "_", collapse="") %>% paste0(., ndatasets, collapse="") %>% paste0(., "_simulated_datasets.Rds")
  
  saveRDS(datasets, paste0(filepath, filename, collapse=""))
} else{
  cat("No simple isomiRs exist in data. Will need to generate them.\n")
}



