##R script for simulating data with no signal to evaluate false positive rate of our method 

#first clear environment
rm(list=ls())

args =  commandArgs(trailingOnly=TRUE)
idx = as.numeric(args[1])

#load necessary packages
library(tidyverse)
library(Biostrings)

#source functions I've written from appropriate file
source("/scratch/hswan/thesis_isomiR_count_denoising/denoise_isomiR_counts_WORKING_FUNCTION.R")

#load initial transition probs matrix we calculated 
transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

#load mouse data 
mousedata = load_mouse_miRNA_data()

#rowdata
rowdata = mousedata$rowdata
#countdata
countdata = mousedata$countdata
#complete pooling across all samples by summing read counts for a given sequence across samples
count_df=rowSums(countdata)

#get list of unique miRNAs 
unique_miRNAs = unique(rowdata$miRNA_name)

#note - some sequences map to multiple miRNAs after alignment via sRNAbench 
#since these guys aren't true miRNAs, we want to leave them out 
#we can tell which ones are which by looking at miRNA names and applying strsplit with pattern = "Mmu"
#if we get back something that has length > 2 then it isn't a true miRNA 
x = sapply(unique_miRNAs, function(x) return(strsplit(x, "Mmu") %>% unlist() %>% length()-1))
true_miRNAs = names(x[x==1])

#function that takes an error sequence and a center sequence and returns the number of differences between 
#error sequence and center sequence

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

#select a miRNA:
#idx = 1
miRNA = true_miRNAs[idx]
cat("simulating counts from error model for miRNA", miRNA, "\n")

#ID center sequence
center_seq = cbind(rowdata, count=count_df) %>% data.frame() %>% filter(., miRNA_name==miRNA) %>%
  filter(.,count==max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()

#get error sequences 
error_seqs = rowdata %>% filter(., miRNA_name == miRNA & uniqueSequence != center_seq) %>%  
  select(., uniqueSequence) %>% unlist()

num_diff = sapply(error_seqs, count_number_differences, center_seq=center_seq)
names(num_diff) = error_seqs
simple_isomiRs = names(num_diff[num_diff==1])

#get alignments 
alignments = lapply(simple_isomiRs, Biostrings::pairwiseAlignment, pattern=center_seq)
names(alignments) = simple_isomiRs

#get transitions
transitions =  lapply(alignments, get_transitions)

#compute lambdas
lambdas = lapply(transitions, compute_lambda, transition_probs=transition_probs) %>% unlist()

#get nj
nj = cbind(mousedata$rowdata, count=count_df) %>% data.frame() %>% filter(., uniqueSequence == center_seq) %>%
  select(., count) %>% unlist() %>% unname()

#calculate error distribution means
error_dist_means = lambdas*nj

ndatasets = 100
datasets  = list()

set.seed(1989)
for(d in 1:ndatasets){
  x = sapply(error_dist_means, function(x) rpois(1, x))
  sd = sapply(error_dist_means, sqrt)
  noise = sapply(sd, function(x) return(sample(-x:x, 1)))
  simulated = round(x+noise,0)
  datasets[[d]] = c(simulated, nj)
}

simulated_rowdata = rbind(filter(mousedata$rowdata, uniqueSequence %in% simple_isomiRs), 
                          filter(mousedata$rowdata, uniqueSequence==center_seq))

partition_objs = list()
start = Sys.time()
for(d in 1:length(datasets)){
  partition_objs[[d]] = denoise_isomiR_counts(simulated_rowdata, datasets[[d]], transition_probs, miRNA, 0.05,
                                              10, "BH")
}
end = Sys.time()-start
cat("Run time:", end, "\n")
num_partitions = lapply(partition_objs, function(x) return(x$partition_df$partition %>% max())) %>% 
  unlist()

fp_rates = num_partitions-1/nrow(rowdata)
filename  = paste0(miRNA, "false_positive_rates.Rds")
saveRDS(fp_rates, paste0("/scratch/hswan/thesis_isomiR_count_denoising/simulate_null_data/", 
                         filename))
