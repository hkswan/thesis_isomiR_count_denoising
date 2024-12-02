##R script for partitioning miRNA sequences into groups based on which isomiR sequence produced that sequence 
#using updated formula for calculating abundance p-values, also looking at the sum of the independent Poisson r.v.'s rather than 
#looking at them individually and then taking the product

#clear environment 
rm(list = ls())

#load necessary packages
library(Biostrings)
library(tidyverse)
library(SummarizedExperiment)

#load some other R scripts that contain some functions that we're going to use
#other functions we need will be written in this script 

source("/scratch/hswan/thesis_isomiR_count_denoising/load_mouse_miRNA_data_function.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/correct_technical_length_variant_functions.R")

#use function to load Ernesto's benchmark dataset to do just that: 
data = load_mouse_miRNA_data()

#list files in transition_count_matrices subdirectory of thesis_isomiR_count_denoising directory to estimate transition probabilities 
matrix_files = list.files("/scratch/hswan/thesis_isomiR_count_denoising/transition_count_matrices")

#load matrices 
matrices = lapply(matrix_files, readRDS)

transition_count_matrix = matrices[[1]]
print(transition_count_matrix)

for(i in 2:length(matrices)){
  transition_count_matrix = transition_count_matrix + matrices[[i]]
}

print(transition_count_matrix)

transition_prob_matrix = transition_count_matrix 
for(i in 1:nrow(transition_prob_matrix)){
  transition_prob_matrix[i,] = transition_prob_matrix[i,]/sum(transition_prob_matrix[i,])  
}

calc_lambda_v2 = function(transition, transition_probs){
  alignment_length = ncol(transition)
  lambda_components  = vector()
  for(i in 1:alignment_length){
    trns = get_transition_at_idx(i, transition)
    x=trns$x
    y=trns$y
    if(x==y){
      lambda_components = c(lambda_components, transition_probs[["same_prob"]])
    }
    if(x != y & x != "-" & y != "-"){
    lambda_components = c(lambda_components, transition_probs[["sub_prob"]])
    }
    if(x == "-" & y != "-"){
      lambda_components = c(lambda_components, transition_probs[["ins_prob"]])
    }
  if(x != "-" & y =="-"){
    lambda_components = c(lambda_components, transition_probs[["del_prob"]])
  }
  }
  p = sum(log(lambda_components)) %>% exp()
  return(p)
}

cat("transition probabilities are:\n")
print(transition_prob_matrix)

##change how we are calculating transition probabilities
A_sub = sum(transition_count_matrix["A", c("C", "G", "T")])
p_A = sum(transition_count_matrix["A",])
C_sub = sum(transition_count_matrix["C", c("A", "G", "T")])
G_sub = sum(transition_count_matrix["G", c("A", "C", "T")])
T_sub = sum(transition_count_matrix["T", c("A", "C", "G")])

#isolate rowdata, countdata 
rowdata = data$rowdata
countdata = data$countdata

#create count_df object that we will be using by calculating rowsums of countdata. going to be looking at the event that the 
#sum of the independent poisson r.v.'s is as extreme or more extreme than the observed sum of read counts across samples (complete pooling)
count_df = rowSums(countdata)
count_df = cbind(rowdata, count=count_df) %>% data.frame()
head(count_df)

#some miRNAs have duplicate sequences mapping to them - this is going to cause problems in our partitioning and Ernesto is looking into it
#Doesn't happen with all the miRNAs so we will just pick a miRNA to work with for now where that doesn't happen 
unique_miRNAs = unique(count_df$miRNA_name)
miRNAs_no_dupe_isomiRs = vector()

for(miRNA in unique_miRNAs){
  isomiRs = filter(count_df, miRNA_name == miRNA) 
  num_rows = nrow(isomiRs)
  num_uniq_seqs = length(unique(isomiRs$uniqueSequence))
  if(num_rows == num_uniq_seqs){
    miRNAs_no_dupe_isomiRs = c(miRNAs_no_dupe_isomiRs, miRNA)
  }
}

test_miRNA = miRNAs_no_dupe_isomiRs[1]
cat("Test miRNA is", test_miRNA, "\n")

cat("Create initial partition\n")
partition_df = filter(count_df, miRNA_name == test_miRNA)
#everyone starts in same partition
partition_df$partition = rep(1, nrow(partition_df))
#id initial center sequence 
initial_center_seq = filter(partition_df, count == max(count)) %>% select(., uniqueSequence) %>% unlist()
cat("Initial center sequence is", initial_center_seq, "\n")
partition_df$center = rep(0, nrow(partition_df))
partition_df$center[partition_df$uniqueSequence == initial_center_seq] = 1
head(partition_df)
filter(partition_df, center == 1)

cat("Getting alignments between isomiR sequences and center sequence \n")
isomiR_seqs = filter(partition_df, center ==0) %>% select(., uniqueSequence) %>% unlist()
alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = initial_center_seq)
names(alignments) = isomiR_seqs

cat("Getting transitions from alignments \n")
transitions = lapply(alignments, get_transitions)

cat("Calculating lambdas\n")
lambdas = lapply(transitions, compute_lambda, transition_probs = transition_prob_matrix)

#ok so this is where a little innovation comes in - our formula for abundance p-value has changed slightly because we're looking at the sum
#and also because we've thought about it a little more carefully 
calculate_abundance_p_value = function(n_j, n_i, lambda){
  num = ppois(n_i-1, n_j*lambda, lower.tail=F)
  denom = 1-dpois(0, n_j*lambda)
  return(num/denom)
}

n_is = sapply(isomiR_seqs, function(x) filter(partition_df, uniqueSequence ==x) %>% select(., count) %>% unlist())
names(n_is) = isomiR_seqs
n_j = filter(partition_df, uniqueSequence == initial_center_seq)%>%select(., count) %>% unlist()

p_values = vector(length = length(isomiR_seqs))
names(p_values) = isomiR_seqs
for(i in 1:length(isomiR_seqs)){
  p = calculate_abundance_p_value(n_j, n_is[isomiR_seqs[i]], lambdas[[isomiR_seqs[i]]])
  p_values[isomiR_seqs[i]] = p
}

#hypothesis testing also has a little bit of innovation - going to use Bonferroni to control false positive rate 
hypothesis_testing = function(p_values, omega_A){
  threshold = omega_A/length(p_values)
  results = data.frame(p_values)
  results$seq = names(p_values)
  results$significant = rep(0, nrow(results))
  results$significant[results$p_values < threshold] = 1
  return(results)
}

results = hypothesis_testing(p_values, 0.05)

#id a new center sequence 
significant_seqs = filter(results, significant == 1) %>% select(., seq) %>% unlist()
new_center_seq = filter(partition_df, uniqueSequence %in% significant_seqs) %>% filter(., count == max(count)) %>% 
  select(., uniqueSequence) %>% unlist()

#copy of partition for updating 
update_df = partition_df
update_df$partition[update_df$uniqueSequence == new_center_seq] = max(update_df$partition) + 1
update_df$center[update_df$uniqueSequence == new_center_seq] = 1

#check to make sure that worked
filter(update_df, center == 1)

#ok now update list of isomiR sequences 
isomiR_seqs = update_df$uniqueSequence[update_df$center != 1]

#align all the isomiR sequences to the new center sequence 
new_alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = new_center_seq)
names(new_alignments) = isomiR_seqs
#get new transitions
new_transitions = lapply(new_alignments, get_transitions)
#get new lambdas - probability of new center sequence being misread as one of isomiR sequences so we can determine which partition 
#they should be in 
new_lambdas = lapply(new_transitions, compute_lambda, transition_probs = transition_prob_matrix)

#now we need a master_lambdas list so we have both old lambdas and these ones so we can determine which center seq is more likely to 
#produce read of given sequence 
master_lambdas = list(lambdas, new_lambdas)
names(master_lambdas) = c("1", "2")

#function that takes update_df partition, master_lambdas list, isomiR sequence, and determines which partition that sequence should join
determine_partition_membership = function(update_df, master_lambdas, sequence){
  lambdas = lapply(master_lambdas, function(x) return(x[[sequence]])) %>% unlist()
  partition = names(lambdas)[lambdas == max(lambdas)]
  update_df$partition[update_df$uniqueSequence == sequence] = as.numeric(partition)
  return(update_df)
}

for(i in isomiR_seqs){
  update_df = determine_partition_membership(update_df, master_lambdas, i)
}


table(update_df$partition)
partition_df = update_df

#check to see proportion of sequences that moved to a new partition:
num_seqs_moved = 0
for(i in partition_df$uniqueSequence){
  if(partition_df$partition[partition_df$uniqueSequence ==i] != update_df$partition[update_df$uniqueSequence==i]){
    num_seqs_moved = num_seqs_moved + 1 
  }
}

proportion_changed = num_seqs_moved/nrow(partition_df)
print(proportion_changed)

p_values = vector(length = length(isomiR_seqs))
names(p_values) = isomiR_seqs

partition_df = update_df
unique_partitions = unique(partition_df$partition)
for(j in unique_partitions){
  center_seq_count = filter(partition_df, partition == j & center == 1) %>% select(., count) %>% unlist(.)
  cat("j:", j, "\n")
  cat("n_j:", center_seq_count, "\n")
  partition_isomiRs = filter(partition_df, partition == j & center == 0) %>% select(., uniqueSequence) %>% unlist()
  for(i in partition_isomiRs){
    n_i = filter(partition_df, uniqueSequence == i) %>% select(., count) %>% unlist()
    lambda = master_lambdas[[j]][[i]]
    p_values[i] = calculate_abundance_p_value(center_seq_count, n_i, lambda)
  }
}

results = hypothesis_testing(p_values, 0.05)
head(results)

#create copy for updating:
update_df = partition_df

for(j in unique_partitions){
  new_partition = max(as.numeric(update_df$partition)) + 1
  new_center_seq = filter(partition_df, center == 0 & partition == j) %>% 
    filter(., uniqueSequence %in% results$seq[results$significant]) %>% filter(., count == max(count)) %>% select(., uniqueSequence) %>%
    unlist()
  update_df$center[update_df$uniqueSequence == new_center_seq] = 1
  update_df$partition[update_df$uniqueSequence == new_center_seq] = new_partition
}

#make sure it worked correctly 
table(update_df$partition)

partitions = unique(update_df$partition)
isomiR_seqs = filter(update_df, center == 0) %>% select(., uniqueSequence) %>% unlist()

for(j in partitions){
  if(!(j %in% names(master_lambdas))){
    center = filter(update_df, partition == j & center == 1) %>% select(., uniqueSequence) %>% unlist()
    new_alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = center)
    names(new_alignments) = isomiR_seqs
    new_transitions = lapply(new_alignments, get_transitions)
    new_lambdas = lapply(new_transitions, compute_lambda, transition_probs = transition_prob_matrix)
    master_lambdas[[j]] = new_lambdas
  }
}

for(seq in isomiR_seqs){
  update_df = determine_partition_membership(update_df, master_lambdas, seq)
}


table(update_df$partition)

partition_df = update_df
unique_partitions = unique(partition_df$partition)

p_values = vector(length = length(isomiR_seqs))
names(p_values) = isomiR_seqs
for(j in unique_partitions){
  cat("j:", j, "\n")
  n_j = filter(partition_df, center == 1 & partition == j) %>% select(., count) %>% unlist()
  cat("n_j:", n_j, "\n")
  partition_isomiRs = filter(partition_df, center == 0 & partition == j) %>% select(., uniqueSequence) %>% unlist()
  for(i in partition_isomiRs){
    lambda = master_lambdas[[j]][[i]]
    n_i = filter(partition_df, uniqueSequence == i) %>% select(., count) %>% unlist()
    p_values[i] = calculate_abundance_p_value(n_j, n_i, lambda)
  }
}

results = hypothesis_testing(p_values, 0.05)
head(results)

##now put everything in a while loop: 
n_iter = 0 
max_iter = 10
epsilon = 0.05
prop_seq_moving = 1
#initialize partitioning:
partition_df = filter(count_df, miRNA_name == test_miRNA) %>% data.frame()
partition_df$partition = rep(1, nrow(partition_df))
partition_df$center = rep(0, nrow(partition_df))
partition_df$center[partition_df$uniqueSequence == initial_center_seq] = 1

filter(partition_df, center == 1)

isomiR_seqs = partition_df$uniqueSequence[partition_df$center != 1]
#get alignments
alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern =  initial_center_seq)
names(alignments) = isomiR_seqs
#get transitions
transitions = lapply(alignments, get_transitions)
#calculate lambdas
lambdas  = lapply(transitions, compute_lambda, transition_probs=transition_prob_matrix)
master_lambdas = list(lambdas)
names(master_lambdas) = unique(partition_df$partition)
#calculate p-values
while(n_iter < max_iter & prop_seq_moving >= epsilon){
  #initialize  p_values vector
  cat("iteration:", n_iter, "\n")
  p_values = vector(length=length(isomiR_seqs))
  names(p_values) = isomiR_seqs
  #get unique partitions
  unique_partitions = unique(partition_df$partition)
  cat("unique partitions:", unique_partitions, "\n")
  for(j in unique_partitions){
    n_j = partition_df$count[partition_df$partition == j &  partition_df$center == 1]
    partition_isomiRs = filter(partition_df, partition == j & center == 0) %>% select(., uniqueSequence) %>%
      unlist()
    for(i in partition_isomiRs){
      n_i = partition_df$count[partition_df$uniqueSequence == i]
      lambda = master_lambdas[[j]][[i]]
      p_values[i] = calculate_abundance_p_value(n_j, n_i, lambda)
    }
  }
  cat("Getting sequences with significant p-values \n")
  significant_seqs = names(p_values)[p_values < 0.05/length(p_values)]
  #for each partition determine if we should start a new one
  #first make a copy for updating
  update_df = partition_df
  for(j in unique_partitions){
    new_partition = max(update_df$partition) + 1
    new_center_seq = filter(partition_df, partition == j  & center == 0) %>% 
      filter(., uniqueSequence %in% significant_seqs) %>% filter(., count == max(count)) %>% 
      select(., uniqueSequence) %>% unlist()
    update_df$partition[update_df$uniqueSequence == new_center_seq]  = new_partition
    update_df$center[update_df$uniqueSequence == new_center_seq] = 1
  }
  print(filter(update_df, center == 1))
  #update isomiR_seqs and unique partitions
  cat("Updating list of isomiR sequences and unique partitions \n")
  isomiR_seqs = update_df$uniqueSequence[update_df$center == 0]
  unique_partitions = unique(update_df$partition)
  # # 
  cat("Getting pairwise alignments between newly identified sequences and new center sequences\n")
  for(j in unique_partitions){
    if(!(j %in% names(master_lambdas))){
      center = filter(update_df, partition == j & center == 1) %>% select(., uniqueSequence)  %>%
        unlist()
      new_alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = center)
      names(new_alignments) = isomiR_seqs
      new_transitions = lapply(new_alignments, get_transitions)
      new_lambdas = lapply(new_transitions, compute_lambda, transition_probs = transition_prob_matrix)
      master_lambdas[[j]] = new_lambdas
    }
  }
  names(master_lambdas) = unique_partitions
  cat("Updating partition membership\n")
  for(i in isomiR_seqs){
    update_df = determine_partition_membership(update_df,master_lambdas, i)
  }
  #get proportion of sequences that moved
  cat("Calculating proportion of sequences that moved to a new proportion \n")
  num_seqs_moved = 0
  for(seq in isomiR_seqs){
    if(update_df$partition[update_df$uniqueSequence == seq] != partition_df$partition[partition_df$uniqueSequence == seq]){
      num_seqs_moved = num_seqs_moved + 1
    }
  }
  prop_seq_moving = num_seqs_moved / nrow(partition_df)
  cat("Proportion of sequences that moved:", prop_seq_moving, "\n")
  #print(prop_seq_moving)
   n_iter = n_iter + 1
  partition_df=update_df
}

#get unique partitions
unique_partitions = unique(partition_df$partition)
for(j in unique_partitions){
  n_j = filter(count_df, partition == j & center == 1) %>% select(., count) %>% unlist()
  partition_isomiRs = filter()
}


#updated function for initializing the partitioning because what we have before wasn't for 
#sums and didn't also return initial center sequence
initialize_partitioning = function(rowdata, count_df, miRNA){
  partition_df = count_df %>% data.frame() %>% filter(., miRNA_name == miRNA)
  partition_df$partition = rep(1, nrow(partition_df))
  partition_df$center = rep(0, nrow(partition_df))
  #id initial center sequence
  initial_center_seq = filter(partition_df, count == max(count)) %>% select(., uniqueSequence) %>%
    unlist()
  partition_df$center[partition_df$uniqueSequence == initial_center_seq] = 1
  return(list(partition_df=partition_df,  initial_center_seq = initial_center_seq))
}


partition_0 = initialize_partitioning(rowdata,count_df, test_miRNA)
partition_df = partition_0$partition_df
head(partition_df)
filter(partition_df, center == 1)
initial_center_seq = partition_0$initial_center_seq

isomiR_seqs = partition_df$uniqueSequence[partition_df$center == 0]
sub_str_search = function(center_seq, isomiR_seq){
  #id which sequence is shorter:
  seq_lengths = sapply(c(center_seq, isomiR_seq), nchar)
  seqs = c(center_seq, isomiR_seq)
  #x1 is longer sequence
  x1 = seqs[which(seq_lengths == max(seq_lengths))]
  x2 = seqs[seqs != x1]
  return(strsplit(x1, x2))
}

alignnments = lapply(isomiR_seqs, pairwiseAlignment, pattern = initial_center_seq)
transitions = lapply(alignments, get_transitions)
lambdas  = lapply(transitions, compute_lambda, transition_probs = transition_prob_matrix)
n_j = sum(partition_df$count)
draws = rpois(10000, n_j*lambdas[[1]])

p_values = vector(length=length(isomiR_seqs))
names(p_values) = isomiR_seqs

for(i in isomiR_seqs){
  n_i = partition_df$count[partition_df$uniqueSequence == i]
  lambda = lambdas[[i]]
  num = ppois(n_i-1, lambda*n_j, lower.tail=F)
  denom = 1-dpois(0, lambda*n_j)
  p_values[i] = num/denom
}

for(i in isomiR_seqs[1:15]){
  cat("Observed read count:", partition_df$count[partition_df$uniqueSequence == i], "\n")
  cat("P-value", p_values[i], "\n")
  cat("P-value is significant:", p_values[i] < 0.05/length(p_values), "\n")
  print(n_j)
}

#ID new center sequence: 
#1. Identify the sequence with the smallest associated p-value
#2. If there is a tie, then identify the most abundant sequence and use that one 

potential_center_seqs = names(p_values)[p_values == min(p_values)]
new_center_seq = filter(partition_df, uniqueSequence %in% potential_center_seqs) %>% filter(., count == max(count)) %>%
  select(., uniqueSequence) %>% unlist()
update_df = partition_df
update_df$partition[update_df$uniqueSequence == new_center_seq] = 2
update_df$center[update_df$uniqueSequence == new_center_seq] = 1

isomiR_seqs = update_df$uniqueSequence[update_df$center == 0]
master_lambdas = list(lambdas)
new_alignments = lapply(isomiR_seqs, pairwiseAlignment, pattern = new_center_seq)
names(new_alignments) = isomiR_seqs
new_transitions = lapply(new_alignments, get_transitions)
new_lambdas = lapply(new_transitions, compute_lambda, transition_probs = transition_prob_matrix)
master_lambdas[[2]] = new_lambdas

names(master_lambdas) = unique(update_df$partition)


for(i in isomiR_seqs){
  l = lapply(master_lambdas, function(x) return(x[[i]])) %>% unlist()
  l = which(l == max(l))
  update_df$partition[update_df$uniqueSequence == i] = l
}

#calculate proportion of sequences that moved 
prop_seqs_moving = 0
for(i in isomiR_seqs){
  if(partition_df$partition[partition_df$uniqueSequence == i] != update_df$partition[update_df$uniqueSequence == i]){
    prop_seqs_moving = prop_seqs_moving + 1 
  }
}
prop_seqs_moving = prop_seqs_moving / nrow(partition_df)

#2nd iteration - calculate p-values 
unique_partitions = unique(update_df$partition) %>% sort()
p_values = vector(length = length(isomiR_seqs))
names(p_values) = isomiR_seqs
for(j in unique_partitions){
  n_j = filter(update_df, partition == j) %>% select(., count) %>% unlist() %>% sum()
  cat("n_j:", n_j, "\n")
  partition_isomiRs = filter(update_df, partition == j & center == 0) %>% select(., uniqueSequence) %>% unlist()
  for(i in partition_isomiRs){
    n_i = update_df$count[update_df$uniqueSequence == i]
    lambda = master_lambdas[[j]][[i]]
    num = ppois(n_i-1, n_j*lambda, lower.tail = F)
    denom = 1-dpois(0, n_j*lambda)
    p_values[[i]] = num/denom
  }
}

#create copy for updating 
partition_df = update_df
update_df = partition_df

#for each existing partition, ID a new center sequence if appropriate 
for(j in unique_partitions){
  cat("j:", j, "\n")
  #step 1 - get partition isomiRs
  partition_isomiRs = filter(partition_df, partition == j & center == 0) %>% select(., uniqueSequence) %>% unlist()
  #step 2 - get p_values for those isomiRs 
  partition_p_values = p_values[partition_isomiRs]
  #step 3 - id partition isomiR with smallest p-value
  new_center_seqs = names(partition_p_values)[partition_p_values == min(partition_p_values)]
  #step 4- not always necessary but of the sequences with the smallest p-values, pick sequence with largest read count to be new center sequence
  #helps break ties 
  new_center_seq = filter(update_df, uniqueSequence %in% new_center_seqs) %>% filter(., count == max(count)) %>% select(., uniqueSequence) %>%
    unlist()
  #create new partition
  new_partition = max(update_df$partition) + 1
  cat("Creating new partition", new_partition, "from partition", j, "\n")
  #now update 
  update_df$partition[update_df$uniqueSequence == new_center_seq] = new_partition
  update_df$center[update_df$uniqueSequence == new_center_seq] = 1
}

#update isomiR sequences list:
isomiR_seqs = update_df$uniqueSequence[update_df$center == 0]
#update unique partitions 
unique_partitions = unique(update_df$partition) %>% sort()
for(j in unique_partitions){
  cat("j:", j, "\n")
   if(!(j %in% names(master_lambdas))){
     cat(j, "\n")
     center = filter(update_df, partition == j & center == 1) %>% select(., uniqueSequence) %>% unlist()
     new_alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = center)
     names(new_alignments) = isomiR_seqs
     new_transitions = lapply(new_alignments, get_transitions)
     new_lambdas = lapply(new_transitions, compute_lambda, transition_probs = transition_prob_matrix)
     master_lambdas[[j]] = new_lambdas 
   }
}

names(master_lambdas) = unique_partitions

registerDoParallel(cores=32)
start = Sys.time()
results = foreach(i=new_partitions, .combine=) %dopar% {
  cat("i:", i,  "\n")
  center = filter(update_df, partition == i & center == 1)  %>% select(., uniqueSequence) %>%
    unlist()
  new_alignments = lapply(isomiR_seqs,  Biostrings::pairwiseAlignment, pattern = center)
  names(new_alignments) = isomiR_seqs
  new_transitions  = lapply(new_alignments, get_transitions)
  new_lambdas = lapply(new_transitions, compute_lambda, transition_probs = transition_prob_matrix)
}
elapsed = Sys.time()-start
cat("Elapsed time:", elapsed, "\n")


##try only letting sequences with a significant p-value move between partitions?? 
#maybe not every sequence needs to move, might save computation time 