rm(list = ls())

library(tidyverse)

matrix_files = list.files("/scratch/hswan/thesis_isomiR_count_denoising/transition_count_matrices")
source("/scratch/hswan/thesis_isomiR_count_denoising/load_mouse_miRNA_data_function.R")
source("/scratch/hswan/thesis_isomiR_count_denoising/correct_technical_length_variant_functions.R")
setwd("/scratch/hswan/thesis_isomiR_count_denoising/transition_count_matrices")

matrices = lapply(matrix_files, readRDS)

correct_NaN = function(matrix){
  for(i in 1:nrow(matrix)){
    for(j in 1:nrow(matrix)){
      if(is.nan(matrix[i,j])){
        matrix[i,j]=0
      }
    }
  }
  return(matrix)
}

matrices = lapply(matrices, correct_NaN)

avg_mat = matrix(0, nrow(matrices[[1]]), ncol(matrices[[1]]))
row.names(avg_mat) = row.names(matrices[[1]])
colnames(avg_mat) = colnames(matrices[[1]])

for(i in 1:length(matrices)){
  avg_mat = avg_mat+matrices[[i]]
}

for(i in 1:nrow(avg_mat)){
  avg_mat[i,] = avg_mat[i,]/sum(avg_mat[i,])
}
avg_mat = avg_mat/length(matrices)

transition_probs = avg_mat
cat("Transition probabilities are\n")
print(transition_probs)

test_miRNA = "Mmu-Mir-362-P3_3p"

mouse_data = load_mouse_miRNA_data()
rowdata = mouse_data$rowdata
countdata = mouse_data$countdata

#make initial partition:
partition_df = cbind(rowdata, partition=rep(1,nrow(rowdata))) %>% filter(., miRNA_name == test_miRNA)
partition_df$center = rep(0, nrow(partition_df))

#use rowsums to identify center sequence for initialized partition 
count_df = rowSums(countdata)
count_df = cbind(rowdata, count=count_df) %>% data.frame()


initial_center_seq = filter(count_df, uniqueSequence %in% partition_df$uniqueSequence) %>% filter(., count == max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()

partition_df$center[partition_df$uniqueSequence == initial_center_seq] = 1 
head(partition_df)

isomiR_seqs = partition_df$uniqueSequence[partition_df$uniqueSequence != initial_center_seq]

alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = initial_center_seq)
names(alignments) = isomiR_seqs
transitions = lapply(alignments, get_transitions)
lambdas = lapply(transitions, compute_lambda, transition_probs=transition_probs)

#computing abundance p-values 
countdata = cbind(rowdata, countdata)

n_js_by_sample = filter(countdata, uniqueSequence == initial_center_seq) %>% select(., colnames(countdata)[!(colnames(.) %in% colnames(rowdata))]) %>%
  unlist()
seq = isomiR_seqs[1]
n_is_by_sample = filter(countdata, uniqueSequence == seq) %>% select(., colnames(countdata)[!(colnames(.) %in% colnames(rowdata))]) %>% unlist()
lambda = lambdas[[seq]]

count_data_cols = colnames(countdata)[!(colnames(countdata) %in% colnames(rowdata))]
num = sapply(count_data_cols, function(x) ppois(n_is_by_sample[[x]], n_js_by_sample[[x]]*lambda, lower.tail=F)) 
denom = sapply(count_data_cols, function(x) return(1-dpois(0, n_js_by_sample[[x]]*lambda)))

##trying a different way of estimating transition probabilities:
transition_probs_matrix = matrix(0, 3, 3)
row.names(transition_probs_matrix) = c("N", "S", "-")
colnames(transition_probs_matrix) = c("N", "S", "-")


for(t in transitions){
  alignment_length = ncol(t)
  for(i in 1:alignment_length){
    trns = get_transition_at_idx(i, t)
    x=trns$x
    y=trns$y
    if(x == y){
      transition_probs_matrix["N", "N"] = transition_probs_matrix["N", "N"]+1
    }
    if(x != y & x != "-" & y != "-"){
      transition_probs_matrix["S", "S"] = transition_probs_matrix["S", "S"] + 1
    }
    if(x == "-" & y != "-"){
      transition_probs_matrix["-", "S"] = transition_probs_matrix["-", "S"]+1
    }
    if(x != "-" & y == "-"){
      transition_probs_matrix["N", "-"] = transition_probs_matrix["N", "-"]+1
    }
  }
}


for(i in 1:nrow(transition_probs_matrix)){
  transition_probs_matrix[i,]=transition_probs_matrix[i,]/sum(transition_probs_matrix[i,])
}


transition = transitions[[1]]
transition_to_type  = function(transition, idx){
  x = transition[1,idx]
  y = transition[2,idx]
  if(x == "-" & y != "-"){
    type = "ins"
  }
  if(x != "-" & y == "-"){
    type = "del"
  }
  if(x != "-"  & y  != "-" & x != y){
    type = "sub"
  }
  if(x == y){
    type = "corr"
  }
  cat("Original nucleotide:", x, "\n")
  cat("Nucleotide in isomiR sequence:",  y, "\n")
  cat("Transition type:", type, "\n")
  return(type)
}

transition_prob_vector = rep(0, 4)
names(transition_prob_vector) = c("ins", "del", "sub", "corr")
for(t in transitions){
  for(i in 1:ncol(t)){
    ttype = transition_to_type(t, i)
    transition_prob_vector[ttype] = transition_prob_vector[ttype]+1
  }
}

transition_prob_vector = transition_prob_vector/sum(transition_prob_vector)

calc_lambda_v2 = function(transition, transition_prob_vector){
  lambda_components = vector()
  for(i in 1:ncol(transition)){
    x = transition[1,i]
    y = transition[2,i]
    if(x == y & x != "-"){
      lambda_components = c(lambda_components, transition_prob_vector["corr"])
    }
    if(x == "-" & y != "-"){
      lambda_components = c(lambda_components, transition_prob_vector["ins"])
    }
    if(x != "-" & y == "-"){
      lambda_components = c(lambda_components, transition_prob_vector["del"])
    }
    if(x != y & x != "-" & y != "-"){
      lambda_components = c(lambda_components, transition_prob_vector["sub"])
    }
  }
  return(lambda_components)
}

lambda_components = calc_lambda_v2(transitions[[seq]], transition_prob_vector)
lambda_2 = prod(lambda_components)

num = ppois(n_is_by_sample[[5]], n_js_by_sample[[5]]*lambda, lower.tail = F) + dpois(n_is_by_sample[[5]], n_js_by_sample[[5]]*lambda)
denom = 1-dpois(0, n_js_by_sample[[5]]*lambda)
p = num / denom
print(p)

n_is_by_sample_2 = filter(countdata, uniqueSequence == seq) %>% select(., colnames(countdata)[(colnames(countdata) %in% count_data_cols)]) %>% unlist()

p_components_num = vector()
p_components_denom = vector()
for(i in 1:length(n_is_by_sample)){
  num = ppois(n_is_by_sample_2[[i]], n_js_by_sample[[i]]*lambda, lower.tail=F)
  denom = 1-dpois(0, n_js_by_sample[[i]]*lambda)
  p_components_num = c(p_components_num, num)
  p_components_denom = c(p_components_denom, denom)
}

## see how it works for a highly expressed isomiR sequence 

seq = df[order(df$count, decreasing=T)[2], 'uniqueSequence']
transition = transitions[[seq]]
lambda_components = calc_lambda_v2(transition)

lambdas_v2 = lapply(transitions, calc_lambda_v2, transition_prob_vector=transition_prob_vector)
lambdas_v2 = lapply(lambdas_v2, prod)


compute_abundance_p_value = function(lambdas, sequence, center_sequence, countdata){
  non_count_cols = c("uniqueSequence", "miRNA_name")
  lambda = lambdas[[sequence]]
  
  n_is_by_sample = filter(countdata, uniqueSequence == sequence) %>% 
    select(., colnames(.)[!(colnames(.) %in% non_count_cols)])  %>% unlist()
  n_js_by_sample = filter(countdata, uniqueSequence == center_sequence) %>% 
    select(., colnames(.)[!(colnames(.) %in% non_count_cols)]) %>% unlist()
  
  p_num_components = vector()
  p_denom_components = vector()
  for(k in 1:length(n_is_by_sample)){
    p_num = ppois(n_is_by_sample[k], n_js_by_sample[k]*lambda, lower.tail=F)
    p_denom = dpois(0, n_js_by_sample[k]*lambda)
    
    p_num_components = c(p_num_components, p_num)
    p_denom_components = c(p_denom_components, p_denom)
  }
  
  p = prod(p_num_components/p_denom_components)
  
  return(list(n_is_by_sample, n_js_by_sample, p_by_sample = p_num_components/p_denom_components, p=p, 
              p_num_components = p_num_components, p_denom_components = p_denom_components)) 
}

n_js_by_sample = filter(countdata, uniqueSequence == initial_center_seq) %>% 
  select(., colnames(.)[!(colnames(.) %in% c("uniqueSequence", "miRNA_name"))]) %>% unlist()
param = lambdas[[isomiR_seqs[1]]] * sum(n_js_by_sample)
n_is_by_sample = filter(countdata, uniqueSequence == isomiR_seqs[1]) %>%
  select(., colnames(.)[!(colnames(.) %in% c("uniqueSequence", "miRNA_name"))]) %>% unlist()

obs_sum = sum(n_is_by_sample)

ppois(obs_sum, param, lower.tail=F)
1-dpois(0, param)


uniq_miRNAs = unique(count_df$miRNA_name)
proportion_of_samples_0_counts = vector(length=length(isomiR_seqs))
total_obs_counts = vector(length = length(isomiR_seqs))

count_cols = colnames(countdata)[3:ncol(countdata)]
for(i in 1:length(proportion_of_samples_0_counts)){
  c = filter(countdata, uniqueSequence == isomiR_seqs[i]) %>% select(., count_cols) %>% unlist()
  proportion_of_samples_0_counts[i]=sum(c == 0)/length(c)
  total_obs_counts[i] = sum(c)
}

plot(log(total_obs_counts), proportion_of_samples_0_counts)

ggplot2::ggplot(data.frame(x=log(total_obs_counts), y=proportion_of_samples_0_counts), aes(x,y)) +
  xlab("Log of total observed read count") + ylab("Proportion of samples with 0 read counts") + geom_point() + ggplot2::theme_classic()


## Make a similar plot for multiple miRNAs 
num_isomiRs = sapply(uniq_miRNAs, function(x) return(filter(count_df, miRNA_name == x) %>% unlist() %>% length()))
for(i in uniq_miRNAs[2:25]){
  cat(i, "\n")
  isomiR_seqs = filter(countdata, miRNA_name==i) %>% select(., uniqueSequence) %>% unlist()
  proportion_of_samples_0_counts = vector(length=num_isomiRs[i])
  total_obs_counts = vector(length=length(isomiR_seqs))
  for(j in 1:length(proportion_of_samples_0_counts)){
    c = filter(countdata, uniqueSequence == isomiR_seqs[j]) %>% select(., count_cols) %>% unlist()
    proportion_of_samples_0_counts[j] = sum(c == 0)/length(c)
    total_obs_counts[j]= sum(c)
  }
  p = ggplot2::ggplot(data.frame(x=log(total_obs_counts), y=proportion_of_samples_0_counts), aes(x,y)) + 
    xlab("Log of total observed read  count") + ylab("Proportion of samples with 0 read counts") + geom_point() +
    ggplot2::theme_classic()
  filename =  paste0("/scratch/hswan/thesis_isomiR_count_denoising/plots/", i, "_plot.png")
  ggsave(filename, p)
}


##sum of poisson random variables is poisson distriubted, parameterized by the sum of the parameters 

sum(n_js_by_sample*lambda)
lambda_sum = sum(n_js_by_sample)*lambda
count_sum = sum(n_is_by_sample)
num = ppois(count_sum, lambda_sum, lower.tail=F) + dpois(count_sum, lambda_sum)
denom = 1-dpois(0, lambda_sum)


compute_denom = function(lambda, n_js_by_sample){
  prob_0_reads = sapply(n_js_by_sample, function(x) return(dpois(0, x*lambda)))
  denom = 1-prod(prob_0_reads)
  return(denom)
}

compute_num = function(lambda, n_js_by_sample, n_is_by_sample){
  prob_greater_reads = sapply(1:length(n_is_by_sample), 
                              function(x) return(ppois(n_is_by_sample[x]-1, 
                                                       n_js_by_sample[x]*lambda, lower.tail=F)))
  num = prod(prob_greater_reads)
  return(num)
}

count_cols = colnames(countdata)[!(colnames(countdata) %in% c("uniqueSequence", 'miRNA_name'))]
n_is_by_sample = filter(countdata, uniqueSequence == isomiR_seqs[1]) %>% 
  select(., colnames(.)[colnames(.) %in% count_cols]) %>% unlist()


p_values = vector()
for(i in isomiR_seqs){
  lambda = lambdas[[i]]
  n_is_by_sample = filter(countdata, uniqueSequence == i) %>% 
    select(., colnames(.)[colnames(.) %in% count_cols])  %>% unlist()
  num = compute_num(lambda, n_js_by_sample, n_is_by_sample)
  denom = compute_denom(lambda, n_js_by_sample)
  p_values = c(p_values, num/denom)
}

adj_p_values = p.adjust(p_values, method="BH")
head(adj_p_values)
length(adj_p_values)


#look at the sum 
n_i = sum(n_is_by_sample)


p_values_sum = vector(length = length(isomiR_seqs))
names(p_values_sum) = isomiR_seqs
for(i in isomiR_seqs){
  lambda = lambdas[[i]]
  n_i = filter(count_df, uniqueSequence == i) %>% select(., count) %>% unlist()
  num = ppois(n_i-1, n_j*lambda, lower.tail=F)
  denom = 1-dpois(0, n_j*lambda)
  p=num/denom
  if(length(p) > 1){
    print(i)
  }
  p_values_sum[i] = p
}

adj_p_values_sum = p.adjust(p_values_sum, method="BH")
head(adj_p_values_sum)
sum(adj_p_values_sum < 0.05)
sum(p_values_sum < 0.05/912)

sum_obs = filter(count_df, uniqueSequence %in% isomiR_seqs) %>%  select(., c(uniqueSequence,count)) 
#names(sum_obs) = isomiR_seqs
sig_seq_obs_sums= sum_obs[p_values_sum < 0.05/912,]
new_center_seq = filter(sig_seq_obs_sums, count == max(count)) %>% select(., uniqueSequence) %>% unlist()

update_df = partition_df
update_df$center[update_df$uniqueSequence == new_center_seq] = 1
update_df$partition[update_df$uniqueSequence == new_center_seq] = 2

master_alignments = list()
master_alignments[['1']] = alignments

master_lambdas = list()
master_lambdas[['1']] = lambdas

isomiR_seqs = isomiR_seqs[isomiR_seqs != new_center_seq]
new_alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = new_center_seq)
new_transitions = lapply(new_alignments, get_transitions)
names(new_transitions) = isomiR_seqs
new_lambdas = lapply(new_transitions, compute_lambda, transition_probs=transition_probs)
master_lambdas[['2']] = new_lambdas

for(i in isomiR_seqs){
  print(i)
  isomiR_lambdas = lapply(master_lambdas, function(x) return(x[[i]])) %>% unlist()
  p = names(isomiR_lambdas)[isomiR_lambdas == max(isomiR_lambdas)]
  print(isomiR_lambdas)
  #p=p[1]
  update_df$partition[update_df$uniqueSequence == i] = p
}

partition_df = update_df
uniq_partitions = unique(partition_df$partition)

p_values = vector(length = length(isomiR_seqs))
names(p_values) = isomiR_seqs
for(j in uniq_partitions){
  center_seq = filter(partition_df, partition == j & center == 1) %>% select(., uniqueSequence) %>%
    unlist()
  n_j = filter(count_df, uniqueSequence == center_seq) %>% select(., count) %>% unlist()
  partition_isomiRs  = filter(partition_df, partition == j & center == 0) %>% select(., uniqueSequence) %>%
    unlist()
  cat("Partition has", length(partition_isomiRs), "isomiRs\n")
  
  #print(partition_isomiRs)
   for(i in partition_isomiRs){
     if(is.na(i)){
       print("is na")
     }
     else{
       print(p_values[i])
     }
     n_i = filter(count_df, uniqueSequence == i) %>% select(., count) %>% unlist()
     cat(i, "\n")
     cat("n_i:", n_i, "\n")
     lambda = master_lambdas[[j]][[i]]
     cat("lambda:", lambda, "\n")
    num = ppois(n_i-1, n_j*lambda, lower.tail=F)
    denom  = 1-dpois(0, n_j*lambda)
    p_values[i] = num/denom

   }
}

new_center_seq = filter(count_df, uniqueSequence %in% names(p_values)[p_values < (0.05/length(p_values))]) %>% 
  filter(., count == max(count))  %>% select(., uniqueSequence) %>% unlist()


miRNAs_no_dupe_isomiRs = vector()
for(miRNA in unique(rowdata$miRNA_name)){
  num_rows = filter(rowdata, miRNA_name  == miRNA) %>% nrow()
  num_unique_isomiRs = filter(rowdata, miRNA_name == miRNA) %>% select(., uniqueSequence) %>%
    unlist() %>% unique() %>% length()
  if(num_rows == num_unique_isomiRs){
    miRNAs_no_dupe_isomiRs = c(miRNAs_no_dupe_isomiRs, miRNA)
  }
}

##picking a new test miRNA because there are duplicate isomiR sequences mapping to some miRNAs 

new_test_miRNA = miRNAs_no_dupe_isomiRs[1]

#create initial partition object
partition_df = filter(rowdata, miRNA_name == new_test_miRNA)
#create center indicator variable
partition_df$center = rep(0, nrow(partition_df))
#choose initial center sequence
initial_center_seq = filter(count_df,  miRNA_name == new_test_miRNA) %>% filter(., count == max(count)) %>%
  select(., uniqueSequence) %>% unlist()
partition_df$center[partition_df$uniqueSequence == initial_center_seq] = 1
#everyone starts in the same partition
partition_df$partition = rep(1, nrow(partition_df))
head(partition_df)
