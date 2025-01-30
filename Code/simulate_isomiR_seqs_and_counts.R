#function that simulates isomiR sequences of any type and read counts drawn from the assumed error distribution
#uses Ernesto's benchmark dataset to select a center sequence for a given miRNA, then generate the isomiR sequences from
#that center sequence 

#to generate counts, align the generated isomiR sequence to the center sequence, get transitions, compute lambda 

source("/scratch/hswan/thesis_isomiR_count_denoising/Code/denoise_isomiR_counts_update_partitions_at_end.R")
transition_probs = readRDS("/scratch/hswan/thesis_isomiR_count_denoising/initial_transition_probs.Rds")

mousedata = load_mouse_miRNA_data()
rowdata = mousedata$rowdata
countdata = mousedata$countdata


generate_simulated_dataset = function(mouse_rowdata, mouse_countdata, miRNA, num_true_partitions=1, 
                                      seed=1989, max_iter=1000, num_isomiRs=10, transition_probs){
  count_df = rowSums(mouse_countdata)
  d = cbind(mouse_rowdata, count=count_df) %>% data.frame() %>% filter(., miRNA_name == miRNA)
  d = arrange(d, desc(count))
  center_seqs = d$uniqueSequence[1:num_true_partitions]
  njs = d$count[1:num_true_partitions]
  names(njs) = center_seqs
  
  isomiR_sequences = list()
  
  #create our list of nucleotides to sample from: 
  nucs = c("A", "C", "G", "T")
  ends = c("3p", "5p")
  typeof_edits = c("seq", "length")
  typeof_length_edit = c("add", "remove")
  num_iter = 0 
  
  cat("Generating isomiR sequences from real data\n")
  set.seed(seed)
  for(j in 1:num_true_partitions){
    #create vector:
    isomiRs_vec = vector()
    while(length(isomiRs_vec) <= num_isomiRs & num_iter < max_iter){
      #starting sequence is center sequence:
      seq = center_seqs[j]
      #sample edit distance:
      edit_dist = sample(1:5, 1)
      #for each of the edits, 
      for(i in 1:edit_dist){
        #we sample the edit type - either length or sequence 
        edit_type = sample(typeof_edits, 1)
        #if its a length edit, then we're going to add or remove a nucleotide 
        if(edit_type == "length"){
          #choose randomly to add or remove 
          typeof_length_change = sample(typeof_length_edit, 1)
          #if we're adding a nucleotide
          if(typeof_length_change == "add"){
            #sample the nucleotide to add
            nuc_to_add = sample(nucs, 1)
            #sample the end we're going to add it to: 
            end_to_add_to = sample(ends, 1)
            if(end_to_add_to == "3p"){
              seq = paste0(nuc_to_add, seq, collapse="")
            } else if(end_to_add_to == "5p"){
              seq = paste0(seq, nuc_to_add, collapse="")
            }
          } else if(typeof_length_change == 'remove'){
            #otherwise we pick an end of the sequence and remove that nucleotide
            end_to_change = sample(ends, 1)
            if(end_to_change == "5p"){
              #if we're changing the 5p end of the sequence, we remove the first nucleotide of the sequence 
              seq = substr(seq, start = 2, stop = nchar(seq))
            } else if(end_to_change == "3p"){
              #if we're changing the 3p end of the sequence, we remove the last nucleotide of the sequence 
              seq = substr(seq, start = 1 , stop = nchar(seq) - 1)
            }
          }
        } else if(edit_type == "seq"){
          #otherwise if the edit type is a sequence then we are going to randomly select a position along the sequence
          #and edit the nucleotide that's currently present at that position 
          nuc_being_subbed = sample(nucs, 1)
          position_to_sub_at = sample(1:nchar(seq), 1)
          substr(seq, position_to_sub_at, position_to_sub_at) = nuc_being_subbed
        }
      }
      #check to see if sequence is already in vector 
      if(!(seq %in% isomiRs_vec) & seq != center_seqs[j]){
        isomiRs_vec = c(isomiRs_vec, seq)
      }
      num_iter = num_iter+1
    }
    isomiR_sequences[[j]] = isomiRs_vec
  }
  names(isomiR_sequences) = center_seqs
  
  unique_isomiRs = unlist(isomiR_sequences) %>% unname() %>% unique()
  
  if(length(unique_isomiRs) != length(unlist(isomiR_sequences))){
    cat("Duplicate sequences present. Need to remove duplicates \n")
    for(seq in unique_isomiRs){
      counter = 0
      idxs = vector()
      for(j in 1:length(isomiR_sequences)){
        if(seq %in% isomiR_sequences[[j]]){
          counter = counter+1
          idxs = c(idxs, j)
        }
      }
      if(counter > 1){
        #randomly select which one its going to stay in 
        keep = sample(idxs, 1)
        remove_frm = !(1:length(isomiR_sequences) %in% keep)
        for(k in remove_frm){
          x = isomiR_sequences[[k]]
          z = which(x == seq)
          x = x[-z]
          isomiR_sequences[[k]] = x
        }
      }
    }
  }
  
  cat("Calculating location parameters of error distributions for generated sequences \n")
  
  lambdas = list()
  
  for(j in 1:length(isomiR_sequences)){
    #alignments
    alignments = lapply(isomiR_sequences[[j]], Biostrings::pairwiseAlignment, pattern = names(isomiR_sequences[j]))
    names(alignments) = isomiR_sequences[[j]]
    #transitions:
    transitions = lapply(alignments, get_transitions)
    #use transition probabilities, calculate lambdas
    lambdas[[j]] = lapply(transitions, compute_lambda, transition_probs=transition_probs) %>% unlist()
  }
  names(lambdas) = center_seqs
  
  mus = list()
  for(j in 1:length(lambdas)){
    mus[[j]] = njs[center_seqs[j]]*lambdas[[j]]
  }
  
  #now draw counts 
  sim_counts = vector()
    for(j in 1:length(mus)){
      mu_vec = mus[[j]]
      for(i in 1:length(mu_vec)){
        sim_counts = c(sim_counts, rpois(1, mu_vec[i]))
      }
    }
  cat("Assembling simulated data for compatability with denoising algorithm \n")
  sim_counts = c(sim_counts, unname(njs))
  uniqueSequences = c(unlist(isomiR_sequences), center_seqs)
  
  sim_rowdata = cbind(uniqueSequences, rep(miRNA, length(uniqueSequences))) %>% data.frame()
  colnames(sim_rowdata) = c("uniqueSequence", "miRNA_name")
  row.names(sim_rowdata) = NULL
  
  return(list(center_seqs=center_seqs, isomiR_sequences=isomiR_sequences, lambdas=lambdas, mus=mus, njs=njs, 
              sim_rowdata=sim_rowdata, sim_counts=sim_counts))
}

tst = generate_simulated_dataset(rowdata, countdata, miRNA, 2, transition_probs=transition_probs)
