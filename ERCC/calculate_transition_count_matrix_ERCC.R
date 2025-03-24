source("/scratch/hswan/thesis_isomiR_count_denoising/Code/denoise_isomiR_counts_update_partitions_at_end.R")

args = commandArgs(trailingOnly = TRUE)
miRNA_idx = as.numeric(args[1])

#load data 
ERCC_SE = readRDS("/scratch/mmccall2_lab/ERCC_UMI/panel_B_SE.rds")
keep_ind = which(rowData(ERCC_SE)$miRNA != "-")

ERCC_SE = ERCC_SE[keep_ind,]

#update column names to what functions are expecting: 
rowdat = rowData(ERCC_SE)
colnames(rowdat) = c('uniqueSequence', 'miRNA_name')

ERCC_SE = SummarizedExperiment(assays = list(counts = ERCC_SE@assays@data$total_counts),
                               rowData = rowdat, colData = colData(ERCC_SE))

unique_mirnas = unique(rowData(ERCC_SE)$miRNA_name)

miRNA = unique_mirnas[miRNA_idx]
cat("miRNA:", miRNA, "\n")

counts = rowSums(assay(ERCC_SE))
count_df = data.frame(count=counts, rowData(ERCC_SE))


isomiR_subset = filter(count_df, miRNA_name == miRNA)

center_seq = filter(isomiR_subset, count == max(count)) %>% select(., uniqueSequence) %>% unlist() %>% unname()
isomiR_seqs = filter(isomiR_subset, uniqueSequence != center_seq) %>% select(., uniqueSequence) %>% unlist() %>% unname()

alignments = lapply(isomiR_seqs, Biostrings::pairwiseAlignment, pattern = center_seq)

transitions = lapply(alignments, get_transitions)

transition_counts = initialize_transition_counts_matrix(5)

if(length(transitions) != 0){
  for(t in transitions){
    alignment_length = ncol(t)
    for(i in 1:alignment_length){
      trns = get_transition_at_idx(i, t)
      x=trns$x
      y=trns$y
      if(i == 1 & x != "-" & y != "-"){
        transition_counts[x,y]=transition_counts[x,y]+1
        transition_counts['-', '-'] = transition_counts['-', '-']+1
      }
      if(i == alignment_length & x != "-" & y !=  "-"){
        transition_counts['-', '-'] = transition_counts['-', '-'] + 1 
        transition_counts[x,y] = transition_counts[x,y]+1
      } else{
        transition_counts[x,y]=transition_counts[x,y]+1
      }
      
    }
  }
}

dir_name = "/scratch/hswan/thesis_isomiR_count_denoising/ERCC/transition_count_matrices/"
fname = paste0(miRNA, "_transition_counts.Rds")
cat("Saving transition count matrix to ", paste0(dir_name, fname), "\n")
saveRDS(transition_counts, paste0(dir_name, fname))
