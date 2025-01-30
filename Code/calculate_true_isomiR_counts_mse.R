suppressMessages(source("/scratch/hswan/thesis_isomiR_count_denoising/Code/denoise_isomiR_counts_WORKING_FUNCTION.R"))

partition_path = "/scratch/hswan/thesis_isomiR_count_denoising/sims/01_24_2025/multi_partitions/partition_objs"
partition_files = list.files(partition_path)

successful_miRNAs = sapply(partition_files, function(x) return(strsplit(x, "_partition")[[1]][[1]]))

tst_miRNA = successful_miRNAs[1]

data_path = "/scratch/hswan/thesis_isomiR_count_denoising/data/simulated_data/length_variants_multi_partition/01_24_2025"
data_files = list.files(data_path)

datafilename = data_files[grepl(tst_miRNA, data_files)]
partitionsfilename = partition_files[grepl(tst_miRNA, partition_files)]

partition_objs = paste0(partition_path, "/", partitionsfilename, collapse="") %>% readRDS
data_list = paste0(data_path, "/", datafilename, collapse="") %>% readRDS()

mses = vector()
for(i in 1:length(partition_objs)){
  df = partition_objs[[i]]$partition_df
  #create true_partition column
  df$true_partition = rep(0, nrow(df))
  #get associated list of isomiR sequences:
  isomiR_seqs = data_list[[i]]$isomiR_seqs
  for(i in 1:length(isomiR_seqs)){
    df$true_partition[df$uniqueSequence == names(isomiR_seqs)[i]] = i
    df$true_partition[df$uniqueSequence %in% isomiR_seqs[[i]]] = i
  }
  true_counts = sapply(1:3, function(x) filter(df, true_partition ==x) %>% select(., count) %>% unlist() %>% sum())
  infer_counts = sapply(1:3, function(x) filter(df, partition ==x) %>% select(., count) %>% unlist() %>% sum())
  mses = c(mses, mean((true_counts-infer_counts)^2))
}


files = list.files("/scratch/hswan/thesis_isomiR_count_denoising/sims/01_24_2025/multi_partitions/upd8_count_estimation/Mmu-Let-7-P1b_3p*")

upd8d_counts=list()
for(i in 1:length(files)){
  upd8d_counts[[i]] = paste0("/gpfs/fs2/scratch/hswan/thesis_isomiR_count_denoising/sims/01_24_2025/multi_partitions/upd8_count_estimation/Mmu-Let-7-P1b_3p*/",
         files[i], collapse="") %>% readRDS()
}

upd8d_mses = vector(length=length(upd8d_counts))
for(i in 1:length(upd8d_counts)){
  upd8d_mses[i] = mean((upd8d_counts[[i]]['true_counts',]-upd8d_counts[[i]]['upd8_counts',])^2)
}

df = data.frame(x=rep(1:length(upd8d_counts),2), y=c(mses, upd8d_mses), z=c(rep("old", length(upd8d_counts)), rep("new", 
                                                                            length(upd8d_counts))))
df$z = as.factor(df$z)

ggplot(df, aes(x=x, y=log(y), col=z)) + geom_point() + geom_abline(slope = 0, intercept = mean(log(upd8d_mses)), col='red') + 
  geom_abline(slope = 0, intercept = mean(log(mses)), col='lightblue3')

sum(upd8d_mses <= mses)
wilcox.test(upd8d_mses-mses)
