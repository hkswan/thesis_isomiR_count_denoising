
ERCC_SE = readRDS("/scratch/mmccall2_lab/ERCC_UMI/panel_B_SE.rds")
keep_ind = which(rowData(ERCC_SE)$miRNA != "-")

ERCC_SE = ERCC_SE[keep_ind,]
ERCC_SE = DESeq2:: collapseReplicates(ERCC_SE, groupby = colData(ERCC_SE)$Run)

library(DESeq2)
uniq_mirnas = unique(rowData(ERCC_SE)$miRNA)
uniq_mirnas

miRNA_counts = matrix(0, nrow = length(uniq_mirnas), ncol = ncol(assay(ERCC_SE)))

miRNA_rowData = uniq_mirnas %>% data.frame()
colnames(miRNA_rowData) = "miRNA"

for(i in 1:length(uniq_mirnas)){
  ind = which(rowData(ERCC_SE)$miRNA == uniq_mirnas[i])
  x = colSums(assay(ERCC_SE)[ind,])
  miRNA_counts[i,] = x
}

colnames(miRNA_counts) = colData(ERCC_SE)$Run

ERCC_miRNA_SE = SummarizedExperiment(assays=list(counts=miRNA_counts), rowData = miRNA_rowData,
                                     colData = colData(ERCC_SE))
ERCC_miRNA_SE = collapseReplicates(ERCC_miRNA_SE, groupby = ERCC_miRNA_SE$GEO_Sample)

design_mat = colData(ERCC_miRNA_SE)$Pool

ERCC_dds = DESeqDataSet(ERCC_miRNA_SE, design = ~ Pool)
ERCC_dds <- DESeq(ERCC_dds)
results(ERCC_dds)

ratio_pool = read.csv(file = "/scratch/mmccall2_lab/ERCC.Rds/FINAL_Ratiometric_SynthA_and_SynthB-1.csv", sep="\t")
ratio_pool$A = sapply(ratio_pool$X.SAMPLE.A., function(x) return(as.numeric(strsplit(x, 'x')[[1]][1])))
ratio_pool$B = sapply(ratio_pool$X.SAMPLE.B., function(x) return(as.numeric(strsplit(x, 'x')[[1]][1])))
ratio_pool$ratio = ratio_pool$B/ratio_pool$A
ratio_pool$true_logFC = log(ratio_pool$ratio)
head(ratio_pool)


results_df = cbind(rowData(ERCC_dds)$miRNA, results(ERCC_dds)$log2FoldChange) %>% data.frame()
colnames(results_df) = c("miRNA", "log2FoldChange")
results_df$log2FoldChange = as.numeric(results_df$log2FoldChange)
results_df = mutate(results_df, est_logFC = log(2^log2FoldChange))

colnames(results_df)[1] = "ratio.seqID"

results_df = merge(ratio_pool, results_df, by = 'ratio.seqID') %>% head()

bias = results_df$est_logFC-results_df$true_logFC
sqerr = bias^2
mse = mean(sqerr)
mse

library(miRglmm)
ERCC_SE = SummarizedExperiment(assays = list(counts = ERCC_SE@assays@data$total_counts), 
                               rowData = rowData(ERCC_SE), colData = colData(ERCC_SE))
ERCC_SE = collapseReplicates(ERCC_SE, groupby = colData(ERCC_SE)$GEO_Sample)

#fit miRglmm with default filtering to ERCC data 
ERCC_mirglmm_fits = miRglmm(ERCC_SE, col_group = colData(ERCC_SE)$Pool, ncores = 48)


