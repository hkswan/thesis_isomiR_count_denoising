
load_mouse_miRNA_data = function(){
  load("/scratch/mmccall2_lab/miRNA_reference/Rdata/summarized_experiment_isomiR.RData")
  isomiR_se_object = se_object
  rm(se_object)
  countdata = isomiR_se_object@assays@data$counts  %>% data.frame()
  colnames(countdata) = str_remove_all(colnames(countdata), "X")
  rowdata = SummarizedExperiment::rowData(isomiR_se_object) %>% data.frame()
  colnames(rowdata) = c("uniqueSequence", "miRNA_name")
  return(list(isomiR_se_object=isomiR_se_object, rowdata=rowdata, countdata=countdata))
}

