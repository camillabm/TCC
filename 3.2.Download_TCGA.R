##########################
### Download TCGA DATA ###
##########################

### Explanation
# Data download directly from GDC Data Portal
# Both expression and clinical data will be downloaded and padronized

####### Libraries
library(SummarizedExperiment)
library(dplyr)
library (biomaRt)
library(TCGAbiolinks)

####### Project identification
project <- "TCGA-LGG" # Changed here to other projects (GBM, BRCA, LUAD, COAD, LGG)

###### Download
# RNA-seq data (transcriptomics)
query_rna <- GDCquery(
  project = project,
  experimental.strategy = "RNA-Seq",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query = query_rna, method = "api", files.per.chunk = 1) 
data_rna <- GDCprepare(query_rna)

##
metadata <- as.data.frame(colData(data_rna)) # extract metadata
metadata$sample_type <- substr(metadata$sample, 14, 15) # extract characters 14 and 15
metadata <- metadata[metadata$sample_type %in% c("01", "11"), ] # filter samples 01 and 11
metadata$condition <- ifelse(metadata$sample_type == "01", "Tumor", "Normal") # define condition
metadata$sample_type <- NULL 

metadata[] <- lapply(metadata, function(x) { 
  if (is.list(x)) {
    return(sapply(x, function(y) paste(y, collapse = ", ")))
  }
  return(x)
})
write.csv(metadata, "metadata_lggR.csv", row.names = TRUE) 

##
matrix <- assay(data_rna) # extract expression matrix

valid_samples <- metadata$barcode # identify samples in metadata
expression_matrix_filtered <- matrix[, colnames(matrix) %in% valid_samples] # filter for samples in metadata

rownames(expression_matrix_filtered) <- sub("\\..*", "", rownames(expression_matrix_filtered)) 
rownames(expression_matrix_filtered) <- as.character(rownames(expression_matrix_filtered)) 
any(duplicated(rownames(expression_matrix_filtered)))  
expression_matrix_filtered <- expression_matrix_filtered[!duplicated(rownames(expression_matrix_filtered)), ] 

write.csv(expression_matrix_filtered, "expression_lgg.csv", row.names = TRUE)

# Clinical Data
options(timeout = 300)
clinical_data <- GDCquery_clinic(
  project = project, 
  type = "clinical")

clinical_data_df <- as.data.frame(clinical_data)
clinical_data_df <- as.data.frame(lapply(clinical_data, function(x) {
  if (is.list(x)) sapply(x, toString) else x
}))
write.csv(clinical_data_df, "FAZER/Clinical data/clinical_data_lgg.csv")


