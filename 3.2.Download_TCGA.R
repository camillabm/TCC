##########################
### Download TCGA DATA ###
##########################

####### Libraries
library(SummarizedExperiment)
library(dplyr)
library (biomaRt)
library(TCGAbiolinks)

####### Project identification
project <- "TCGA-LGG" # Changed here to other projects (GBM, BRCA, LUAD, COAD)

###### Download
# Mutation Data (genomics)
query_mutation <- GDCquery(
  project = project,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )

GDCdownload(query_mutation, method = "api", files.per.chunk = 1)
mutation_data <- GDCprepare(query_mutation)
write.csv(mutation_data, "Mutation/tables/mutation_coad.csv", row.names = TRUE)
  
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

# Methylation Data (epigenomics)
BiocManager::install(c("sesameData", "sesame", "GenomeInfoDb", "GenomicRanges"))

library (sesameData)
library (sesame)
sesameDataCache() #rodar apenas na primeira vez que instala o pacote
library (GenomeInfoDb)
library (GenomicRanges)

query_met <- GDCquery(
  project = project,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

GDCdownload(query = query_met, method = "api", files.per.chunk = 1) #deixou mais leve pra baixar
data_met <- GDCprepare(query_met)

##
metadata <- as.data.frame(colData(data_met)) # extrai metadados
metadata$sample_type <- substr(metadata$sample, 14, 15) #extrai os caracteres 14 e 15
metadata <- metadata[metadata$sample_type %in% c("01", "11"), ] #filtra apenas amostras 01 e 11
metadata$condition <- ifelse(metadata$sample_type == "01", "Tumor", "Normal") #atribui condição
metadata$sample_type <- NULL #remove a coluna sample_type

metadata[] <- lapply(metadata, function(x) { # Converte listas para caracteres (ou outro tipo, conforme necessário)
  if (is.list(x)) {
    return(sapply(x, function(y) paste(y, collapse = ", ")))
  }
  return(x)
})
write.csv(metadata, "Methylation/metadata/metadata_coadM.csv", row.names = TRUE) #salva metadados em csv

##
beta_values <- assay(data_met) #extrai os beta values
beta_values <- beta_values[, metadata$barcode] #filtra só as amostras que tem nos metadados
rownames(beta_values) <- sub("\\..*", "", rownames(beta_values)) #normaliza sufixos dos genes

write.csv(beta_values, "Methylation/met/methylation_coad.csv")

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


