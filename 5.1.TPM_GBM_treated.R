#########################
### Conversion to TPM ###
#########################

### Explanation
# Here, we eill manually convert expression values of GEO data when needed.
# Data will padronizated with TPM round values an HGNC symbol

##### Libraries
library(readxl)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)

##### Obtaining gene lengths
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths_all <- sum(width(reduce(exons_by_gene)))

gene_lengths_df <- data.frame(
  entrez_id = names(gene_lengths_all),
  gene_length = as.numeric(gene_lengths_all)
)

entrez_to_ensembl <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_lengths_df$entrez_id,
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENTREZID"
)

gene_lengths_df <- left_join(gene_lengths_df, entrez_to_ensembl, by = c("entrez_id" = "ENTREZID"))

gene_lengths_df <- gene_lengths_df %>% filter(!is.na(ENSEMBL)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

gene_lengths_df <- gene_lengths_df %>%
  mutate(gene_length_kb = gene_length / 1000)

### Expression matrix #Change
RAPA <- read_excel("FAZER/GBM/Analysis/before normalization/cgga_exp.xlsx")
write.csv(RAPA, "FAZER/GBM/Analysis/expression_u251.csv")
RAPA <- read.csv("FAZER/GBM/Analysis/before normalization/cgga_exp.csv")
RAPA$Gene <- RAPA$X
RAPA <- RAPA[!apply(is.na(RAPA), 1, any), ]
RAPA$Gene <- sub("\\..*", "", RAPA$Gene)
RAPA <- RAPA[!duplicated(RAPA$Gene), ]
rownames(RAPA) <- RAPA$Gene
RAPA$Gene <- NULL
RAPA$X <- NULL

### Integrating with expression matrix
counts_matrix <- RAPA
gene_lengths_filtered <- gene_lengths_df %>%
  filter(ENSEMBL %in% rownames(counts_matrix))#Change
gene_ids_valid <- gene_lengths_df$ENSEMBL #Change
counts_matrix <- counts_matrix[rownames(counts_matrix) %in% gene_ids_valid, ]

counts_matrix <- counts_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Gene_ID")
counts_with_lengths <- left_join(counts_matrix, gene_lengths_filtered, by = c("Gene_ID" = "ENSEMBL")) #Change
counts_with_lengths <- counts_with_lengths %>% filter(!is.na(gene_length_kb))
counts_with_lengths <- as.data.frame(counts_with_lengths)
counts_with_lengths <- counts_with_lengths[complete.cases(counts_with_lengths), ]

count_data <- counts_with_lengths[, !(names(counts_with_lengths) %in% c("Gene_ID", "entrez_id", "gene_length", "gene_length_kb", "SYMBOL"))] #Change

gene_lengths_kb <- counts_with_lengths$gene_length_kb
gene_lengths_kb <- as.numeric(gene_lengths_kb)
gene_ids <- counts_with_lengths$Gene_ID

### Calculating FPKM
count_data <- as.matrix(count_data)
fpkm_matrix <- sweep(count_data, 1, gene_lengths_kb, FUN = "/")     
fpkm_matrix <- sweep(fpkm_matrix, 2, colSums(count_data, na.rm = TRUE) / 1e6, FUN = "/")  

fpkm_matrix <- cbind(Gene_ID = gene_ids, fpkm_matrix)

########## for RPKM directly ###
RAPA$Gene_ID <- rownames(RAPA)
RAPA <- RAPA[, c("Gene_ID", setdiff(names(RAPA), "Gene_ID"))]
RAPA <- RAPA[complete.cases(RAPA), ]
gene_ids <- rownames(RAPA)
rpkm_matrix <- RAPA
########## for RPKM directly ###

### Calculating TPM
fpkm_only <- as.matrix(fpkm_matrix[, -1])
fpkm_only <- fpkm_only[!apply(fpkm_only == "#####", 1, any), ]
fpkm_only <- apply(fpkm_only, 2, as.numeric)
tpm_matrix <- apply(fpkm_only, 2, function(x) x / sum(x) * 1e6)

tpm_matrix <- cbind(Gene_ID = gene_ids, as.data.frame(tpm_matrix))
tpm_matrix <- tpm_matrix[!duplicated(tpm_matrix$Gene_ID), ]

########## for TPM directly ###
RAPA$Gene_ID <- rownames(RAPA)
tpm_matrix <- RAPA
########## for TPM directly ###

### Converting to HGNC (if necessary)
#Optional\
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", # Converting gene info to ensembl # ARRUMARRRRRRRRRR
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://grch37.ensembl.org")  
hgnc_mapping <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),  
  filters = "ensembl_gene_id",  
  values = tpm_matrix$Gene_ID,  
  mart = ensembl
)
tpm_matrix <- merge(tpm_matrix, hgnc_mapping, by.x = "Gene_ID", by.y = "ensembl_gene_id", all.x = TRUE)
tpm_matrix <- tpm_matrix[!duplicated(tpm_matrix$hgnc_symbol), ] #CHANGE
#Optional\
rownames(tpm_matrix) <- NULL
tpm_matrix <- tpm_matrix[!apply(is.na(tpm_matrix), 1, any), ]
tpm_matrix <- tpm_matrix %>%
  tibble::column_to_rownames(var = "hgnc_symbol") #CHANGE
tpm_matrix$Gene_ID <- NULL
tpm_matrix <- round(tpm_matrix)

write.csv(tpm_matrix, file = "FAZER/GBM/Analysis/TPM/expression_u251.csv", row.names = TRUE)
