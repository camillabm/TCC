##################
### LGG Cohort ###
##################

### Explanation
# Other glioma grades samples will be used for a comparative analysis
# There will be ramdomly selected 200 samples from each of them
# Normal tissue expression values will come from GTEx (there will be selected 200 samples from the previous 367 used)

###### Libraries
library(readxl)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)

###### ARGs
autophagy_genes <- read.csv("FAZER/ARGs/segunda parte/genes.csv") 
genes <- autophagy_genes$ensembl_gene_id

###### Cohorts
### TCGA GTEx
tcga <- read.csv("FAZER/Expression/results/matriX_gbmR_comb/TPM_HGNC.csv")
md_tcga <- read.csv("FAZER/Expression/tables/metadata/cohort/metadata_gbmR_comb.csv")

gene_col <- tcga$X
gtex_cols <- grep("^GTEX", colnames(tcga), value = TRUE)
gtex_sample <- sample(gtex_cols, 200)

tcga_gtex <- tcga[, colnames(tcga) %in% gtex_sample] # Normal tissues expression
tcga_gtex <- cbind(gene_col, tcga_gtex)
md_gtex_filtered <- md_tcga[md_tcga$barcode %in% colnames(tcga_gtex), ] # Normal tissues metadata

### LGG
lgg <- read.csv("FAZER/GBM/LGG/expression_lgg.csv")
md_lgg <- read.csv("FAZER/GBM/LGG/metadata_lggR.csv")

lgg_ARG  <- lgg[lgg$X %in% genes, ] # Filtering for selected genes

gene_values <- lgg_ARG$X
lgg_names <- setdiff(colnames(lgg_ARG), "X")
lgg_cols <- lgg_ARG[, colnames(lgg_ARG) %in% lgg_names]
colnames(lgg_cols) <- gsub("\\.", "-", colnames(lgg_cols))

G2 <- md_lgg %>%
  filter(paper_Grade == "G2") %>%
  slice_sample(n = min(200, nrow(.)))
G3 <- md_lgg %>%
  filter(paper_Grade == "G3") %>%
  slice_sample(n = min(200, nrow(.)))

barcodes_g2 <- G2$barcode
barcodes_g3 <- G3$barcode
barcodes_grade <- c(barcodes_g2, barcodes_g3)

lgg_sample <- lgg_cols[, colnames(lgg_cols) %in% barcodes_grade]
lgg_with_genes <- cbind(gene_values, lgg_sample)

### TPM normalization
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

counts_matrix <- lgg_with_genes
rownames(counts_matrix) <- counts_matrix$gene_values
counts_matrix$gene_values <- NULL

gene_lengths_filtered <- gene_lengths_df %>%
  filter(ENSEMBL %in% rownames(counts_matrix))
gene_ids_valid <- gene_lengths_df$ENSEMBL
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

count_data <- as.matrix(count_data)
fpkm_matrix <- sweep(count_data, 1, gene_lengths_kb, FUN = "/")     
fpkm_matrix <- sweep(fpkm_matrix, 2, colSums(count_data, na.rm = TRUE) / 1e6, FUN = "/")  

fpkm_matrix <- cbind(Gene_ID = gene_ids, fpkm_matrix)

fpkm_only <- as.matrix(fpkm_matrix[, -1])
fpkm_only <- fpkm_only[!apply(fpkm_only == "#####", 1, any), ]
fpkm_only <- apply(fpkm_only, 2, as.numeric)
tpm_matrix <- apply(fpkm_only, 2, function(x) x / sum(x) * 1e6)

tpm_matrix <- cbind(Gene_ID = gene_ids, as.data.frame(tpm_matrix))
tpm_matrix <- tpm_matrix[!duplicated(tpm_matrix$Gene_ID), ]

rownames(tpm_matrix) <- NULL
tpm_matrix <- tpm_matrix[!apply(is.na(tpm_matrix), 1, any), ]
tpm_matrix <- tpm_matrix %>%
  tibble::column_to_rownames(var = "Gene_ID") #CHANGE
tpm_matrix$Gene_ID <- NULL
tpm_matrix <- round(tpm_matrix)
### TPM normalization

### ENSEMBL to HGNC
tpm_matrix <- tpm_matrix %>%
  tibble::rownames_to_column("ensembl_gene_id")
df_merged <- tpm_matrix %>%
  left_join(autophagy_genes, by = "ensembl_gene_id")
df_merged <- df_merged[!duplicated(df_merged$HGNC), ]
rownames(df_merged) <- NULL
tpm_matrix <- df_merged %>%
  tibble::column_to_rownames("HGNC")

tpm_matrix$gene_values <- rownames(tpm_matrix)  
tpm_matrix <- tpm_matrix[, !(names(tpm_matrix) %in% c("ensembl_gene_id", "HGNC", "group","name", "function.", "influence.in.autophagy","reference", "Pathway"))]
### ENSEMBL to HGNC

lgg_with_genes <- tpm_matrix
gene_names <- rownames(lgg_with_genes)
lgg_with_genes$gene_values <- rownames(lgg_with_genes) 

lgg_normal <- tcga_gtex %>%
  filter(gene_col %in% gene_names)
new_genes <- lgg_normal$gene_col 
lgg_normal$gene_col <- NULL
md_lgg_normal <- md_tcga[md_tcga$barcode %in% colnames(lgg_normal), ]

lgg_actual <- lgg_with_genes %>%
  filter(gene_values %in% new_genes) %>%
  arrange(factor(gene_values, levels = new_genes))
lgg_actual$gene_values <- NULL
rownames(lgg_actual) <- NULL

lgg_final <- cbind(new_genes, lgg_normal, lgg_actual)

rownames(lgg_final) <- lgg_final$new_genes
lgg_final$new_genes <- NULL

md_lgg_filtered <- md_lgg[md_lgg$barcode %in% colnames(lgg_final), ]
md_lgg_filtered$project_id <- "TCGA-LGG"
md_lgg_filtered$condition <- "Tumor"
md_lgg_filtered <- dplyr::select(
  md_lgg_filtered,
  condition, project_id, tissue_or_organ_of_origin, barcode, paper_Grade
)

missing_cols <- setdiff(colnames(md_lgg_filtered), colnames(md_lgg_normal))

for (col in missing_cols) {
  md_lgg_normal[[col]] <- NA
}

md_lgg_normal <- md_lgg_normal[, colnames(md_lgg_filtered)]

md_final <- rbind(md_lgg_normal, md_lgg_filtered)

write.csv(lgg_final, "FAZER/GBM/LGG/lgg_exp.csv", row.names = TRUE)
write.csv(md_final, "FAZER/GBM/LGG/lgg_md.csv")

