#########################
### ARG GBM treatment ###
#########################

### Explanation
# We will use the previous TPM converted GEO datasets with chemotherapy glioma treated samples
# STEPS:
# (1) Filtering for ARGs
# (2) Analysis of ARGs 
# (3) Analysis of gene sets with GSVA for profiles definitions

###### Libraries
library(limma)
library(EnhancedVolcano)
library(umap)
library(uwot)
library(GSVA)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(tidyr)
library(biomaRt)
library(sva)
library(readxl)
library(grid)
library(gridExtra)

##### Paths and data
input_file <- "FAZER/GBM/Analysis/TPM/expression_u251.csv" #CHANGE ########
output_dir <- "FAZER/GBM/Results/U251"
arg_file <- "FAZER/ARGs/segunda parte/genes.csv"

###### Metadata
metadata_file <- read_excel("FAZER/GBM/Data/GEO/metadata.xlsx")
write.csv(metadata_file, "FAZER/GBM/Data/GEO/metadata.csv")
md <- read.csv("FAZER/GBM/Data/GEO/metadata.csv")
md_fil <- md[md$cell == "U251", ]
#md_fil <- md_fil[md_fil$treatment == "Radiation", ]
md_fil$X <- NULL

###### ARGs
ARG <- read.csv(arg_file)
ARG$Pathway <- "Autophagy"
hgnc_genes <- unique(ARG$HGNC)

###### Analysis
csv_file <- read.csv(input_file)
csv_file_fil <- csv_file[csv_file$X %in% hgnc_genes, ]
rownames(csv_file_fil) <- csv_file_fil$X
csv_file_fil$X <- NULL
csv_file_fil <- as.data.frame(csv_file_fil)

n_colunas <- ncol(csv_file_fil)
tamanho_bloco <- 4

inicios <- seq(1, n_colunas, by = tamanho_bloco)

for (i in seq_along(inicios)) {
  ini <- inicios[i]
  fim <- min(ini + tamanho_bloco - 1, n_colunas)
  bloco <- csv_file_fil[, ini:fim]
  
  assign(paste0("df_", i), bloco, envir = .GlobalEnv)
}

n_blocos <- length(seq(1, ncol(csv_file_fil), by = 6))

for (i in seq_len(n_blocos)) {
  df_i <- get(paste0("df_", i))
  ids_para_filtrar <- colnames(df_i)
  md_fil_i <- md_fil[md_fil$ID %in% ids_para_filtrar, , drop = FALSE]
  assign(paste0("md_fil_", i), md_fil_i, envir = .GlobalEnv)
  assign(paste0(i, "tc"), t(md_fil_i), envir = .GlobalEnv)
}

run_analysis_block <- function(df_expr, md_fil, prefix, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  log_expr <- log2(df_expr + 1)
  group <- factor(md_fil$condition)
  design <- model.matrix(~ group)
  coef_name <- colnames(design)[2]
  
  # DEG analysis
  fit <- lmFit(log_expr, design)
  fit <- eBayes(fit)
  res_n1 <- topTable(fit, coef = coef_name, number = Inf, adjust = "BH")
  res_n1$gene <- rownames(res_n1)
  write.csv(res_n1, file.path(output_dir, paste0(prefix, "_results_n1.csv")), row.names = FALSE)
  
  # log2FC individual da expressão
  control_samples <- md_fil$ID[md_fil$condition == "Control"]
  control_means <- rowMeans(log_expr[, control_samples, drop = FALSE])
  log2fc_indiv <- sweep(log_expr, 1, control_means, "-")
  log2fc_indiv$gene <- rownames(log2fc_indiv)
  write.csv(log2fc_indiv, file.path(output_dir, paste0(prefix, "_n1_individual.csv")), row.names = FALSE)
  
  # Preparar dados para GSVA
  RAPA <- data.frame(Gene = rownames(df_expr), df_expr)
  rownames(RAPA) <- RAPA$Gene
  RAPA$Gene <- NULL
  RAPA[is.na(RAPA)] <- 0
  
  RAPA_log2 <- log2(RAPA + 1)
  data <- as.matrix(RAPA_log2)
  data_filtered <- data[apply(data, 1, var) > 0, , drop = FALSE]
  
  md_filtered <- md_fil[match(colnames(data_filtered), md_fil$ID), , drop = FALSE]
  group_gsva <- ifelse(md_filtered$condition == "Treated", 1, 0)
  design_gsva <- model.matrix(~ group_gsva)
  
  # Usar sua lista de genes aqui (substitua ARG e listas correspondentes)
  autophagy_genes <- list(
    Autophagy = autophagy,
    Initiation = initiation,
    Nucleation = nucleation,
    Elongation = elongation,
    Maturation = maturation,
    Lisosomal_activity = lisosomal_activity,
    Enzymes = enzimes,
    Membrane = membrane,
    Transporters = transporters,
    Reg_pos = reg_pos,
    Reg_neg = reg_neg
  )
  
  gsvapar <- gsvaParam(data_filtered, autophagy_genes)
  gsva_es <- suppressMessages(suppressWarnings(gsva(gsvapar)))
  gsva_es_df <- as.data.frame(gsva_es)
  gsva_es_df$GeneSet <- rownames(gsva_es_df)
  write.csv(gsva_es_df, file.path(output_dir, paste0(prefix, "_GSVA_results.csv")), row.names = FALSE)
  
  # log2FC individual - GSVA para amostras normais
  control_samples_gsva <- md_filtered$ID[md_filtered$condition == "Control"]
  log_tpm_control <- gsva_es[, colnames(gsva_es) %in% control_samples_gsva, drop = FALSE]
  control_means <- rowMeans(log_tpm_control, na.rm = TRUE)
  log2fc_indiv_gsva <- sweep(gsva_es, 1, control_means, "-")
  log2fc_indiv_gsva_df <- as.data.frame(log2fc_indiv_gsva)
  log2fc_indiv_gsva_df$GeneSet <- rownames(log2fc_indiv_gsva_df)
  write.csv(log2fc_indiv_gsva_df, file.path(output_dir, paste0(prefix, "_GSVA_results_normal.csv")), row.names = FALSE)
  
  # Grouped GSVA
  fit_gsva <- lmFit(gsva_es, design_gsva)
  fit_gsva <- eBayes(fit_gsva)
  res_gsva <- topTable(fit_gsva, coef = "group_gsva", number = Inf)
  res_gsva$label <- rownames(res_gsva)
  write.csv(res_gsva, file.path(output_dir, paste0(prefix, "_GSVA_grouped.csv")), row.names = FALSE)
  
  message("Análise concluída e arquivos salvos para ", prefix)
}

run_analysis_block(df_1, md_fil_1, prefix = "df1", output_dir)
