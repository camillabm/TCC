#####################
### TCGA analysis ###
#####################

### Explanation
# We will use the previous cohort to the the TCGA analysis.
# STEPS:
  # (1) Conversion of raw counts to TPM
  # (2) Conversion of gene names from ENSEMBL to HGNC
  # (3) Analysis of ARGs 
  # (4) Analysis of gene sets with GSVA for profiles definitions
  # (5) Data plots as volcanos 

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
input_file <- "FAZER/Expression/tables/rna/cohort/matrix_coadR_comb.csv" #CHANGE
output_dir <- "FAZER/Expression/results/"
metadata_dir <- "FAZER/Expression/tables/metadata/cohort/"
arg_file <- "FAZER/ARGs/segunda parte/genes.csv"

###### ARGs
ARG <- read.csv(arg_file)
ARG$Pathway <- "Autophagy"
hgnc_genes <- unique(ARG$HGNC)

###### Colors for plottings
ordered_labels <- c("Initiation", "Nucleation", "Elongation", "Maturation", "Lysosomal transporters", 
                    "Lysosomal membrane", "Lysosomal enzimes", "CAMKK2", "TSC1", "TSC2", "AMPK", "PKA", "MTOR")
custom_colors <- c(
  "Initiation" = "#B39DDB",
  "Nucleation" = "#B39DDB",
  "Elongation" = "#B39DDB",
  "Maturation" = "#B39DDB",
  "Lysosomal transporters" = "#AED581",
  "Lysosomal membrane" = "#AED581",
  "Lysosomal enzimes" = "#AED581",
  "CAMKK2" = "#EF9A9A",
  "TSC1" = "#EF9A9A",
  "TSC2" = "#EF9A9A",
  "AMPK" = "#EF9A9A",
  "PKA" = "#90CAF9",
  "MTOR" = "#90CAF9"
)

###### Funtions to data conversions
  # Raw counts to TPM
  convert_to_fpkm <- function(data, gene_mapping_file) {
    data <- data.frame(Gene = rownames(data), data)
    rownames(data) <- NULL
    
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
      columns = c("ENSEMBL"),
      keytype = "ENTREZID"
    )
    
    gene_lengths_df <- left_join(gene_lengths_df, entrez_to_ensembl, by = c("entrez_id" = "ENTREZID"))
    gene_lengths_df <- gene_lengths_df %>% filter(!is.na(ENSEMBL)) %>%
      distinct(ENSEMBL, .keep_all = TRUE)
    
    gene_lengths_df <- gene_lengths_df %>%
      mutate(gene_length_kb = gene_length / 1000)
    
    data_hgnc <- merge(data, gene_lengths_df, by.x = "Gene", by.y = "ENSEMBL", all.x = TRUE)
    data_hgnc <- data_hgnc[!is.na(data_hgnc$gene_length_kb), ]
    data_hgnc$entrez_id <- NULL
    data_hgnc$gene_length <- NULL
    
    expr_matrix <- data_hgnc[ , -1]
    expr_matrix$gene_length_kb <- NULL
    rownames(expr_matrix) <- data_hgnc$Gene
    expr_matrix <- as.matrix(expr_matrix)
  
    gene_lengths_kb <- data_hgnc$gene_length_kb
    fpkm_matrix <- sweep(expr_matrix, 1, gene_lengths_kb, FUN = "/")     
    fpkm_matrix <- sweep(fpkm_matrix, 2, colSums(expr_matrix, na.rm = TRUE) / 1e6, FUN = "/")
    fpkm_matrix <- na.omit(fpkm_matrix)
    
    tpm_matrix <- apply(fpkm_matrix, 2, function(x) x / sum(x) * 1e6)
    tpm_matrix <- round(tpm_matrix)
    message("Genes na matriz de expressão: ", nrow(data))
    message("Genes com comprimento conhecido: ", nrow(gene_lengths_df))
    
    return(list(fpkm = fpkm_matrix, tpm = tpm_matrix))
  }

  # ENSEMBL to HGNC
  convert_to_hgnc <- function(tpm_matrix) {
    message("Iniciando conversão de ENSEMBL para HGNC...")
    ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                       dataset = "hsapiens_gene_ensembl",
                       host = "https://grch37.ensembl.org") 
    message("Conexão com BioMart estabelecida.")
    
    tpm_matrix <- as.data.frame(tpm_matrix)
    tpm_matrix$ensembl_gene_id <- rownames(tpm_matrix)
    rownames(tpm_matrix) <- NULL
    
    hgnc_mapping <- getBM(
      attributes = c("hgnc_symbol", "ensembl_gene_id"),
      filters = "ensembl_gene_id",
      values = tpm_matrix$ensembl_gene_id,
      mart = ensembl
    )
    message("Mapeamento obtido: ", nrow(hgnc_mapping), " genes mapeados.")
  
    tpm_matrix_hgnc <- merge(tpm_matrix, hgnc_mapping, by = "ensembl_gene_id", all.x = TRUE)
    message("Mapeamento realizado: ", nrow(tpm_matrix_hgnc), " genes merged.")
    
    tpm_matrix_hgnc <- tpm_matrix_hgnc[!is.na(tpm_matrix_hgnc$hgnc_symbol) & tpm_matrix_hgnc$hgnc_symbol != "", ]
    tpm_matrix_hgnc <- tpm_matrix_hgnc[!duplicated(tpm_matrix_hgnc$hgnc_symbol), ]
    
    rownames(tpm_matrix_hgnc) <- tpm_matrix_hgnc$hgnc_symbol
    tpm_matrix_hgnc$ensembl_gene_id <- NULL
    tpm_matrix_hgnc$hgnc_symbol <- NULL
    
    message("Conversão finalizada. Genes com símbolo HGNC: ", nrow(tpm_matrix_hgnc))
    return(tpm_matrix_hgnc)
  }
  
###### Data processing and analysis
  sample_name <- tools::file_path_sans_ext(basename(input_file))
  message("Processando: ", sample_name)
    
  sample_outdir <- file.path(output_dir, sample_name)
  dir.create(sample_outdir, recursive = TRUE, showWarnings = FALSE)
  
  csv_file <- input_file
  
  data <- read.csv(csv_file)
  rownames(data) <- data$X
  data$X <- NULL
  colnames(data) <- ifelse(    #Sample ID padronization
    grepl("^TCGA", colnames(data)),        
    gsub("\\.", "-", colnames(data)),      
    colnames(data)                         
  ) 
  
  results <- convert_to_fpkm(data, arg_file)
  fpkm_matrix <- results$fpkm
  tpm_matrix <- results$tpm
  
  write.csv(fpkm_matrix, file.path(sample_outdir, "FPKM_results.csv"))
  write.csv(tpm_matrix, file.path(sample_outdir, "TPM_results.csv"))
  
  tpm_matrix_hgnc <- convert_to_hgnc(tpm_matrix)
  colnames(tpm_matrix_hgnc) <- ifelse(    #Sample ID padronization
    grepl("^TCGA", colnames(tpm_matrix_hgnc)),        
    gsub("\\.", "-", colnames(tpm_matrix_hgnc)),      
    colnames(tpm_matrix_hgnc)                         
  ) 
  
  write.csv(tpm_matrix_hgnc, file.path(sample_outdir, "TPM_HGNC.csv"))
  
  cohort_code_expr <- substr(sample_name, 8, 11)
  
  # Metadata defining
  metadata_files <- list.files(metadata_dir, pattern = "\\.csv$", full.names = TRUE)
  matched_meta <- metadata_files[substr(basename(metadata_files), 10, 13) == cohort_code_expr]
  
  if (length(matched_meta) != 1) {
    warning("Metadado não encontrado ou múltiplos encontrados para: ", sample_name)
    next
  }
  
  md <- read.csv(matched_meta)
  md_filtered <- md %>% filter(barcode %in% colnames(tpm_matrix_hgnc))
  md_filtered <- md_filtered[match(colnames(tpm_matrix_hgnc), md_filtered$barcode), ]
  group <- md_filtered$condition
  
  # Design for DEGs analysis
  design <- model.matrix(~ group)
  coef_name <- colnames(design)[2]
  
  # DEG analysis
  log_tpm_n1 <- log2(tpm_matrix_hgnc + 1)
  fit_n1 <- lmFit(log_tpm_n1, design)
  fit_n1 <- eBayes(fit_n1)
  res_n1 <- topTable(fit_n1, coef = coef_name, number = Inf, adjust = "BH")
  res_n1$gene <- rownames(res_n1)
  write.csv(res_n1, file.path(sample_outdir, "results_n1.csv"), row.names = FALSE)
  
  control_samples <- md_filtered$barcode[group == "Normal"]
  control_means <- rowMeans(log_tpm_n1[, control_samples, drop = FALSE])
  log2fc_indiv <- sweep(log_tpm_n1, 1, control_means, "-")
  write.csv(log2fc_indiv, file.path(sample_outdir, "n1_individual.csv"))
  
  # GSVA
  RAPA <- read.csv("FAZER/Expression/results/matrix_coadR_comb/TPM_HGNC.csv") #CHANGE
  names(RAPA)[names(RAPA) == "X"] <- "Gene"
  RAPA <- merge(data.frame(Gene = unique(hgnc_genes)), RAPA, by = "Gene", all.x = TRUE)
  rownames(RAPA) <- RAPA$Gene
  RAPA$Gene <- NULL
  RAPA[is.na(RAPA)] <- 0
  colnames(RAPA) <- ifelse(    #Sample ID padronization
    grepl("^TCGA", colnames(RAPA)),        
    gsub("\\.", "-", colnames(RAPA)),      
    colnames(RAPA)                         
  ) 
  
  RAPA_log2 <- log2(RAPA + 1) 
  data <- RAPA_log2
  
  data <- as.matrix(data)
  gene_var <- apply(data, 1, var)
  data_filtered <- data[gene_var > 0, ]
  
  md_filtered <- md_filtered[match(colnames(data_filtered), md_filtered$barcode), ]
  
  group_gsva <- ifelse(md_filtered$condition == "Tumor", 1, 0)
  design <- cbind(Intercept = 1, Tumor_vs_Control = group)
  
  initiation <- ARG %>% filter(group == "Initiation") %>% pull(HGNC)
  nucleation <- ARG %>% filter(group == "Nucleation") %>% pull(HGNC)
  elongation <- ARG %>% filter(group == "Elongation") %>% pull(HGNC)
  maturation <- ARG %>% filter(group == "Maturation") %>% pull(HGNC)
  enzimes <- ARG %>% filter(group == "Lysosomal enzymes") %>% pull(HGNC)
  membrane <- ARG %>% filter(group == "Lysosomal membrane") %>% pull(HGNC)
  transporters <- ARG %>% filter(group == "Lysosomal transporters") %>% pull(HGNC)
  reg_pos <- ARG %>% filter(influence.in.autophagy == "positive") %>% pull(HGNC)
  reg_neg <- ARG %>% filter(influence.in.autophagy == "negative") %>% pull(HGNC)
  
  autophagy <- c(initiation, nucleation, elongation, maturation)
  lisosomal_activity <- c(enzimes, membrane, transporters)
  
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
  gsva_es <- as.data.frame(gsva_es)
  write.csv(gsva_es, file = file.path(sample_outdir, "GSVA_results.csv"))
  
     # Individual
  control_samples <- md_filtered$barcode[md_filtered$condition == "Normal"]
  log_tpm_control <- gsva_es[, colnames(gsva_es) %in% control_samples, drop = FALSE]
  control_means <- rowMeans(log_tpm_control, na.rm = TRUE)
  log2fc_indiv_gsva <- sweep(gsva_es, 1, control_means, "-")
  tcga_columns <- grep("^GTEX", colnames(log2fc_indiv_gsva), value = TRUE)
  log2fc_indiv_tcga <- log2fc_indiv_gsva[, tcga_columns]
  write.csv(log2fc_indiv_tcga, file.path(sample_outdir, "GSVA_results_normal.csv"))
  
    # Grouped
  design <- model.matrix(~ group_gsva)
  fit <- lmFit(gsva_es, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "group_gsva", number = Inf)
  results$label <- rownames(results)
  write.csv(results, file.path(sample_outdir, "GSVA_grouped.csv"))

###### Data plots
  # Volcano
  res_n1$source <- "ARG"
  res_n1 <- res_n1 %>%
    as.data.frame() %>%
    tibble::rownames_to_column("HGNC")
  res_n1 <- res_n1 %>%
    left_join(ARG[, c("HGNC", "group","influence.in.autophagy")], by = "HGNC")
  res_n1 <- res_n1 %>% filter(!is.na(group))
  res_n1 <- res_n1 %>% distinct(HGNC, .keep_all = TRUE)
  res_n1 <- res_n1 %>%
    mutate(group = ifelse(group == "Regulation" & influence.in.autophagy == "positive", "Reg_pos",
                          ifelse(group == "Regulation" & influence.in.autophagy == "negative", "Reg_neg", group)))
  rownames(res_n1) <- res_n1$HGNC
  res_n1$HGNC <- NULL
  res_n1$influence.in.autophagy <- NULL
  res_n1$group <- gsub(" ", "_", res_n1$group)
  res_n1$gene <- NULL
  res_n1$gene <- rownames(res_n1)
  
  results$source <- "GSVA"
  results$group <- NULL
  rownames(results) <- rownames(results) %>%
    sub("^L\\.", "Lysosomal_", .) %>%
    sub("Enzymes", "Lysosomal_enzymes", .) %>%
    sub("Membrane", "Lysosomal_membrane", .) %>%
    sub("Transporters", "Lysosomal_transporters", .)
  results$group <- rownames(results)
  results$gene <- NULL
  results$gene <- rownames(results)
  
  res_combined <- res_n1[, c("gene", "logFC", "P.Value", "source", "group")]
  results_combined <- results[, c("gene", "logFC", "P.Value", "source", "group")]
  
  combined_data <- rbind(res_combined, results_combined)
  
  group_colors <- c(
    Autophagy = "#1B9E77",           
    Initiation = "#D95F02",          
    Nucleation = "#7570B3",          
    Maturation = "#E7298A",         
    Elongation = "#66A61E",          
    
    Lisosomal_activity = "#E6AB02",  
    Lysosomal_enzymes = "#A6761D",   
    Lysosomal_membrane = "#666666",  
    Lysosomal_transporters = "#1F78B4", 
    
    Reg_pos = "#B22222",             
    Reg_neg = "#00008B"              
  )
  
  combined_data$color <- group_colors[combined_data$group]
  combined_data$lab <- combined_data$gene
  
  point_colors <- setNames(
    group_colors[combined_data$group], 
    rownames(combined_data)
  )
  
  hline <- ifelse(all(is.na(combined_data$P.Value)) || all(combined_data$P.Value == 1), 
                  NA, -log10(0.05))

  volcano_plot <- EnhancedVolcano(combined_data,
                                  lab = rep("", nrow(combined_data)),
                                  lab = ifelse(combined_data$source == "GSVA", combined_data$gene, ""),
                                  x = 'logFC',
                                  y = 'P.Value',
                                  pCutoff = 0.05,
                                  FCcutoff = 0,
                                  pointSize = 3,
                                  labSize = 5,
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.5, 
                                  min.segment.length = 0.001,
                                  colAlpha = 0.8,
                                  colCustom = point_colors,
                                  vlineWidth = 1,
                                  hline = 0,
                                  hlineCol = 'black',
                                  hlineType = 'solid',
                                  hlineWidth = 1.2,
                                  legendPosition = 'none',
                                  max.overlaps = Inf
  ) + 
    labs(
      x = "logFC", 
      y = "-log10(p-value)"
    ) +
    theme_minimal() + 
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA,linewidth = 1.2),
      axis.line = element_line(size = 1.2),
      axis.ticks = element_line(size = 1.2),
      axis.text = element_text(size = 14, color = "black"),        
      axis.title = element_text(size = 16, face = "plain", color = "black"),
      legend.position = "none"
    )
  
  volcano_grob <- ggplotGrob(volcano_plot)
  
  labels <- unique(combined_data$group[combined_data$source == "GSVA"])
  colors <- unique(point_colors[combined_data$source == "GSVA"])
  
  legend_grob <- legendGrob(
    labels = labels,
    pch = 16,
    gp = gpar(col = colors, fontsize = 12),
    nrow = length(labels)
  )
  
  png(file.path(sample_outdir, "volcano_correct.png"), width = 3000, height = 2000, res = 300)
  

  grid.arrange(volcano_grob, legend_grob,
               ncol = 2,
               widths = c(4, 1))  
  
  dev.off()
  
  # # UMAP
  # MD <- read.csv("FAZER/Expression/results/survival/gbm_states.csv")
  # MD$X <- NULL
  # 
  # data <- gsva_es
  # 
  # data <- ComBat(dat = data, batch = MD$project_id)
  # 
  # set.seed(125)
  # 
  # umap_input <- t(data)
  # 
  # umap_result <- umap(
  #   umap_input,
  #   n_neighbors = 15,  # Ajusta o número de vizinhos (valores menores = mais locais, valores maiores = mais globais)
  #   min_dist = 0.001,  # Controla a dispersão dos pontos (menor = clusters mais compactos, maior = mais espalhado)
  #   spread =  2,  # Controla a separação dos grupos
  #   metric = "correlation",  # Pode testar "euclidean", "correlation" ou "cosine"
  #   n_components = 2,  # Número de dimensões do UMAP
  #   learning_rate = 0.1,
  #   fast_sgd = TRUE,  # Torna o cálculo mais rápido
  #   verbose = TRUE  # Mostra informações durante o processo
  # )
  # 
  # umap_df <- data.frame(
  #   UMAP1 = umap_result[, 1],
  #   UMAP2 = umap_result[, 2],
  #   Condition = MD$condition,
  #   Project = MD$project_id,
  #   State = MD$state,  
  #   Treatment = MD$treatments_pharmaceutical_treatment_or_therapy
  # )
  # 
  # 
  # umap_df$ColorState <- "gray"
  # 
  # umap_df$ColorState[umap_df$State == 2] <- "red"
  # 
  # umap_df$ColorState[umap_df$State != 2 & umap_df$Treatment == "yes"] <- "blue"
  # 
  # umap_df$ColorState[umap_df$State == 2 & umap_df$Treatment == "yes"] <- "brown"
  # 
  # umap_df$ColorState <- factor(umap_df$ColorState, levels = c("gray", "red", "brown"))
  # 
  # colors <- c(
  #   "gray" = "gray",
  #   "red" = "red",
  #   "brown" = "brown",
  #   "blue" = "blue"
  # )
  # 
  # umap_df$State <- as.factor(umap_df$State)
  # 
  # UMAP <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = ColorState, shape = Condition)) +
  #   geom_point(size = 2, alpha = 0.7) +
  #   scale_color_manual(values = colors) +
  #   scale_shape_manual(values = c("Tumor" = 16, "Normal" = 15)) +  # shapes personalizados
  #   theme_minimal(base_family = "Arial") +
  #   theme(
  #     plot.title = element_text(size = 16, face = "bold"),
  #     legend.title = element_text(size = 12),
  #     legend.text = element_text(size = 10),
  #     panel.grid = element_blank()
  #   ) +
  #   labs(
  #     color = "State & Treatment",
  #     shape = "Condition"
  #   )
  # 
  # ggsave(
  #   filename = file.path(sample_outdir, "UMAP_log2fc_indiv.png"),
  #   plot = UMAP,
  #   width = 8, height = 6, dpi = 300, bg = "white"
  # )
  