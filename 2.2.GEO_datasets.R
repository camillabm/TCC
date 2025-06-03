####################
### GEO Datasets ###
####################

##### Libraries
library(readxl)
library(limma)
library(ggplot2)

##### Data download
# The GEO Data were directly downloaded from the webserver and the expession tables were adapted as needed.

# The metadata table contains the necessary metadata for all projects used.
#md <- read_excel("FAZER/Autophagy states datasets/metadata/metadata_v2.xlsx")
#write.csv(md, "FAZER/Autophagy states datasets/metadata/metadata_v2.csv", row.names = FALSE)
md <- read.csv("FAZER/Autophagy states datasets/metadata/metadata_v2.csv")

# The Autophagy genes were previous established and divided into groups that will be a parameter for the differential gene expression analysis.
ARG <- read.csv("FAZER/ARGs/segunda parte/genes.csv")
ARG$Pathway <- "Autophagy"
hgnc_genes <- unique(ARG$HGNC)

#The files will be simultaneously from the directory and the analysis will be automatic
# Paths
input_dir <- "FAZER/Autophagy states datasets/expression_matrix/tpm/Inhibition/"
output_dir <- "FAZER/Autophagy states datasets/results/Inhibition/"
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# ggplot info
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

##### Data analysis
# There will be 3 mains analysis: total genes (n3), ARG genes(n1) and ARG genes per group(n2, this analysis will be also made by GSVA). 
# The Data will be respectively named: treatment_cell, treatment_cell_ARG and treatment_cell_ARG_g.
# The raw counts will be transformed into TPM to avoid batch effects.

for (csv_file in csv_files) {
  sample_name <- tools::file_path_sans_ext(basename(csv_file))
  message("Processando: ", sample_name)
  
  # Criating new directory for each sample
  sample_outdir <- file.path(output_dir, sample_name)
  dir.create(sample_outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Organizing data tables
  data <- read.csv(csv_file)
  names(data)[names(data) == "X"] <- "Gene"
  
  data_ARG <- merge( # Analysis n1
    data.frame(Gene = unique(hgnc_genes)), 
    unique(data),
    by = "Gene",
    all.x = TRUE
  )
  rownames(data_ARG) <- data_ARG$Gene 
  data_ARG$Gene <- NULL
  data_ARG[is.na(data_ARG)] <- 0
  data_ARG$Gene_ID <- NULL
  
  data_ARG_g <- data_ARG # Analysis n2
  data_ARG_g$Gene <- rownames(data_ARG_g) 
  data_ARG_g$group <- ARG$group[match(data_ARG_g$Gene, ARG$HGNC)] 
  
  regulation <- data_ARG_g[data_ARG_g$group == "Regulation",]
  regulation$group <- NULL
  regulation$group <- regulation$Gene 
  regulation$Gene <- NULL
  regulation$hgnc_symbol <- NULL
  regulation <- regulation[complete.cases(regulation), ]
  regulation <- regulation[, c("group", setdiff(names(regulation), "group"))]
  
  prkac_rows <- grepl("^PRKAC", regulation$group)
  pka_row <- colSums(regulation[prkac_rows, sapply(regulation, is.numeric)])
  pka_row <- c(group = "PKA", pka_row)
  regulation <- regulation[!prkac_rows, ]
  regulation <- rbind(regulation, pka_row)
  
  new_rows <- grepl("^PRKA", regulation$group)
  regulation[-which(names(regulation) == "group")] <- lapply(
    regulation[-which(names(regulation) == "group")],
    function(x) suppressWarnings(as.numeric(x))
  )
  ampk_row <- colSums(regulation[new_rows, sapply(regulation, is.numeric)])
  ampk_row <- c(group = "AMPK", ampk_row)
  regulation <- regulation[!new_rows, ]
  regulation <- rbind(regulation, ampk_row)
  
  others <- data_ARG_g[data_ARG_g$group != "Regulation", ]
  numeric_cols <- sapply(others, is.numeric)
  numeric_cols["group"] <- TRUE
  others_sum <- aggregate(. ~ group, data = others[, numeric_cols], FUN = sum)
  data_ARG_g <- rbind(
    others_sum,
    regulation
  )
  rownames(data_ARG_g) <- make.unique(as.character(data_ARG_g$group))  
  data_ARG_g$group <- NULL
  data_ARG_g <- as.data.frame(data_ARG_g)
  data_ARG_g[] <- lapply(data_ARG_g, as.numeric)
  data_ARG_g$Gene_ID <- NULL
  
  data <- data[!is.na(data$Gene), ] # Analysis n3
  data <- data[!duplicated(data$Gene), ]
  rownames(data) <- data$Gene
  data$Gene <- NULL
  data$Gene_ID <- NULL
  
  # Filtering metadata
  md_filtered <- subset(md, group == sample_name)
  
  # Design
  treatments <- unique(md_filtered$treatment)
  if (length(treatments) != 2) stop(paste("Esperado exatamente dois tratamentos, mas encontrado:", paste(treatments, collapse = ", ")))
  group <- factor(md_filtered$treatment, levels = treatments)
  design <- model.matrix(~ group)
  coef_name <- colnames(design)[2] 
  
  ### Analysis N1
  log_tpm_n1 <- log2(data_ARG + 1)
  fit_n1 <- lmFit(log_tpm_n1, design)
  fit_n1 <- eBayes(fit_n1)
  res_n1 <- topTable(fit_n1, coef = coef_name, number = Inf, adjust = "BH")
  res_n1$gene <- rownames(res_n1)
  write.csv(res_n1, file.path(sample_outdir, "results_n1.csv"), row.names = FALSE)
  
  # Log2FC individual
  control_samples <- md_filtered$ID[group == "Control"]
  control_means <- rowMeans(log_tpm_n1[, control_samples, drop = FALSE])
  log2fc_indiv <- sweep(log_tpm_n1, 1, control_means, "-")
  write.csv(log2fc_indiv, file.path(sample_outdir, "n1_individual.csv"))
  
  ### Analysis N2 (gene sets)
  log_tpm_n2 <- log2(data_ARG_g + 1)
  fit_n2 <- lmFit(log_tpm_n2, design)
  fit_n2 <- eBayes(fit_n2)
  res_n2 <- topTable(fit_n2, coef = coef_name, number = Inf, adjust = "BH")
  res_n2$label <- rownames(res_n2)
  write.csv(res_n2, file.path(sample_outdir, "results_n2.csv"), row.names = FALSE)
  
  # ggplot n2
  res_n2$label <- factor(res_n2$label, levels = rev(ordered_labels))
  p <- ggplot(res_n2, aes(x = label, y = logFC, fill = label)) +
    geom_bar(stat = "identity", width = 0.6) +
    coord_flip() +
    labs(x = NULL, y = "log2FC") +
    scale_fill_manual(values = custom_colors) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  ggsave(file.path(sample_outdir, "barplot_n2.png"), plot = p, width = 6, height = 4, dpi = 300, bg = "white")
  
  ### Analysis N3 (DEGs, |logFC| > 1 and p < 0.05)
  log_tpm_n3 <- log2(data + 1)
  write.csv(log_tpm_n3, file.path(sample_outdir, "results_log2.csv"), row.names = TRUE)
  fit_n3 <- lmFit(log_tpm_n3, design)
  fit_n3 <- eBayes(fit_n3)
  res_n3 <- topTable(fit_n3, coef = coef_name, number = Inf, adjust = "BH")
  res_n3_t <- res_n3
  res_n3_t$gene <- rownames(res_n3_t)
  write.csv(res_n3_t, file.path(sample_outdir, "results_n3.csv"), row.names = TRUE)
  res_n3 <- res_n3[res_n3$P.Value < 0.05 & abs(res_n3$logFC) > 1.0, ]
  res_n3$gene <- rownames(res_n3)
  write.csv(res_n3, file.path(sample_outdir, "results_n3_filt.csv"), row.names = FALSE)
}
