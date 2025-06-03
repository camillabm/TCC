############
### GSVA ###
############

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("dplyr")
BiocManager::install(c("GSVA","limma","GSEABase"))
BiocManager::install("BiocParallel")
BiocManager::install("BiocParallel")
writeLines(
  'PATH="C:\\rtools43\\usr\\bin;C:\\rtools43\\x86_64-w64-mingw32.static.posix\\bin;${PATH}"',
  con = "~/.Renviron"
)
Sys.which("g++")

###### Libraries
library(GSVA)
library(limma)
library(dplyr)
library(ggplot2)
library(GSEABase)

###### Paths
input_dir <- "FAZER/Autophagy states datasets/expression_matrix/tpm/Induction/"
output_dir <- "FAZER/Autophagy states datasets/results/Induction/"
metadata_file <- "FAZER/Autophagy states datasets/metadata/metadata_v2.csv"
arg_file <- "FAZER/ARGs/segunda parte/genes.csv"
results_file <- file.path(output_dir, "inhibition_GSVA_results.csv")

###### Data loading
md <- read.csv(metadata_file)
ARG <- read.csv(arg_file)
ARG$Pathway <- "Autophagy"
ensembl_genes <- ARG$ensembl_gene_id
hgnc_genes <- ARG$HGNC

csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

###### GSVA analysis
for (csv_file in csv_files) {
  sample_name <- tools::file_path_sans_ext(basename(csv_file))
  
  # Data reading
  RAPA <- read.csv(csv_file)
  names(RAPA)[names(RAPA) == "X"] <- "Gene"
  RAPA <- merge(
    data.frame(Gene = unique(hgnc_genes)), 
    unique(RAPA),
    by = "Gene",
    all.x = TRUE
  )
  rownames(RAPA) <- RAPA$Gene 
  RAPA$Gene <- NULL
  RAPA[is.na(RAPA)] <- 0
  RAPA$Gene_ID <- NULL
  
  RAPA_log2 <- log2(RAPA + 1) 
  data <- RAPA
  data <- as.matrix(data)
  
  gene_var <- apply(data, 1, var)
  data_filtered <- data[gene_var > 0, ]
  
  samples <- colnames(data)
  group <- ifelse(grepl("^C", samples), 0, 1)
  order <- c(grep("^C", samples, value = TRUE), grep("^[^C]", samples, value = TRUE))
  data <- data[, order]
  
  # Subgroups
  initiation <- ARG %>% filter(group == "Initiation") %>% pull(HGNC)
  nucleation <- ARG %>% filter(group == "Nucleation") %>% pull(HGNC)
  elongation <- ARG %>% filter(group == "Elongation") %>% pull(HGNC)
  maturation <- ARG %>% filter(group == "Maturation") %>% pull(HGNC)
  enzimes <- ARG %>% filter(group == "Lysosomal enzimes") %>% pull(HGNC)
  membrane <- ARG %>% filter(group == "Lysosomal membrane") %>% pull(HGNC)
  transporters <- ARG %>% filter(group == "Lysosomal transporters") %>% pull(HGNC)
  reg_pos <- ARG %>% filter(influence.in.autophagy == "positive") %>% pull(HGNC)
  reg_neg <- ARG %>% filter(influence.in.autophagy == "negative") %>% pull(HGNC)
  
  autophagy <- c(initiation, nucleation, elongation, maturation)
  lisosomal_activity <- c(enzimes, membrane, transporters)
  
  autophagy_genes <- list(
    Autophagy = autophagy,
    Lisosomal_activity = lisosomal_activity,
    Initiation = initiation,
    Nucleation = nucleation,
    Elongation = elongation,
    Maturation = maturation,
    Enzymes = enzimes,
    Membrane = membrane,
    Transporters = transporters,
    Reg_pos = reg_pos,
    Reg_neg = reg_neg
  )
  
  # GSVA
  gsvapar <- gsvaParam(data, autophagy_genes)
  gsva_es <- suppressMessages(suppressWarnings(gsva(gsvapar)))
  write.csv(gsva_es, file = file.path(output_dir, sample_name, "GSVA_results.csv"))
  
  order_labels <- c("Autophagy", "Initiation", "Nucleation", "Elongation", "Maturation",
                    "Lisosomal_activity", "Transporters", "Membrane", "Enzymes",
                    "Reg_pos", "Reg_neg")
  
    #Individual
  control_samples <- md$ID[md$treatment == "Control"]
  log_tpm_control <- gsva_es[, colnames(gsva_es) %in% control_samples, drop = FALSE]
  control_means <- rowMeans(log_tpm_control, na.rm = TRUE)
  log2fc_indiv <- sweep(gsva_es, 1, control_means, "-")
  log2fc_indiv<- log2fc_indiv[match(order_labels, rownames(log2fc_indiv)), , drop = FALSE]
  sample_outdir <- file.path(output_dir, sample_name)
  dir.create(sample_outdir, recursive = TRUE, showWarnings = FALSE)
  write.csv(log2fc_indiv, file = file.path(sample_outdir, "GSVA_individual.csv"))
  
    #Grouped
  design <- cbind(Intercept = 1, Treated_vs_Control = group)
  fit <- lmFit(gsva_es, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "Treated_vs_Control", number = Inf)
  results$label <- rownames(results)
  results$label <- factor(results$label, levels = order_labels)
  results_1 <- results[match(order_labels, as.character(results$label)), ]

  results$label <- factor(results$label, levels = rev(order_labels))

  custom_colors <- c(
    "Autophagy" = "#B39DDB",
    "Initiation" = "#B39DDB",
    "Nucleation" = "#B39DDB",
    "Elongation" = "#B39DDB",
    "Maturation" = "#B39DDB",
    "Lisosomal_activity" = "#AED581",
    "Transporters" = "#AED581",
    "Membrane" = "#AED581",
    "Enzymes" = "#AED581",
    "Reg_pos" = "#EF9A9A",
    "Reg_neg" = "#90CAF9"
  )

  p <- ggplot(results, aes(x = label, y = logFC, fill = label)) +
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

  # Saving plot
  ggsave(
    filename = file.path(output_dir, sample_name, "barplot_GSVA.png"),
    plot = p,
    width = 6, height = 4, dpi = 300, bg = "white"
  )

  # Saving table
  cat(paste0("=====  ", sample_name, "  =====\n"), file = results_file, append = TRUE)
  write.table(
    results_1,
    file = results_file,
    sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE,
    append = TRUE
  )
  cat("\n", file = results_file, append = TRUE)
}

