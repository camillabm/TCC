###############
### ARG GBM ###
###############

### Explanation
# We will use the previous TPM converted GEO datasets
# STEPS:
# (1) Filtering for ARGs
# (2) Analysis of ARGs 
# (3) Analysis of gene sets with GSVA for profiles definitions
# (4) Data plots as volcanos 
# (5) Data separation for differential clinical conditions

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
input_file <- "FAZER/GBM/LGG/lgg_exp.csv" #CHANGE
output_dir <- "FAZER/GBM/Results/LGG" #CHANGE
arg_file <- "FAZER/ARGs/segunda parte/genes.csv"
md <- read.csv("FAZER/GBM/LGG/lgg_md.csv") #CHANGE

###### ARGs
ARG <- read.csv(arg_file)
ARG$Pathway <- "Autophagy"
hgnc_genes <- unique(ARG$HGNC)

###### Analysis
csv_file <- read.csv(input_file)
csv_file_fil <- csv_file[csv_file$X %in% hgnc_genes, ]
rownames(csv_file_fil) <- csv_file_fil$X
csv_file_fil$X <- NULL
cols_to_change <- grepl("^TCGA", colnames(csv_file_fil))
colnames(csv_file_fil)[cols_to_change] <- gsub("\\.", "-", colnames(csv_file_fil)[cols_to_change])

group <- md$condition

design <- model.matrix(~ group)
coef_name <- colnames(design)[2]

# GSVA
RAPA <- csv_file_fil
RAPA <- cbind(Gene = rownames(RAPA), RAPA)
rownames(RAPA) <- NULL
RAPA <- merge(data.frame(Gene = unique(hgnc_genes)), RAPA, by = "Gene", all.x = TRUE)
rownames(RAPA) <- RAPA$Gene
RAPA$Gene <- NULL
RAPA[is.na(RAPA)] <- 0

RAPA_log2 <- log2(RAPA + 1) 
data <- RAPA_log2

data <- as.matrix(data)
gene_var <- apply(data, 1, var)
data_filtered <- data[gene_var > 0, ]

group_gsva <- ifelse(md$condition == "Tumor", 1, 0)
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
write.csv(gsva_es, file = file.path(output_dir, "GSVA_results.csv"))

# Individual
control_samples <- md$barcode[md$condition == "Normal"]
log_tpm_control <- gsva_es[, colnames(gsva_es) %in% control_samples, drop = FALSE]
control_means <- rowMeans(log_tpm_control, na.rm = TRUE)
log2fc_indiv_gsva <- sweep(gsva_es, 1, control_means, "-")
write.csv(log2fc_indiv_gsva, file.path(output_dir, "GSVA_results_individual.csv"))

# Grouped
## G3
G3_barcodes <- md %>%
  filter(paper_Grade == "G3") %>%
  pull(barcode)
barcodes_gtex <- md %>%
  filter(project_id == "GTEx") %>%
  pull(barcode)
grouped_g3 <- c(barcodes_gtex,G3_barcodes)

md_g3 <- md[md$barcode %in% grouped_g3,]
gsva_es_g3 <- gsva_es[, colnames(gsva_es) %in% grouped_g3]
group_gsva_g3 <- ifelse(md_g3$condition == "Tumor", 1, 0)

design <- model.matrix(~ group_gsva_g3)
fit <- lmFit(gsva_es_g3, design)
fit <- eBayes(fit)
results_g3 <- topTable(fit, coef = "group_gsva_g3", number = Inf)
results_g3$label <- rownames(results_g3)
write.csv(results_g3, file.path(output_dir, "GSVA_tot_g3.csv"))

## G2
G2_barcodes <- md %>%
  filter(paper_Grade == "G2") %>%
  pull(barcode)
barcodes_gtex <- md %>%
  filter(project_id == "GTEx") %>%
  pull(barcode)
grouped_g2 <- c(barcodes_gtex,G2_barcodes)

md_g2 <- md[md$barcode %in% grouped_g2,]
gsva_es_g2 <- gsva_es[, colnames(gsva_es) %in% grouped_g2]
group_gsva_g2 <- ifelse(md_g2$condition == "Tumor", 1, 0)

design <- model.matrix(~ group_gsva_g2)
fit <- lmFit(gsva_es_g2, design)
fit <- eBayes(fit)
results_g2 <- topTable(fit, coef = "group_gsva_g2", number = Inf)
results_g2$label <- rownames(results_g2)
write.csv(results_g2, file.path(output_dir, "GSVA_tot_g2.csv"))

# DEG analysis
log_tpm_n1 <- log2(csv_file_fil + 1)

## G3
log_tpm_n1_g3 <- log_tpm_n1[, colnames(log_tpm_n1) %in% grouped_g3]
group_log_tpm_n1_g3 <- ifelse(md_g3$condition == "Tumor", 1, 0)

design <- model.matrix(~ group_log_tpm_n1_g3)
fit <- lmFit(log_tpm_n1_g3, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "group_log_tpm_n1_g3", number = Inf)
res_n1_g3 <- results
res_n1_g3$label <- rownames(res_n1_g3)
write.csv(res_n1_g3, file.path(output_dir, "group_log_tpm_n1_g3.csv"))

## G2
log_tpm_n1_g2 <- log_tpm_n1[, colnames(log_tpm_n1) %in% grouped_g2]
group_log_tpm_n1_g2 <- ifelse(md_g2$condition == "Tumor", 1, 0)

design <- model.matrix(~ group_log_tpm_n1_g2)
fit <- lmFit(log_tpm_n1_g2, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "group_log_tpm_n1_g2", number = Inf)
res_n1_g2 <- results
res_n1_g2$label <- rownames(res_n1_g2)
write.csv(res_n1_g3, file.path(output_dir, "group_log_tpm_n1_g2.csv"))


###### Data plots
res_n1 <- res_n1_g2
results <- results_g2

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
                                #lab = ifelse(combined_data$source == "GSVA", combined_data$gene, ""),
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

png(file.path(output_dir, "volcano_g2_correct.png"), width = 3000, height = 2000, res = 300)


grid.arrange(volcano_grob, legend_grob,
             ncol = 2,
             widths = c(4, 1))  

dev.off()

###### Data saving for differential clinical conditions
split_gsva_by_group <- function(metadata, gsva, group_col, group_levels = NULL) {
  library(dplyr)
  
  if (!is.null(group_levels)) {
    metadata <- metadata %>%
      mutate(!!group_col := factor(.data[[group_col]], levels = group_levels))
  }

  groups <- if (!is.null(group_levels)) group_levels else sort(unique(metadata[[group_col]]))
  
  result <- list()
  
  for (grp in groups) {
    md_sub <- metadata %>% filter(.data[[group_col]] == grp)
    barcodes_grp <- md_sub$barcode
    barcodes_grp <- intersect(barcodes_grp, colnames(gsva))
    gsva_sub <- gsva[, barcodes_grp, drop = FALSE]
    result[[as.character(grp)]] <- list(
      metadata = md_sub,
      gsva = gsva_sub
    )
  }
  
  return(result)
}

res_by_grade <- split_gsva_by_group(
  metadata = md,
  gsva = log2fc_indiv_gsva,
  group_col = "paper_Grade", #CHANGE
  group_levels = c("G2", "G3")
)

for (grp in names(res_by_grade)) {
  write.csv(res_by_grade[[grp]]$metadata, file = file.path(output_dir, paste0("metadata_", grp, ".csv")), row.names = FALSE)
  write.csv(res_by_grade[[grp]]$gsva, file = file.path(output_dir, paste0("gsva_", grp, ".csv")), row.names = TRUE)
}
