############
### UMAP ###
############

###### Libraries
library(sva)
library(umap)
library(uwot)
library(readxl)

###### Data
gbm <- read.csv("FAZER/Expression/results/matrix_gbmR_comb/GSVA_results_summary.csv")
rownames(gbm) <- gbm$X
gbm$X <- NULL

brca <- read.csv("FAZER/Expression/results/matrix_brcaR_comb/GSVA_results_summary.csv")
rownames(brca) <- brca$X
brca$X <- NULL

luad <- read.csv("FAZER/Expression/results/matrix_luadR_comb/GSVA_results_summary.csv")
rownames(luad) <- luad$X
luad$X <- NULL

coad <- read.csv("FAZER/Expression/results/matrix_coadR_comb/GSVA_results_summary.csv")
rownames(coad) <- coad$X
coad$X <- NULL

data <- cbind(gbm, brca, luad, coad)

###### Metadata
MD <- read_excel("FAZER/Expression/results/survival/umap.xlsx")
write.csv(MD, "FAZER/Expression/results/survival/umap.csv")
MD <- read.csv("FAZER/Expression/results/survival/umap.csv")
MD$X <- NULL

###### UMAP
data <- ComBat(dat = gbm, batch = MD$state)

set.seed(125)

umap_input <- t(gbm)

umap_result <- umap(
  umap_input,
  n_neighbors = 15,
  min_dist = 0.001,
  spread = 10,
  metric = "cosine",        # ou "euclidean", se preferir
  n_components = 2,
  learning_rate = 1,
  fast_sgd = TRUE,
  verbose = TRUE
)

umap_df <- data.frame(
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2],
  Project = MD$project,
  State = MD$state
)

colors <- c(
  "1" = "gray",
  "2" = "red")

shape_values <- c("GBM" = 16)

umap_df$State <- as.factor(umap_df$State)

UMAP <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = State)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shape_values) +
  theme_classic(base_family = "Arial") +
  theme(
    text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8),  # mesmo valor da borda
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    plot.title = element_blank()
  ) +
  labs(
    x = "UMAP_1",
    y = "UMAP_2",
    color = "State",
    shape = "Project"
  )
