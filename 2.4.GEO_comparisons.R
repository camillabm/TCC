#######################
### GEO Comparisons ###
#######################

### Explanation
# Here, the expression data from the different GEO Datasets will be compared.
# In order to try defining the autophagy cell states, there will be done a Heatmap with complete clusterization with log2FC values from GSVA analysis.
# Data will also be compared in level of common DEGs and expression of specific autophagy markers.

###### Libraries 
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(EnhancedVolcano)
library(grid)
library(limma)
library(sva)
library (VennDiagram)
library(UpSetR)
library(forcats)
library(purrr)

###### Heatmap
  # Loading data
  data <- read_excel("FAZER/Autophagy states datasets/results/Heatmap/data.xlsx")
  write.csv(data, "FAZER/Autophagy states datasets/results/Heatmap/data.csv")
  data <- read.csv("FAZER/Autophagy states datasets/results/Heatmap/data.csv")
  data$X <- NULL
  data <- data %>%
    tibble::column_to_rownames(var = "...1") #CHANGE
  
  md <- read_excel("FAZER/Autophagy states datasets/results/Heatmap/metadata.xlsx")
  write.csv(md, "FAZER/Autophagy states datasets/results/Heatmap/metadata.csv")
  md <- read.csv("FAZER/Autophagy states datasets/results/Heatmap/metadata.csv")
  md$X <- NULL
  names(md)[names(md) == "...1"] <- "SampleId"

  # Annotation
  annotation_col <- data.frame(
    Treatment = md$Tretament[match(colnames(data), md$SampleId)],
    Modulation = md$Modulation[match(colnames(data), md$SampleId)]
  )
  annotation_col$Treatment <- as.factor(annotation_col$Treatment)
  annotation_col$Modulation <- as.factor(annotation_col$Modulation)
  
  annotation_colors <- list(
    Treatment = c(
    "Rapamycin"   = "#7A4DA3",  
    "Starvation"  = "#9B59B6",  
    "Imperatorin" = "#B388EB",  
    "Z36"         = "#D6B3FF",  
    
    "Chloroquine" = "#1B9E77",  
    "Bafilomycin" = "#33A02C",  
    "Elaiophylin" = "#66C2A5",  
    "Ambroxol"    = "#A6D854",  
    "siBRDA"      = "#BCE8A3",  
    "Spautin-1"   = "#E6F5C9"   
    ),
    Modulation = c(
        "Induction"  = "#F5B7B1",  
        "Inhibition" = "#AED6F1"  
      )
    )
  
  mat_filteredA <- as.matrix(data)
  
  all(colnames(data) %in% rownames(annotation_col))  # Should return TRUE
  rownames(annotation_col) <- colnames(data)

  symbol_cols <- setdiff(colnames(md), c("SampleId", "Tretament", "Modulation"))
  md[symbol_cols] <- lapply(md[symbol_cols], function(col) {
    col[is.na(col)] <- ""
    return(col)
  })
  symbols_matrix <- md[, c("SampleId", symbol_cols)] %>%
    pivot_longer(cols = all_of(symbol_cols), names_to = "Gene", values_to = "Symbol") %>%
    pivot_wider(names_from = SampleId, values_from = Symbol) %>%
    column_to_rownames("Gene")
  symbols_matrix <- symbols_matrix[rownames(mat_filteredA), colnames(mat_filteredA)]
  
  # Heatmap plotting
  png("FAZER/Autophagy states datasets/results/Heatmap/Heatmap.png", width = 10, height = 7, units = "in", res = 1000)
  heatmap <- pheatmap(mat_filteredA,
                      scale = "row",
                      cluster_rows = TRUE,         
                      cluster_cols = TRUE,
                      treeheight_row = 10,           
                      treeheight_col = 10,
                      color = colorRampPalette(c("blue","white", "red"))(100),
                      display_numbers = as.matrix(symbols_matrix),
                      number_color = "black", 
                      show_rownames = TRUE,
                      show_colnames = TRUE,
                      annotation_col = annotation_col,
                      annotation_colors = annotation_colors,
                      legend = TRUE,
                      cellwidth = 15,
                      cellheight = 15,
                      fontsize = 10,
                      border_color = "black"
  )
  grid.text("* - p-value < 0.05", x = 0.9, y = 0.2, gp = gpar(fontsize = 10, col = "black"))
  grid.text("** - p-value < 0.01", x = 0.9, y = 0.17, gp = gpar(fontsize = 10, col = "black"))
  grid.text("*** - p-value < 0.001", x = 0.9, y = 0.14, gp = gpar(fontsize = 10, col = "black"))
  grid.text("# - p-value < 0.0001", x = 0.9, y = 0.11, gp = gpar(fontsize = 10, col = "black"))
  dev.off()

###### Commom overall DEGs
  # Data
  main_dir <- file.path("FAZER", "Autophagy states datasets", "results", "Inhibition")
  file_paths <- list.files(path = main_dir, 
                           pattern = "results_n3_filt\\.csv$", 
                           recursive = TRUE, 
                           full.names = TRUE)
  if (length(file_paths) == 0) {
    stop("Nenhum arquivo encontrado. Verifique o caminho e o padrão de busca.")
  } else if (length(file_paths) < 2) {
    warning("Apenas um arquivo encontrado. O gráfico UpSet requer múltiplos conjuntos.")
  }
  
  # Genes
  gene_lists <- lapply(file_paths, function(file) {
    tryCatch({
      df <- read.csv(file, stringsAsFactors = FALSE)
      if (!"gene" %in% colnames(df)) {
        stop(paste("Arquivo sem coluna 'gene':", file))
      }
      unique(trimws(as.character(df$gene)))  # remove espaços e garante caracteres
    }, error = function(e) {
      stop(paste("Erro ao processar arquivo", file, ":", e$message))
    })
  })
  names(gene_lists) <- make.names(basename(dirname(file_paths)))
  all_genes <- unique(unlist(gene_lists))
  
  # Matrix
  binary_matrix <- do.call(cbind, lapply(gene_lists, function(glist) {
    as.integer(all_genes %in% glist)
  }))
  rownames(binary_matrix) <- all_genes
  binary_df <- as.data.frame(binary_matrix)
  cat("\n=== Resumo dos dados ===\n")
  cat("Número de genes únicos:", length(all_genes), "\n")
  cat("Número de conjuntos:", length(gene_lists), "\n")
  cat("Tamanho de cada conjunto:\n")
  print(sapply(gene_lists, length))
  filtered_genes <- binary_df[rowSums(binary_df) >= 2, ]
  
  # UpSet plot
  png("FAZER/Autophagy states datasets/results/Intersections/inhibitors_upset_plot.png", width=2000, height=2100, res=300)
  upset(
    filtered_genes,
    nsets = ncol(filtered_genes),
    nintersects = 30,
    order.by = "degree",
    sets.bar.color = "steelblue"
  )
  dev.off()
  
  # Intersections
  all_combos <- list()
  
  for (i in 2:ncol(filtered_genes)) {
    combos <- combn(colnames(filtered_genes), i, simplify = FALSE)
    for (combo in combos) {
      genes <- rownames(filtered_genes)[
        rowSums(filtered_genes[, combo, drop = FALSE] == 1) == length(combo) &
          rowSums(filtered_genes[, !colnames(filtered_genes) %in% combo, drop = FALSE]) == 0
      ]
      if (length(genes) > 0) {
        nome <- paste(combo, collapse = "_AND_")
        all_combos[[nome]] <- genes
      }
    }
  }
  
  info <- data.frame(
    name = names(all_combos),
    count = sapply(all_combos, length),  # Quantidade de genes na interseção
    degree = sapply(names(all_combos), function(x) length(strsplit(x, "_AND_")[[1]])),  # Número de conjuntos
    stringsAsFactors = FALSE
  )

  top <- head(info[order(-info$degree), ], 30)
  
  result_df <- data.frame(
    intersection = top$name,
    genes = sapply(top$name, function(name) paste(all_combos[[name]], collapse = ", "))
  )
  
  write.csv(result_df, "FAZER/Autophagy states datasets/results/Intersections/intersections_inh.csv", row.names = FALSE)
  
###### Common Genes
  # Gene defining
  fats <- c("TFEB","ZKSCAN2","FOXO1","FOXO3","GATA1","GATA3")
  reg <- c("RPTOR","DEPTOR")
  atgs <- c("ATG1","BECN1","ATG7","MAP1LC3B","MAP1LC3B2","MAP1LC3C","MAP1LC3A","GABARAPL1","ATG9A","ATG9B","ATG12","ATG13","ATG14","WIPI1","SQSTM1")
  
  # Directory
  main_dir <- "FAZER/Autophagy states datasets/results/Inhibition/"

  # log2FC
  file_paths <- list.files(path = main_dir, pattern = "results_n3\\.csv$", recursive = TRUE, full.names = TRUE)
  result_list <- list()

  for (file in file_paths) {
    df <- read.csv(file, stringsAsFactors = FALSE)
    experiment_name <- basename(dirname(file))
    selected_genes <- df[df$gene %in% c(fats, reg, atgs), ]
    selected_genes$Class <- NA
    selected_genes$Class[selected_genes$gene %in% fats] <- "Transcription factors"
    selected_genes$Class[selected_genes$gene %in% reg] <- "Regulatory genes"
    selected_genes$Class[selected_genes$gene %in% atgs] <- "ATGs"
    selected_genes$Values <- paste(selected_genes$logFC, selected_genes$P.Value, sep = ", ")
    selected_genes$Experiment <- experiment_name
    result_list[[experiment_name]] <- selected_genes
  }
  final_result <- do.call(rbind, result_list)
  final_result_wide <- final_result %>%
    select(Class, gene, Experiment, Values) %>%
    pivot_wider(names_from = Experiment, values_from = Values)

  exp_cols <- setdiff(colnames(final_result_wide), c("Class", "gene"))
  
  for (col in exp_cols) {
    final_result_wide <- final_result_wide %>%
      separate(
        col,
        into = c(paste0(col, "_logFC"), paste0(col, "_P.Value")),
        sep = ",\\s*",
        extra = "merge",
        fill = "right",
        convert = TRUE
      )
  }
  
  # for (col in exp_cols) {
  #   logfc_col <- paste0(col, "_logFC")
  #   pval_col <- paste0(col, "_P.Value")
  #   
  #   final_result_wide[[logfc_col]][final_result_wide[[pval_col]] > 0.05] <- NA
  #   final_result_wide[[pval_col]][final_result_wide[[pval_col]] > 0.05] <- NA
  # }
  final_result_wide <- final_result_wide %>%
    arrange(Class, gene)

  write.csv(final_result_wide, "FAZER/Autophagy states datasets/results/specific_genes/genes_inducers_log2fc.csv", row.names = FALSE)

  # log2
  genes_interesse <- c(fats, reg, atgs)
  
  file_paths <- list.files(path = main_dir, pattern = "results_log2\\.csv$", recursive = TRUE, full.names = TRUE)
  
  result_list <- list()
  
  for (file in file_paths) {
    experiment_name <- basename(dirname(file))
    df <- read.csv(file, stringsAsFactors = FALSE)

    if ("X" %in% colnames(df)) {
      df <- df %>% rename(gene = X)
    }

    df <- df %>%
      filter(gene %in% genes_interesse) %>%
      select(-starts_with("C"))

    df <- df %>%
      mutate(
        Class = case_when(
          gene %in% fats ~ "Transcription factors",
          gene %in% reg  ~ "Regulatory genes",
          gene %in% atgs ~ "ATGs",
          TRUE ~ NA_character_
        )
      )

    df_renamed <- df %>%
      select(-Class) %>%
      rename_with(~ paste0(experiment_name, "_", .x), -gene)

    df_final <- df_renamed %>%
      left_join(df %>% select(gene, Class) %>% distinct(), by = "gene")
    
    result_list[[experiment_name]] <- df_final
  }

  final_result <- reduce(result_list, full_join, by = c("gene", "Class")) %>%
    relocate(Class, gene) %>%
    arrange(Class, gene)
  
  write.csv(final_result, "FAZER/Autophagy states datasets/results/specific_genes/genes_log2_inhibitors.csv", row.names = FALSE)