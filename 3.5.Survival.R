#####################
### TCGA Survival ###
#####################

### Explanation
# Autophagy states were defined by state plot of AutoIndex and LysoIndex. 
# The inhibition samples were selected to match the number of induction samples.

###### Libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)

##### Paths and data
md_dir <- "FAZER/Expression/tables/metadata/cohort"
clin_dir <- "FAZER/Clinical data"

###### Clinical data binding
md_files <- list.files(md_dir, full.names = TRUE, pattern = "\\.csv$")
clin_files <- list.files(clin_dir, full.names = TRUE, pattern = "\\.csv$")

get_md_code <- function(path) str_sub(basename(path), 10, 13)
get_clin_code <- function(path) str_sub(basename(path), 15, 18)

file_pairs <- expand.grid(md = md_files, clin = clin_files, stringsAsFactors = FALSE) %>%
  mutate(md_code = get_md_code(md),
         clin_code = get_clin_code(clin)) %>%
  filter(md_code == clin_code)

metadata_list <- list()

for (i in seq_len(nrow(file_pairs))) {
  md_path <- file_pairs$md[i]
  clin_path <- file_pairs$clin[i]
  code <- file_pairs$md_code[i]
  
  md <- read_csv(md_path, show_col_types = FALSE)
  clin <- read_csv(clin_path, show_col_types = FALSE)
  
  if (!"barcode_12" %in% names(md)) {
    md <- md %>%
      mutate(barcode_12 = str_sub(barcode, 1, 12))
  }
  
  merged <- left_join(md, clin, by = c("barcode_12" = "submitter_id"))
  if (!inherits(merged, "data.frame")) {
    warning(paste("Objeto 'merged' não é um data.frame para:", code))
    next
  }
  
  cols_to_keep <- c("condition", "project_id", "tissue_or_organ_of_origin", "barcode", "days_to_death",
                    "treatments_pharmaceutical_treatment_or_therapy",
                    "treatments_radiation_treatment_or_therapy")
  
  merged_clean <- dplyr::select(merged, any_of(cols_to_keep))
  
  metadata_list[[code]] <- merged_clean
}

####### States
gbm <- metadata_list[["gbmR"]]
brca <- metadata_list[["brca"]]
luad <- metadata_list[["luad"]]
coad <- metadata_list[["coad"]]

states_xlsx <- read_excel("FAZER/Expression/results/Survival/luad_all.xlsx") #CHANGE
write.csv(states_xlsx, "FAZER/Expression/results/Survival/luad_all.csv", row.names = FALSE)
states_df <- read.csv("FAZER/Expression/results/Survival/luad_all.csv", stringsAsFactors = FALSE)

states_df$barcode <- ifelse(
  grepl("^TCGA", states_df$barcode),
  gsub("\\.", "-", states_df$barcode),
  states_df$barcode
)

adiciona_state <- function(metadata, states_df) {
  
  metadata$barcode <- as.character(metadata$barcode)
  states_df$barcode <- as.character(states_df$barcode)
  
  metadata <- merge(metadata, states_df[, c("barcode", "A", "B", "c_status","m_status","state")],
                    by = "barcode", all.x = TRUE)

  
  metadata$Inhibitors <- NULL
  metadata$Inducers <- NULL
  
  return(metadata)
}

dataset <- adiciona_state(luad, states_df) #CHANGE

###### Analysis
dataset$Status_OS <- ifelse(is.na(dataset$days_to_death), 0, 1)

dataset$terc <- ifelse(
  dataset$A == 2, "a",
  ifelse(dataset$B == 1, "c", "b")
)

write.csv(dataset, "FAZER/Expression/results/survival/luad_ready.csv")