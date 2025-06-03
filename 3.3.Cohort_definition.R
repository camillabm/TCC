#########################
### Cohort definition ###
#########################

### Warning!
# We will use a N = 367 for both normal and tumor samples for each tumor type (maximum of GBM Tumor samples).
# Therefore, every TCGA normal sample will be depleted and the 367 normal samples will come only from GTEx projects from same tissue type.
# When there are more than 367 samples available, they will be randomly selected from available samples. 

# Libraries
library(SummarizedExperiment)

### GBM ###
gbm_cohort <- read.csv("FAZER/Expression/tables/rna/others/expression_gbm.csv", row.names = 1) #Previous expression table (genes x samples)
md_gbm_cohort <- read.csv("FAZER/Expression/tables/metadata/others/metadata_gbmR_f.csv") #Metadata table after manual curation (in .csv) (sample x metadata)

valid_samples <- md_gbm_cohort$barcode #Identification of samples in metadata after manual curation
colnames(gbm_cohort) <- as.character(gsub("\\.", "-", colnames(gbm_cohort))) #Sample ID padronization
gbm_cohort <- gbm_cohort[, colnames(gbm_cohort) %in% valid_samples] #Expression table filtered

gbm_normal_cohort <- read.csv("FAZER/Expression/tables/rna/others/matrix_gbmR_normal.csv", row.names = 1)
md_gbm_normal_cohort <- read.csv("FAZER/Expression/tables/metadata/others/metadata_gbmR_normal.csv")

  # Normal samples selection
cortex_samples <- 184 
frontal_cortex_samples <- 183 #This will provide 367 normal samples 

md_gbm_cortex <- md_gbm_normal_cohort[md_gbm_normal_cohort$tissue_or_organ_of_origin == "Brain cortex", ]
cortex_samples_sel <- sample(md_gbm_cortex$barcode, cortex_samples)
matrix_cortex_sel <- gbm_normal_cohort[, colnames(gbm_normal_cohort) %in% cortex_samples_sel]
md_cortex_sel <- md_gbm_normal_cohort[md_gbm_normal_cohort$barcode %in% colnames(matrix_cortex_sel), ]

md_gbm_frontal <- md_gbm_normal_cohort[md_gbm_normal_cohort$tissue_or_organ_of_origin == "Brain frontal cortex", ]
frontal_cortex_samples_sel <- sample(md_gbm_frontal$barcode, frontal_cortex_samples)
matrix_frontal_sel <- gbm_normal_cohort[, colnames(gbm_normal_cohort) %in% frontal_cortex_samples_sel]
md_frontal_sel <- md_gbm_normal_cohort[md_gbm_normal_cohort$barcode %in% colnames(matrix_frontal_sel), ]

matrix_comb <- cbind(matrix_cortex_sel, matrix_frontal_sel)
metadata_comb <- rbind(md_cortex_sel, md_frontal_sel)
metadata_comb$X <- NULL

  # Tumor samples selection
normal_del <- md_gbm_cohort$barcode[md_gbm_cohort$condition == "Normal"]
gbm_cohort <- gbm_cohort[, !(colnames(gbm_cohort) %in% normal_del)] #Should rest only 367 variables

md_gbm_cohort <- md_gbm_cohort[, c("barcode", "tissue_or_organ_of_origin", "condition", "project_id")] #Should rest only 4 variables
md_gbm_cohort <- md_gbm_cohort[md_gbm_cohort$barcode %in% colnames(gbm_cohort), ] #Should rest only 367 variables

 # Data combining
matrix_comb <- SummarizedExperiment(
  assays = list(counts = as.matrix(matrix_comb)),
  rowData = DataFrame(geneID = rownames(matrix_comb))
)
matrix_comb <- assay(matrix_comb)
rownames(matrix_comb) <- gsub("\\..*$", "", rownames(matrix_comb))
matrix_comb <- matrix_comb[!duplicated(rownames(matrix_comb)), ] #Remove double lines

genes_normal <- rownames(matrix_comb)
genes_tumor <- rownames(gbm_cohort)
genes_comum <- intersect(genes_normal, genes_tumor)

gbm_cohort <- gbm_cohort[rownames(gbm_cohort) %in% genes_comum, ] #Tumor
matrix_comb <- matrix_comb[rownames(matrix_comb) %in% genes_comum, ] #Normal

   # Filtering by autophagy genes 
autophagy_genes <- read.csv("FAZER/ARGs/segunda parte/genes.csv") 
head(autophagy_genes)
ensembl_genes <- autophagy_genes$ensembl_gene_id  #Columns in gene table with ensembl ID

gbm_cohort  <- gbm_cohort [rownames(gbm_cohort) %in% ensembl_genes, ] 
matrix_comb  <- matrix_comb[rownames(matrix_comb) %in% ensembl_genes, ] 

  # Binding tables
matrix_comb <- matrix_comb[rownames(gbm_cohort), ] #Reorder the rownames to match
matrix_gbmR <- cbind(matrix_comb, gbm_cohort)
metadata_gbmR <- rbind(metadata_comb, md_gbm_cohort)

  #Saving tables
write.csv(matrix_gbmR, "FAZER/Expression/tables/rna/cohort/matrix_gbmR_comb.csv", row.names = TRUE) #expression matrix
write.csv(metadata_gbmR, "FAZER/Expression/tables/metadata/cohort/metadata_gbmR_comb.csv", row.names = TRUE) #metadata

### BRCA ###
brca_cohort <- read.csv("FAZER/Expression/tables/rna/others/expression_brca.csv", row.names = 1) #Previous expression table (genes x samples)
md_brca_cohort <- read.csv("FAZER/Expression/tables/metadata/others/metadata_brcaR_f.csv") #Metadata table after manual curation (in .csv) (sample x metadata)

valid_samples <- md_brca_cohort$barcode #Identification of samples in metadata after manual curation
colnames(brca_cohort) <- as.character(gsub("\\.", "-", colnames(brca_cohort))) #Sample ID padronization
brca_cohort <- brca_cohort[, colnames(brca_cohort) %in% valid_samples] #Expression table filtered

brca_normal_cohort <- read.csv("FAZER/Expression/tables/rna/others/matrix_brcaR_normal.csv", row.names = 1)
md_brca_normal_cohort <- read.csv("FAZER/Expression/tables/metadata/others/metadata_brcaR_normal.csv")

  # Normal samples selection
breast_samples <- 367 #This will provide 367 samples 

breast_samples_sel <- sample(md_brca_normal_cohort$barcode, breast_samples)
matrix_breast_sel <- brca_normal_cohort[, colnames(brca_normal_cohort) %in% breast_samples_sel]
md_breast_sel <- md_brca_normal_cohort[md_brca_normal_cohort$barcode %in% colnames(matrix_breast_sel), ]
md_breast_sel$X <- NULL

  # Tumor samples selection
normal_del <- md_brca_cohort$barcode[md_brca_cohort$condition == "Normal"]
brca_cohort <- brca_cohort[, !(colnames(brca_cohort) %in% normal_del)] 

md_brca_cohort <- md_brca_cohort[, c("barcode", "tissue_or_organ_of_origin", "condition", "project_id")] #Should rest only 4 variables
md_brca_cohort <- md_brca_cohort[md_brca_cohort$barcode %in% colnames(brca_cohort), ] 

brca_samples_sel <- sample(md_brca_cohort$barcode, breast_samples)
brca_cohort <- brca_cohort[, colnames(brca_cohort) %in% brca_samples_sel]
md_brca_cohort_sel <- md_brca_cohort[md_brca_cohort$barcode %in% colnames(brca_cohort), ] #Should rest only 367 variables

  # Data combining
matrix_breast_sel <- SummarizedExperiment(
  assays = list(counts = as.matrix(matrix_breast_sel)),
  rowData = DataFrame(geneID = rownames(matrix_breast_sel))
)
matrix_breast_sel <- assay(matrix_breast_sel)
rownames(matrix_breast_sel) <- gsub("\\..*$", "", rownames(matrix_breast_sel))
matrix_breast_sel <- matrix_breast_sel[!duplicated(rownames(matrix_breast_sel)), ] #Remove double lines

genes_normal <- rownames(matrix_breast_sel)
genes_tumor <- rownames(brca_cohort)
genes_comum <- intersect(genes_normal, genes_tumor)

brca_cohort <- brca_cohort[rownames(brca_cohort) %in% genes_comum, ] #Tumor
matrix_breast_sel <- matrix_breast_sel[rownames(matrix_breast_sel) %in% genes_comum, ] #Normal

  # Filtering by autophagy genes 
autophagy_genes <- read.csv("FAZER/ARGs/segunda parte/genes.csv") 
head(autophagy_genes)
ensembl_genes <- autophagy_genes$ensembl_gene_id  #Columns in gene table with ensembl ID

brca_cohort  <- brca_cohort[rownames(brca_cohort) %in% ensembl_genes, ] 
matrix_breast_sel  <- matrix_breast_sel[rownames(matrix_breast_sel) %in% ensembl_genes, ] 

  # Binding tables
matrix_breast_sel <- matrix_breast_sel[rownames(brca_cohort), ] #Reorder the rawnames to match
matrix_brcaR <- cbind(matrix_breast_sel, brca_cohort)
metadata_brcaR <- rbind(md_breast_sel, md_brca_cohort_sel)

  #Saving tables
write.csv(matrix_brcaR, "FAZER/Expression/tables/rna/cohort/matrix_brcaR_comb.csv", row.names = TRUE) #expression matrix
write.csv(metadata_brcaR, "FAZER/Expression/tables/metadata/cohort/metadata_brcaR_comb.csv", row.names = TRUE) #metadata

### LUAD ###
luad_cohort <- read.csv("FAZER/Expression/tables/rna/others/expression_luad.csv", row.names = 1) #Previous expression table (genes x samples)
md_luad_cohort <- read.csv("FAZER/Expression/tables/metadata/others/metadata_luadR_f.csv") #Metadata table after manual curation (in .csv) (sample x metadata)

valid_samples <- md_luad_cohort$barcode #Identification of samples in metadata after manual curation
colnames(luad_cohort) <- as.character(gsub("\\.", "-", colnames(luad_cohort))) #Sample ID padronization
luad_cohort <- luad_cohort[, colnames(luad_cohort) %in% valid_samples] #Expression table filtered

luad_normal_cohort <- read.csv("FAZER/Expression/tables/rna/others/matrix_luadR_normal.csv", row.names = 1)
md_luad_normal_cohort <- read.csv("FAZER/Expression/tables/metadata/others/metadata_luadR_normal.csv")

  # Normal samples selection
lung_samples <- 367 #This will provide 367 samples 

lung_samples_sel <- sample(md_luad_normal_cohort$barcode, lung_samples)
matrix_lung_sel <- luad_normal_cohort[, colnames(luad_normal_cohort) %in% lung_samples_sel]
md_lung_sel <- md_luad_normal_cohort[md_luad_normal_cohort$barcode %in% colnames(matrix_lung_sel), ]
md_lung_sel$X <- NULL

  # Tumor samples selection
normal_del <- md_luad_cohort$barcode[md_luad_cohort$condition == "Normal"]
luad_cohort <- luad_cohort[, !(colnames(luad_cohort) %in% normal_del)] 

md_luad_cohort <- md_luad_cohort[, c("barcode", "tissue_or_organ_of_origin", "condition", "project_id")] #Should rest only 4 variables
md_luad_cohort <- md_luad_cohort[md_luad_cohort$barcode %in% colnames(luad_cohort), ] 

luad_samples_sel <- sample(md_luad_cohort$barcode, lung_samples)
luad_cohort <- luad_cohort[, colnames(luad_cohort) %in% luad_samples_sel]
md_luad_cohort_sel <- md_luad_cohort[md_luad_cohort$barcode %in% colnames(luad_cohort), ] #Should rest only 367 variables

  # Data combining
matrix_lung_sel <- SummarizedExperiment(
  assays = list(counts = as.matrix(matrix_lung_sel)),
  rowData = DataFrame(geneID = rownames(matrix_lung_sel))
)
matrix_lung_sel <- assay(matrix_lung_sel)
rownames(matrix_lung_sel) <- gsub("\\..*$", "", rownames(matrix_lung_sel))
matrix_lung_sel <- matrix_lung_sel[!duplicated(rownames(matrix_lung_sel)), ] #Remove double lines

genes_normal <- rownames(matrix_lung_sel)
genes_tumor <- rownames(luad_cohort)
genes_comum <- intersect(genes_normal, genes_tumor)

luad_cohort <- luad_cohort[rownames(luad_cohort) %in% genes_comum, ] #Tumor
matrix_lung_sel <- matrix_lung_sel[rownames(matrix_lung_sel) %in% genes_comum, ] #Normal

  # Filtering by autophagy genes 
autophagy_genes <- read.csv("FAZER/ARGs/segunda parte/genes.csv") 
head(autophagy_genes)
ensembl_genes <- autophagy_genes$ensembl_gene_id  #Columns in gene table with ensembl ID

luad_cohort  <- luad_cohort[rownames(luad_cohort) %in% ensembl_genes, ] 
matrix_lung_sel  <- matrix_lung_sel[rownames(matrix_lung_sel) %in% ensembl_genes, ] 

  # Binding tables
matrix_lung_sel <- matrix_lung_sel[rownames(luad_cohort), ] #Reorder the rownames to match
matrix_luadR <- cbind(matrix_lung_sel, luad_cohort)
metadata_luadR <- rbind(md_lung_sel, md_luad_cohort_sel)

  #Saving tables
write.csv(matrix_luadR, "FAZER/Expression/tables/rna/cohort/matrix_luadR_comb.csv", row.names = TRUE) #expression matrix
write.csv(metadata_luadR, "FAZER/Expression/tables/metadata/cohort/metadata_luadR_comb.csv", row.names = TRUE) #metadata

### COAD ###
coad_cohort <- read.csv("FAZER/Expression/tables/rna/others/expression_coad.csv", row.names = 1) #Previous expression table (genes x samples)
md_coad_cohort <- read.csv("FAZER/Expression/tables/metadata/others/metadata_coadR_f.csv") #Metadata table after manual curation (in .csv) (sample x metadata)

valid_samples <- md_coad_cohort$barcode #Identification of samples in metadata after manual curation
colnames(coad_cohort) <- as.character(gsub("\\.", "-", colnames(coad_cohort))) #Sample ID padronization
coad_cohort <- coad_cohort[, colnames(coad_cohort) %in% valid_samples] #Expression table filtered

coad_normal_cohort <- read.csv("FAZER/Expression/tables/rna/others/matrix_coadR_normal.csv", row.names = 1)
md_coad_normal_cohort <- read.csv("FAZER/Expression/tables/metadata/others/metadata_coadR_normal.csv")

  # Normal samples selection
sigmoid_samples <- 184 
transverse_samples <- 183 #This will provide 367 normal samples 

md_coad_sigmoid <- md_coad_normal_cohort[md_coad_normal_cohort$tissue_or_organ_of_origin == "Colon Sigmoid", ]
sigmoid_samples_sel <- sample(md_coad_sigmoid$barcode, sigmoid_samples)
matrix_sigmoid_sel <- coad_normal_cohort[, colnames(coad_normal_cohort) %in% sigmoid_samples_sel]
md_sigmoid_sel <- md_coad_normal_cohort[md_coad_normal_cohort$barcode %in% colnames(matrix_sigmoid_sel), ]

md_coad_transverse <- md_coad_normal_cohort[md_coad_normal_cohort$tissue_or_organ_of_origin == "Colon transverse", ]
transverse_samples_sel <- sample(md_coad_transverse$barcode, transverse_samples)
matrix_transverse_sel <- coad_normal_cohort[, colnames(coad_normal_cohort) %in% transverse_samples_sel]
md_transverse_sel <- md_coad_normal_cohort[md_coad_normal_cohort$barcode %in% colnames(matrix_transverse_sel), ]

matrix_comb <- cbind(matrix_sigmoid_sel, matrix_transverse_sel)
metadata_comb <- rbind(md_sigmoid_sel, md_transverse_sel)
metadata_comb$X <- NULL

  # Tumor samples selection
colon_samples <- 367

normal_del <- md_coad_cohort$barcode[md_coad_cohort$condition == "Normal"]
coad_cohort <- coad_cohort[, !(colnames(coad_cohort) %in% normal_del)] 

md_coad_cohort <- md_coad_cohort[, c("barcode", "tissue_or_organ_of_origin", "condition", "project_id")] #Should rest only 4 variables
md_coad_cohort <- md_coad_cohort[md_coad_cohort$barcode %in% colnames(coad_cohort), ] 

coad_samples_sel <- sample(md_coad_cohort$barcode, colon_samples)
coad_cohort <- coad_cohort[, colnames(coad_cohort) %in% coad_samples_sel]
md_coad_cohort_sel <- md_coad_cohort[md_coad_cohort$barcode %in% colnames(coad_cohort), ] #Should rest only 367 variables

  # Data combining
matrix_comb <- SummarizedExperiment(
  assays = list(counts = as.matrix(matrix_comb)),
  rowData = DataFrame(geneID = rownames(matrix_comb))
)
matrix_comb <- assay(matrix_comb)
rownames(matrix_comb) <- gsub("\\..*$", "", rownames(matrix_comb))
matrix_comb <- matrix_comb[!duplicated(rownames(matrix_comb)), ] #Remove double lines

genes_normal <- rownames(matrix_comb)
genes_tumor <- rownames(coad_cohort)
genes_comum <- intersect(genes_normal, genes_tumor)

coad_cohort <- coad_cohort[rownames(coad_cohort) %in% genes_comum, ] #Tumor
matrix_comb <- matrix_comb[rownames(matrix_comb) %in% genes_comum, ] #Normal

  # Filtering by autophagy genes 
autophagy_genes <- read.csv("FAZER/ARGs/segunda parte/genes.csv") 
head(autophagy_genes)
ensembl_genes <- autophagy_genes$ensembl_gene_id  #Columns in gene table with ensembl ID

coad_cohort  <- coad_cohort[rownames(coad_cohort) %in% ensembl_genes, ] 
matrix_comb  <- matrix_comb[rownames(matrix_comb) %in% ensembl_genes, ] 

  # Binding tables
matrix_comb <- matrix_comb[rownames(coad_cohort), ] #Reorder the rawnames to match
matrix_coadR <- cbind(matrix_comb, coad_cohort)
metadata_coadR <- rbind(metadata_comb, md_coad_cohort_sel)

  #Saving tables
write.csv(matrix_coadR, "FAZER/Expression/tables/rna/cohort/matrix_coadR_comb.csv", row.names = TRUE) #expression matrix
write.csv(metadata_coadR, "FAZER/Expression/tables/metadata/cohort/metadata_coadR_comb.csv", row.names = TRUE) #metadata

