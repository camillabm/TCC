#######################
### ARGs definition ###
#######################

### Explanation
# Here, we will extract ensembl gene and transcript annotations for the genes previous selected.

###### Libraries
install.packages("readxl")
library(readxl)

####### ARGs table as csv
df <- read_excel("FAZER/ARGs/segunda parte/genes.xlsx")
write.csv(df, "FAZER/ARGs/segunda parte/genes.csv", row.names = FALSE)
ARG <- read.csv("FAZER/ARGs/segunda parte/genes.csv")

####### Ensemble ID annotation
install.packages("biomaRt")
library (biomaRt)
library (dplyr)

  # ENSG (gene)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://grch37.ensembl.org") #Connection to ensembl

hgnc_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"), 
  filters = "hgnc_symbol",  # Agora filtramos pelo HGNC
  values = ARG$HGNC,  # Passamos os HGNC codes
  mart = ensembl
)

ARG <- merge(ARG, hgnc_mapping, by.x = "HGNC", by.y = "hgnc_symbol", all.x = TRUE) #Add HGNC symbol to results

ARG <- ARG %>%
  mutate(ensembl_gene_id = if_else(HGNC == "ATG101", "ENSG00000123395", ensembl_gene_id))
ARG <- ARG %>%
  mutate(ensembl_gene_id = if_else(HGNC == "GBA1", "ENSG00000177628", ensembl_gene_id))
ARG <- ARG %>%
  mutate(ensembl_gene_id = if_else(HGNC == "PIP4P1", "ENSG00000165782", ensembl_gene_id))
ARG <- ARG %>%
  mutate(ensembl_gene_id = if_else(HGNC == "RAB7B", "ENSG00000276600", ensembl_gene_id))

write.csv(ARG, "FAZER/ARGs/segunda parte/genes.csv", row.names = FALSE)

  # ENST (transcript)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://grch37.ensembl.org") #Connection to ensembl

result <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
                filters = "ensembl_gene_id",
                values = ARG$ensembl_gene_id,
                mart = ensembl)

ARG <- merge(ARG, result, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE) #Add HGNC symbol to results


