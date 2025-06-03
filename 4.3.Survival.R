#####################
### LGG Survival ###
#####################

### Explanation
# Autophagy states were defined by state plot of AutoIndex and LysoIndex. 
# The inhibition samples were selected to match the number of induction samples.

###### Libraries
library(readxl)

###### Paths and data
states_g2 <- read_excel("FAZER/GBM/Survival/g2_all.xlsx")
write.csv(states_g2, "FAZER/GBM/Survival/g2_all.csv")
states_g2 <- read.csv("FAZER/GBM/Survival/g2_all.csv")

states_g3 <- read_excel("FAZER/GBM/Survival/g3_all.xlsx")
write.csv(states_g3, "FAZER/GBM/Survival/g3_all.csv")
states_g3 <- read.csv("FAZER/GBM/Survival/g3_all.csv")

metadata <- read.csv("FAZER/GBM/LGG/lgg_md.csv")
clinical <- read.csv("FAZER/GBM/LGG/metadata_lggR.csv")

###### Clinical data biding
clinical_subset <- clinical %>% dplyr::select(barcode, paper_Survival..months.)
merged <- left_join(metadata, clinical_subset, by = "barcode")
merged$paper_Survival..months. <- round(merged$paper_Survival..months. * 30)
merged$Status <- "1"

####### States
states <- rbind (states_g2, states_g3)
md_final <- left_join(merged, states, by = "barcode")

write.csv(md_final, "FAZER/GBM/Survival/grades_ready.csv")

