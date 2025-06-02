library(dplyr)
library(stringr)

peri_cnv <- read.csv("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/CIN_Features/CNV_Data/TCGA_pericentro_cnv.csv")

print(head(peri_cnv[,1:10]))


peri_cnv$sampleID <- peri_cnv$sampleID %>% str_replace_all(.,"-",".")

write.csv( peri_cnv, "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/CIN_Features/CNV_Data/TCGA_pericentro_cnv_hpc.csv", row.names = FALSE)