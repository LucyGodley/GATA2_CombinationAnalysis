#load libraries

library("readxl")
#

combinations_df_2<-read_xlsx("Combinationsof2_FrequencyOnly_022526.xlsx")
combinations_df_2<-combinations_df_2[,-c(2)]
colnames(combinations_df_2)[2]<-"Counts_UKB"
combinations_df_2$Counts_GATA2<-rep(NA,nrow(combinations_df_2))
colnames(phenotype_eid_df) <- trimws(colnames(phenotype_eid_df))


for (j in 1:nrow(combinations_df_2)) {
  print(j)
  
  for (k in 1:2) {
    pheno_idx <- trimws(strsplit(combinations_df_2$Combination[j], "_")[[1]][k])
    print(pheno_idx)
    if (k == 1) {
      pheno_1 <- list(phenotype_eid_df[!is.na(phenotype_eid_df[,pheno_idx]),pheno_idx])
    } else {
      pheno_2 <- list(phenotype_eid_df[!is.na(phenotype_eid_df[,pheno_idx]),pheno_idx])
    }
  }
  
  common <- Reduce(intersect, list(unlist(pheno_1), unlist(pheno_2)))
  combinations_df_2$Counts_GATA2[j] <- length(common)
  print(paste0("Number of participants=", length(common)))
}

write.csv(combinations_df_2,"Combinationsof2_GATA2_022526.csv",row.names = FALSE)


