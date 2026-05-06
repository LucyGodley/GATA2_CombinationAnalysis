#load libraries

library("readxl")
#

phenotype_eid_df<-read.csv("phenotype_eid_df_Adults_v8.csv")
combinations_df_3<-read_xlsx("Combinationsof3_FrequencyOnly_022526.xlsx")
combinations_df_3<-combinations_df_3[,-2]
colnames(combinations_df_3)[2]<-"Counts_UKB"
combinations_df_3$Counts_GATA2<-rep(NA,nrow(combinations_df_3))
colnames(phenotype_eid_df) <- trimws(colnames(phenotype_eid_df))


for (j in 1:nrow(combinations_df_3)) {
  print(j)
  
  for (k in 1:3) {
    pheno_idx <- trimws(strsplit(combinations_df_3$Combination[j], "_")[[1]][k])
    print(pheno_idx)
    if (k == 1) {
      pheno_1 <- list(phenotype_eid_df[!is.na(phenotype_eid_df[,pheno_idx]),pheno_idx])
    }
    else if(k==2){
      pheno_2<-list(phenotype_eid_df[!is.na(phenotype_eid_df[,pheno_idx]),pheno_idx])
      #print(pheno_1)
    }
    else {
      pheno_3 <- list(phenotype_eid_df[!is.na(phenotype_eid_df[,pheno_idx]),pheno_idx])
    }
  }
  
  common <- Reduce(intersect, list(unlist(pheno_1), unlist(pheno_2),unlist(pheno_3)))
  combinations_df_3$Counts_GATA2[j] <- length(common)
  print(paste0("Number of participants=", length(common)))
}

write.csv(combinations_df_3,"Combinationsof3_GATA2_022526.csv",row.names = FALSE)


