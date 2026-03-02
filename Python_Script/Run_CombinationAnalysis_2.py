import pandas as pd
import ast
import itertools
import numpy as np


phenotype_common=pd.read_excel("GATA2_Phenotypes_Final_v8.xlsx")

phenotype_common = phenotype_common[(phenotype_common['Present_in_GATA2'] == 'Yes') & (phenotype_common['Present_in_UKB'] == 'Yes') & (phenotype_common['Keep_Adults']=='Yes')]
phenotype_common

pheno_participant=phenotype_common[["Phenotype","Phenotype_Description","Phenotype_Index"]]

###----------------- Generate Combinations of 2 for the SerialNumber column
combinations = list(itertools.combinations(phenotype_common['Phenotype_Index'].astype(str), 2))



# Create a new DataFrame with combinations
combinations_df = pd.DataFrame(['_'.join(combination) for combination in combinations], columns=['Combination'])
# Create a new DataFrame with combinations
combinations_df = pd.DataFrame({
    'Combination': ['_'.join(combination) for combination in combinations],
    'PhenotypeSet': [combination for combination in combinations]  # Store original combinations
})
combinations_df





#PID_Dataframe
pid_dataframe=pd.read_excel('Phenotype_pid_final_022526_v8.xlsx')






### Adding columns with number of each phenotypes
combinations_df["Counts_UKB"]=None #IMPORTANT





for idx,row in combinations_df.iterrows():
  print(f"Combination Number {idx+1}")
  print(f"Combination {combinations_df['Combination'][idx]}")
  combo=combinations_df.loc[idx,"Combination"] 
  pid_pheno1_index=combo.split("_")[0]
  pid_pheno1=pid_dataframe[pid_pheno1_index][pid_dataframe[pid_pheno1_index].astype(bool)]
  pid_pheno2_index=combo.split("_")[1]
  pid_pheno2=pid_dataframe[pid_pheno2_index].astype(bool)
  pid_pheno2=pid_dataframe[pid_pheno2_index][pid_dataframe[pid_pheno2_index].astype(bool)]
  intersection = set(pid_pheno1).intersection(set(pid_pheno2))
  len(intersection)
  #combinations_df.loc[idx,"Combined_Unique_PIDs"]=[[list(intersection)]]
  combinations_df.loc[idx,"Counts_UKB"]=len(intersection)
  
  #combinations_df.loc[idx,"Combined_Unique_PIDs"]=[list(set(pid_pheno1).intersection(pid_pheno2))]
  #combinations_df.loc[idx,"Final_Frequency_of_Combination"]=len(list(set(pid_pheno1).intersection(pid_pheno2)))

 
combinations_df.to_excel('Combinationsof2_FrequencyOnly_022526.xlsx', index=False)
