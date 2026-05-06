# %%
import pandas as pd

# %%
gata2_phenotypes = pd.read_excel(r"GATA2_Phenotypes_Final_v8.xlsx")
gata2_phenotypes

# %%
# Step 2: Load the ICD10 code descriptions file
icd10_df = pd.read_csv(r"coding19.tsv", sep='\t')
icd10_df


# Assume: 'coding' is the column with code, 'meaning' is the description
icd10_dict = dict(zip(icd10_df['coding'], icd10_df['meaning']))

# %%
gata2_phenotypes['ICD10_code'] = gata2_phenotypes['ICD10_code'].astype(str)

# Remove brackets and quotes, then split
gata2_phenotypes['ICD10_List'] = (
    gata2_phenotypes['ICD10_code']
    .str.replace(r"[\[\]']", "", regex=True)
    .str.split(',')
)

# Strip whitespace
gata2_phenotypes['ICD10_List'] = gata2_phenotypes['ICD10_List'].apply(lambda x: [code.strip() for code in x if code.strip() != ''])

# Explode to separate rows
exploded_df = gata2_phenotypes.explode('ICD10_List')
exploded_df['ICD10_List'] = exploded_df['ICD10_List'].str.strip()

# %%
# Map ICD10 code to meaning
exploded_df['ICD10_Description'] = exploded_df['ICD10_List'].map(icd10_dict)
exploded_df

# %%
exploded_df.to_csv(r"GATA2_Phenotypes_ICD10_Descriptions_Final_v8.csv", index=False)

# %%
icd10_df=exploded_df.copy()
icd10_df

# %%
#--------------------------Loading Dataset
#Load the ICD10 Diagnosis data 
data=pd.read_csv("DiagnosisICD10_41270.csv",sep =",",low_memory=False)
print(data.shape) 
data=data.set_index('eid')

# %%
#Remove all the withdrawn participants from the dataset
withdrawn_pid=pd.read_csv(r"withdraw83200_480_20260109.txt",sep="\t",header=None)

# %%
data=data[~data.index.isin(withdrawn_pid[0])]
print("UKBiobank Control Dataset after removing withdrawn PID=",len(data.axes[0]))
#-----------------------------Clean Data
#Set PID as row index

data_df=data.dropna(inplace=False,how='all')
# ^ this line removes all the blank rows in UKBiobank database

print(len(data.axes[0])) # This command tells you the number of patients in the entire dataset
print(len(data_df.axes[0]))# Number of rows which has at least one value in one of the columns.
#N=440017/New number 446837


no_diagnosis=list(set(data.index)-set(data_df.index))
print("Number of rows/Participants with no ICD10 codes",len(no_diagnosis))


# %% [markdown]
# Remove participants who withdrew

# %%
variant_pid_list=pd.read_excel(r"UKB_Variants_PID_227Genes_Final.xlsx")

print("Number of Participants with P/LP variant",variant_pid_list['ParticipantID'].nunique(dropna=True))

# %%
#Remove all rows with blank ParticipantID
variant_pid_list = variant_pid_list.dropna(subset=["ParticipantID"])

# %%
#Read all the GATA2 participants in UKBiobank
gata2_variants_pid=pd.read_csv(r"GATA2_Variants_UKB.csv")

# %%
gata2_variants_pid['Source_File'].value_counts()

# %%
variant_pid = list(
    set(
        variant_pid_list['ParticipantID'].dropna().tolist()
        + gata2_variants_pid['column'].dropna().tolist()
    )
)
print("Number of Participants with P/LP variant in 227 genes including GATA2",len(variant_pid))

# %% [markdown]
# Create a control for Combination analysis
# Remove all the participants with P/LP variants in 227 genes

# %%
control_df=data_df.loc[~data_df.index.isin(variant_pid)]
print("Number of Control Participants without P/LP variant",len(control_df.axes[0]))
control_df

# %% [markdown]
# Subset GATA2 Cases

# %%
# get ParticipantIDs for rows where Gene == 'GATA2'
gata2_pids = variant_pid_list.loc[variant_pid_list['Gene'] == 'GATA2', 'ParticipantID'].dropna()

# ensure IDs match the dtype of data.index (convert to int where possible)
try:
	gata2_pids = gata2_pids.astype(int).tolist()
except Exception:
	gata2_pids = [int(float(x)) for x in gata2_pids]

#gata2_pids

# 3. Combine, normalize dtype, and keep unique (order preserved)
# Combine with gata2_variants_pid['column'], drop NA, normalize dtype, keep unique (order preserved)
gata2_pids = (
    pd.Series(gata2_pids + gata2_variants_pid['column'].dropna().tolist())
      .astype(float)
      .astype(int)
      .unique()
      .tolist()
)

print("Number of GATA2 PLP Carriers:",len(gata2_pids))

# subset data by those participant IDs (eids)
gata2_cases = data_df.loc[data_df.index.isin(gata2_pids)]
print("Number of GATA2 cases with non missing ICD10 dataset:", len(gata2_cases.axes[0]))
gata2_cases.to_csv("GATA2_177Cases_ICD10.csv")

# %% [markdown]
# Remove participants who have any cancer

# %%
#Read all the participants who do not have any the cancer ICD10 code in p40006,p41270,p41271
noncancer_ukb=pd.read_csv(r"Negative_Cancer_Phenotype_ICD10_ICD9_Cancer.csv")
noncancer_ukb

# %%
control_df=control_df[control_df.index.isin(noncancer_ukb['eid'])]
print("UKBiobank Control Dataset=",len(control_df.axes[0]))

# %%
# Select columns correctly using a list
data_df_merge = control_df.merge(
	noncancer_ukb[['eid', 'Sex', 'Age_at_recruitment']],
	how='left',
	on='eid'
)

# %%
print(data_df_merge['Sex'].value_counts())



# %%
# Select columns correctly using a list
gata2_df_merge = gata2_cases.merge(
	noncancer_ukb[['eid', 'Sex', 'Age_at_recruitment']],
	how='left',
	on='eid'
)

len(gata2_cases)


# %%
print(gata2_df_merge['Sex'].value_counts())

# %%
import pandas as pd

# Assume data_df looks like this:
# Columns: eid, diag_1, diag_2, diag_3, ...
# Example ICD10 codes are present in diag_* columns

# STEP 1: Define columns to search (all columns except 'eid')
icd10_cols = [col for col in control_df.columns if col != 'eid']

# STEP 2: Flatten the unique ICD10 codes from phenotype_df
all_icd10_codes = set(icd10_df['ICD10_List'].dropna().unique())

# STEP 3: For each ICD10 code, find eids
icd10_to_eids = {}

for icd10 in all_icd10_codes:
    # Check for presence in any diag_* column
    mask = control_df[icd10_cols].isin([icd10]).any(axis=1)
    eids = control_df.index[mask].tolist()
    icd10_to_eids[icd10] = eids

# OPTIONAL: save as dictionary or dataframe
# Example: convert to DataFrame
icd10_eid_df = pd.DataFrame([
    {'ICD10': code, 'EIDs': eids, 'N': len(eids)}
    for code, eids in icd10_to_eids.items()
])





# %%
# Find max list length for padding
max_len = max(len(eids) for eids in icd10_to_eids.values())

# Pad lists with NaN and create DataFrame
icd10_eid_df = pd.DataFrame({
    code: eids + [pd.NA] * (max_len - len(eids))
    for code, eids in icd10_to_eids.items()
})




# %%
# 'coding' is the column with code, 'meaning' is the description
icd10_dict = dict(zip(icd10_df['coding'], icd10_df['meaning']))


# %%
icd10_eid_df_labels=icd10_eid_df.copy()
# After icd10_eid_df has been created (columns are ICD10 codes), map meanings into column names
def make_col_label(code):
    meaning = icd10_dict.get(code)
    if meaning is None or (isinstance(meaning, float) and pd.isna(meaning)):
        return str(code)
    # short-clean the meaning to avoid very long names if desired; here keep full meaning
    return f"{meaning}"

# apply mapping
new_columns = [make_col_label(col) for col in icd10_eid_df.columns]

icd10_eid_df_labels.columns = new_columns
icd10_eid_df_labels
# Optional: if you need to also keep the original codes, create a mapping DataFrame
# colmap_df = pd.DataFrame({'code': icd10_eid_df.columns, 'label': new_columns})

# Save to CSV/Excel with new descriptive headers
icd10_eid_df.to_csv("ICD10_EID_Match_Columns_with_meaning_030226.csv", index=False)


# %%
import pandas as pd

phenotype_df = pd.read_excel("GATA2_Phenotypes_Final_v8.xlsx")
phenotype_df

# %%

# STEP 1: Create a dictionary to hold lists of eid for each Phenotype_Index
phenotype_index_to_eids = {}

for idx, row in phenotype_df.iterrows():
    phenotype_index = row['Phenotype_Index']
    
    # Clean ICD10 code list for this phenotype
    icd10_codes = row['ICD10_code']
    if isinstance(icd10_codes, str):
        icd10_codes = icd10_codes.replace("[", "").replace("]", "").replace("'", "").split(",")
        icd10_codes = [code.strip() for code in icd10_codes if code.strip() != '']
    else:
        icd10_codes = []
    
    # Collect all eids who have at least one of these ICD10 codes
    eid_set = set()
    for code in icd10_codes:
        if code in icd10_eid_df.columns:
            eids = icd10_eid_df[code].dropna().tolist()
            eid_set.update(eids)
    
    phenotype_index_to_eids[phenotype_index] = list(eid_set)

# STEP 2: Determine max length for padding
max_len = max(len(eids) for eids in phenotype_index_to_eids.values())

# Create DataFrame with columns = Phenotype_Index, values = eid lists (padded with NA)
phenotype_eid_wide_df = pd.DataFrame({
    phenotype_index: eids + [pd.NA] * (max_len - len(eids))
    for phenotype_index, eids in phenotype_index_to_eids.items()
})






# %%
import pandas as pd

# Read the EIDs from the text file
with open("participant_id_low_monocytes.txt", "r") as file:
    low_mono_eids = [line.strip() for line in file if line.strip() != '']

# Ensure they are unique (optional)
low_mono_eids = list(set(low_mono_eids))

# Replace/add them in P37 column:
max_len = max(len(low_mono_eids), len(phenotype_eid_wide_df))

# Ensure phenotype_eid_wide_df can accommodate the new length
if len(phenotype_eid_wide_df) < max_len:
    # Pad existing DataFrame with rows
    additional_rows = max_len - len(phenotype_eid_wide_df)
    padding_df = pd.DataFrame([ [pd.NA] * len(phenotype_eid_wide_df.columns) ] * additional_rows, columns=phenotype_eid_wide_df.columns)
    phenotype_eid_wide_df = pd.concat([phenotype_eid_wide_df, padding_df], ignore_index=True)

# Now replace the P37 column
phenotype_eid_wide_df['P37'] = low_mono_eids + [pd.NA] * (max_len - len(low_mono_eids))

# View result
print(phenotype_eid_wide_df)

# OPTIONAL: Save updated DataFrame
phenotype_eid_wide_df.to_excel("Phenotype_EID_Wide_v8_030226.xlsx", index=False)


# %%
# 1. Remove withdrawn participants
withdrawn_set = set(withdrawn_pid[0].astype(str).str.strip().values)
low_mono_eids_trim = [
    eid.strip() for eid in low_mono_eids 
    if str(eid).strip() not in withdrawn_set
]

# 2. Remove cancer-risk eids (pid_cancer_risk_all)
if 'eid' in variant_pid_list.columns:
    cancer_risk_set = set(variant_pid_list['ParticipantID'].astype(str).str.strip().values)
else:
    cancer_risk_set = set(variant_pid_list.iloc[:, 0].astype(str).str.strip().values)

low_mono_eids_trim = [
    eid for eid in low_mono_eids_trim 
    if str(eid).strip() not in cancer_risk_set
]

# 3. Remove eids present in no_diagnosis
no_diag_set = {str(x).strip() for x in no_diagnosis}

low_mono_eids_trim = [
    eid for eid in low_mono_eids_trim 
    if str(eid).strip() not in no_diag_set
]

print("Number of low monocyte eids after trimming:", len(low_mono_eids_trim))

# %%
import pandas as pd

gata2_pid_clean = (
    pd.Series(gata2_pids)
      .astype(str)
      .str.strip()
      .pipe(pd.to_numeric, errors="coerce")
      .dropna()
      .astype(int)
      .tolist()
)

low_mono_eids_trim_clean = (
    pd.Series(low_mono_eids_trim)
      .astype(str)
      .str.strip()
      .pipe(pd.to_numeric, errors="coerce")
      .dropna()
      .astype(int)
      .tolist()
)

# overlap check
overlap = list(set(gata2_pid_clean) & set(low_mono_eids_trim_clean))



# %%
low_mono_eids_trim_clean = list(
    set(low_mono_eids_trim_clean) - set(gata2_pid_clean)
)


# %%
# ...existing code...
# Ensure consistent string types and trim whitespace
eid_set = set(str(x).strip() for x in low_mono_eids_trim_clean)

selected = [['eid', 'Sex', 'Age_at_recruitment']].copy()
selected = selected[selected['eid'].astype(str).str.strip().isin(eid_set)]

# View or save
print(selected.shape)

# ...existing code...

# %%
# ...existing code...
# Create data_merge_final from data_df_merge then add selected and remove duplicate eids
data_merge_final = data_df_merge[['eid', 'Sex', 'Age_at_recruitment']].copy()

# Normalize eid strings to avoid mismatches
data_merge_final['eid'] = data_merge_final['eid'].astype(str).str.strip()
selected['eid'] = selected['eid'].astype(str).str.strip()

# Append selected and drop duplicated eids (keep first occurrence)
data_merge_final = pd.concat([data_merge_final, selected], ignore_index=True)
data_merge_final = data_merge_final.drop_duplicates(subset='eid', keep='first').reset_index(drop=True)

# Inspect / save
print(data_merge_final.shape)
data_merge_final.head()
# Optionally save:
data_merge_final.to_csv("Control_GATA2_EID_Final_030226.csv", index=False)
# ...existing code...

# %%
# Ensure all columns are padded to max rows (optional but safe step)
max_len = max(phenotype_eid_wide_df[col].notna().sum() for col in phenotype_eid_wide_df.columns)

if len(phenotype_eid_wide_df) < max_len:
    additional_rows = max_len - len(phenotype_eid_wide_df)
    padding_df = pd.DataFrame(
        [[pd.NA] * len(phenotype_eid_wide_df.columns)] * additional_rows,
        columns=phenotype_eid_wide_df.columns
    )
    phenotype_eid_wide_df = pd.concat([phenotype_eid_wide_df, padding_df], ignore_index=True)

# Fill blanks (NaN) with 0 (as integer or string depending on your need)
phenotype_eid_wide_df = phenotype_eid_wide_df.fillna(0)

# OPTIONAL: Convert to integer type if desired (only works if all values are numeric after fillna)
# phenotype_eid_wide_df = phenotype_eid_wide_df.astype(int)

# Save to Excel
phenotype_eid_wide_df.to_excel("Phenotype_pid_final_022526_v8.xlsx", index=False)



