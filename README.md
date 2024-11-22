# ukbbPheSurv
This repositary is used for generating of survival phenotype.

## Required R packages
- data.table
- tidyr
- parallel
- intervals
- XML
- RCurl
- RCurl
- htmltab
- bitops
- stringr
- stringr
- zoo


## Step 0: ./script/function.collectFields.r   
Collect fields data of UKBB   
output: ./data/FieldListing.txt    

## Step 1: ./Step1_createPheCode.R   
Create phenome and phecode files of UKBB    
**input**:     
UK Biobank main data (*.csv)    
UK Biobank genome fam data (*.fam)    
**output**:
./results/UKB_PHENOME_DESCRIPTION.txt    
./results/UKB_PHENOME.txt    
./results/UKB_PHENOME_ICD9_PHECODE_MAP.txt    
./results/UKB_PHENOME_ICD10_PHECODE_MAP.txt    
./results/UKB_PHENOME_UNMAPPED_ICD9_CODES.txt    
./results/UKB_PHENOME_UNMAPPED_ICD10_CODES.txt    

## Step 2: ./Step2_mapPheCode.R    
Map ICD9 and ICD10    
**input**:    
./results/UKB_PHENOME_UNMAPPED_ICD9_CODES.txt    
./results/UKB_PHENOME_UNMAPPED_ICD10_CODES.txt    
**output**:    
/results/UKB_PHECODE_v1.2b1_ICD_Mapping.txt    

## Step 3: ./Step3_survPheCode.R    
Create survival phenotype    
**input**:    
UK Biobank main data (*.csv)    
UK Biobank genome fam data (*.fam)    
./data/FieldListing.txt    
./results/UKB_PHECODE_v1.2b1_ICD_Mapping.txt    
./results/UKB_PHENOME_DESCRIPTION.txt    
./results/UKB_PHENOME.txt    
**output**:
./results/*_pheSurv.txt    
