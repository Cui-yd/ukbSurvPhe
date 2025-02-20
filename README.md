# ukbSurvPhe
This repositary is used for generating the time-to-event phenotypes in UK Biobank based on PheCodes (using ICD9 and ICD10).

- The **failure time** was defined as the interval, rounded to the nearest full month, from their birth year and month to the diagnosis date of the relevant PheCode-specific PheCodes.

- Subjects with relevant PheCodes but without corresponding diagnosis dates were excluded from the analysis.

- Subjects without relevant PheCodes and without exclusion PheCodes as defined by the PheCode were classified as **censored**. For these subjects, the censoring time was defined as the interval, rounded to the nearest full month, from their birth year and month to the latest of the following: the most recent recorded diagnosis date of any PheCodes, the lost-to-follow-up date, or the recorded date of death.



## Required R packages
<a href="https://cran.r-project.org/web/packages/data.table/index.html">data.table</a>,
<a href="https://cran.r-project.org/web/packages/tidyr/index.html">tidyr</a>, 
<a href="https://cran.r-project.org/web/packages/parallelly/index.html">parallel</a>, 
<a href="https://cran.r-project.org/web/packages/intervals/index.html">intervals</a>, 
<a href="https://cran.r-project.org/web/packages/XML/index.html">XML</a>, 
<a href="https://cran.r-project.org/web/packages/RCurl/index.html">RCurl</a>, 
<a href="https://cran.r-project.org/src/contrib/Archive/htmltab/">htmltab</a>, 
<a href="https://cran.r-project.org/web/packages/bitops/index.html">bitops</a>, 
<a href="https://cran.r-project.org/web/packages/stringr/index.html">stringr</a>, and 
<a href="https://cran.r-project.org/web/packages/zoo/index.html">zoo</a>, 



## Usage

### Step 0: Run -> ./script/function.collectFields.r
Collect fields data of UKBB   
output: ./data/FieldListing.txt    

### Step 1: Run -> ./Step1_createPheCode.R    
Create phenome and phecode files of UKBB    
**input**:     
UK Biobank main data (\*.csv)    
UK Biobank genome fam data (\*.fam)    
**output**:    
./results/UKB_PHENOME_DESCRIPTION.txt    
./results/UKB_PHENOME.txt    
./results/UKB_PHENOME_ICD9_PHECODE_MAP.txt    
./results/UKB_PHENOME_ICD10_PHECODE_MAP.txt    
./results/UKB_PHENOME_UNMAPPED_ICD9_CODES.txt    
./results/UKB_PHENOME_UNMAPPED_ICD10_CODES.txt    

### Step 2: Run -> ./Step2_mapPheCode.R    
Map ICD9 and ICD10    
**input**:    
./results/UKB_PHENOME_UNMAPPED_ICD9_CODES.txt    
./results/UKB_PHENOME_UNMAPPED_ICD10_CODES.txt    
**output**:    
/results/UKB_PHECODE_v1.2b1_ICD_Mapping.txt    

### Step 3: Run -> ./Step3_survPheCode.R    
Create time-to-event phenotypes    
**input**:    
UK Biobank main data (\*.csv)    
UK Biobank genome fam data (\*.fam)    
./data/FieldListing.txt    
./results/UKB_PHECODE_v1.2b1_ICD_Mapping.txt    
./results/UKB_PHENOME_DESCRIPTION.txt    
./results/UKB_PHENOME.txt    
**output**:    
./phenotypes/*_pheSurv_useAttendTime.txt    



## Reference
https://github.com/umich-cphds/createUKBphenome
https://github.com/atgu/ukbb_pan_ancestry/blob/master/create_phecode_mapping.R



## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
