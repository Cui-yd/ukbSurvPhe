##########################################################
# Generate time-to-event phenotypes in UK Biobank data
# Yidan Cui
# Initiate date: 2024/05/17
# Current date: 2025/02/20
##########################################################

options(stringsAsFactors=F)
suppressPackageStartupMessages(library(data.table, quietly = T))
suppressPackageStartupMessages(library(parallel, quietly = T))
suppressPackageStartupMessages(library(tidyverse, quietly = T))
suppressPackageStartupMessages(library(stringr, quietly = T))
suppressPackageStartupMessages(library(zoo, quietly = T))


setDTthreads(detectCores()/2)


### change this path by your own UKBB data
UKB_mainDat = "PATH/TO/*.csv"


field_list = fread("./data/FieldListing.txt")
uk_column = colnames(fread(UKB_mainDat, header = T, sep=",", nrows = 1))[-1]
uk_fieldID = str_sub(uk_column, 1, str_locate(uk_column, "-")-1) %>% unique() %>% as.numeric()
uk_fieldID = str_extract(uk_column, "\\d+") %>% unique() %>% as.numeric()



## phecode mapped icd code
phecode_map = fread("./results/UKB_PHECODE_v1.2b1_ICD_Mapping.txt", data.table = F)
## phecode dictionary
phecode_dic = fread("./results/UKB_PHENOME_DESCRIPTION.txt", data.table = F)
## individuals' phenotype
phenome = fread("./results/UKB_PHENOME.txt", data.table = F)
rownames(phenome) = phenome$IID
phenome = phenome[,-1]


uk_fieldlist = field_list[na.omit(match(uk_fieldID, field_list$`Field ID`)),]
field_ICD9 = uk_fieldlist[grepl("ICD9", Description), c(1:2)]
field_ICD10 = uk_fieldlist[grepl("ICD10", Description), c(1:2)]
field_birth = uk_fieldlist[grepl("birth", Description), c(1:2)][1:2,]
field_date = uk_fieldlist[grepl("Date", Description), c(1:2)][c(2,5,6),]
field_use = rbind(field_ICD9, field_ICD10, field_birth, field_date)
field_use = field_use[order(field_use$`Field ID`), ]
rm(field_birth, field_date, field_ICD10, field_ICD9)


## read ukbb main data and match with genetic data
## 34-0.0 Year of birth;  52-0.0 Month of birth;  53-0.0 Attendance;  191-0.0 Date lost to follow-up
uk_extract = grep(paste(field_use$`Field ID`[-3:-1], sep = "", collapse = "|"), uk_column)
uk_extractName = c("eid", "34-0.0", "52-0.0", "53-0.0", "191-0.0", uk_column[uk_extract])
uk_extractData = fread(UKB_mainDat, header = T, select = uk_extractName, data.table = F)
uk_extractData = uk_extractData[na.omit(match(rownames(phenome), uk_extractData$eid)),]



## match phenome data and phecode data
phenome_colname = gsub("X", "", colnames(phenome)) %>% as.numeric()
summary(phenome_colname == phecode_dic$phecode)
length(which(is.na(match(phenome_colname, phecode_dic$phecode))))
phecode_dic = phecode_dic[na.omit(match(phenome_colname, phecode_dic$phecode)),]
phecode_map = phecode_map[na.omit(match(phenome_colname, phecode_map$phecode)),]
summary(phenome_colname == phecode_dic$phecode)


uk_use = data.frame("eid" = uk_extractData$eid)
uk_use$birth_date = as.Date(NA)
for (j in 1:nrow(uk_use)) uk_use$birth_date[j] = as.Date(paste(uk_extractData$`34-0.0`[j],uk_extractData$`52-0.0`[j],"15", sep = "-"))
uk_use$attend_date = uk_extractData$`53-0.0`
uk_use$death_date = uk_extractData$`40000-0.0`
uk_use$lost_date = uk_extractData$`191-0.0`

### change this date by your own UKBB application date
uk_use$app_date = as.Date("YYYY-MM-DD")



uk_module = list()
for (p in 1:16) {
  extract_p = grep(field_use$`Field ID`[-4:-1][p], colnames(uk_extractData))
  uk_module[[p]] = uk_extractData[, extract_p]
  print(dim(uk_module[[p]]))
}
names(uk_module) = c("primary_death_icd10", "secondry_death_icd10",
                     "cancer_icd10_date", "cancer_icd10",
                     "cancer_icd9", "external_cause_icd10",
                     "diagnose_main_icd10", "diagnose_main_icd9",
                     "diagnose_2nd_icd10", "diagnose_2nd_icd9",
                     "diagnose_main_icd10_date", "diagnose_main_icd9_date",
                     "diagnose_icd10", "diagnose_icd9",
                     "diagnose_icd10_date", "diagnose_icd9_date")


uk_module_icd = cbind(uk_module[[4]], uk_module[[7]], uk_module[[8]], uk_module[[13]], uk_module[[14]])
uk_module_date = cbind(uk_module[[3]], uk_module[[11]], uk_module[[12]], uk_module[[15]], uk_module[[16]])





## input your interested phenotype
phecode_id = "X0000"

phecode_index = which(colnames(phenome) == phecode_id)
phenome_event_index = which(phenome[, phecode_index] == 1)  # extract case
phenome_censor_index = which(phenome[, phecode_index] == 0)  # extract case

## earliest diagnosed date of event individuals
event_date = rep(NA, nrow(phenome))
recent_diag_date = rep(NA, nrow(phenome))
## output survival phenotype
output_path = "./phenotypes/"


if (length(phenome_event_index) == 0) {

  stop("\n", paste0("Phecode[", phecode_index, "] ", colnames(phenome)[phecode_index], ": number of diagnosed individuals = ", length(phenome_event_index), "   SKIP!"))

} else {

  message("\n", paste0("Phecode[", phecode_index, "] ", colnames(phenome)[phecode_index], ": number of diagnosed individuals = ", length(phenome_event_index)))
  icd10_i = strsplit(phecode_map$icd10_code[phecode_index], ",")[[1]]    ## phecode icd10
  icd9_i = strsplit(phecode_map$icd9_code[phecode_index], ",")[[1]]      ## phecode icd9


  ### for event samples
  for (index_j in phenome_event_index) {                                   ## for each diagnosed individual

    sample_icd10_j = which(uk_module_icd[index_j,] %in% icd10_i)
    sample_icd9_j = which(uk_module_icd[index_j,] %in% icd9_i)

    sample_index_j = c(sample_icd10_j, sample_icd9_j)
    sample_date_j = as.numeric(uk_module_date[index_j, sample_index_j])

    if(length(sample_index_j) >= 1) event_date[index_j] = min(sample_date_j)
  }

  event_date = as.Date(event_date, origin = "1970-01-01")
  status_ind = rep(NA, length(event_date))
  status_ind[which(!is.na(event_date))] = 1


  ### for censored samples
  pb <- txtProgressBar(min = 1, max = length(phenome_censor_index), style = 3)
  for (j in 1:length(phenome_censor_index)) {                                   ## for each diagnosed individual

    index_j = phenome_censor_index[j]
    row_data <- (unlist(as.vector(uk_module_date[index_j, ])))
    valid_dates <- as.Date(row_data, format = "%Y-%m-%d")
    recent_diag_date[index_j] <- max(valid_dates, na.rm = TRUE)

    setTxtProgressBar(pb, j)

  }

  recent_diag_date = as.Date(recent_diag_date, origin = "1970-01-01")
  status_ind[which(!is.na(recent_diag_date))] = 0

  ## input phecode time (for max censor time)
  surv_phe = cbind(uk_use, event_date, recent_diag_date, status_ind)
  print(head(surv_phe))
  surv_phe = surv_phe[which(!is.na(surv_phe[, "status_ind"])),]
  surv_phe$app_month = (as.yearmon(surv_phe$app_date)-as.yearmon(surv_phe$birth_date))*12                 ## birth to application pause time
  surv_phe$life_month = (as.yearmon(surv_phe$death_date)-as.yearmon(surv_phe$birth_date))*12              ## birth to death time
  surv_phe$lost_month = (as.yearmon(surv_phe$lost_date)-as.yearmon(surv_phe$birth_date))*12               ## birth to lost time
  surv_phe$recent_month = (as.yearmon(surv_phe$recent_diag_date)-as.yearmon(surv_phe$birth_date))*12      ## birth to the most recent diagnosed time (censored)
  surv_phe$attend_month = (as.yearmon(surv_phe$attend_date)-as.yearmon(surv_phe$birth_date))*12           ## birth to attend assessment time
  surv_phe$event_month = (as.yearmon(surv_phe$event_date)-as.yearmon(surv_phe$birth_date))*12             ## birth to event diagnosed time

  print(summary(surv_phe))

  ## input phecode event time
  input_event = which(surv_phe$status_ind == 1)
  summary(surv_phe$app_month[input_event] > surv_phe$event_month[input_event])                       ## censor time > event time
  summary(surv_phe$attend_month[input_event] < surv_phe$event_month[input_event])                    ## attend time > event time (maybe some individuals were diagnosed before entry)


  all.equal(which(!is.na(surv_phe$event_date)), which(is.na(surv_phe$recent_diag_date)))
  all.equal(which(is.na(surv_phe$event_date)), which(!is.na(surv_phe$recent_diag_date)))
  input_censor = which(surv_phe$status_ind == 0)
  surv_phe$event_month[input_censor] = apply(surv_phe[input_censor, 10:13], 1, function(x){min(x, na.rm = T)})

  ## some individuals have no diagnosis/death/lost-to-follow-up records
  ## use the time from birth to attendance as their survival time
  inf_index = which(is.infinite(surv_phe$recent_diag_date) & is.na(surv_phe$death_date) & is.na(surv_phe$lost_date))
  all.equal(surv_phe$event_month[inf_index], surv_phe$app_month[inf_index])
  surv_phe$event_month[inf_index] = surv_phe$attend_month[inf_index]

  ## or, just remove them
  # surv_phe = surv_phe[-inf_index,]

  surv_phe = surv_phe[, c(1, 1, 9, 15)]
  colnames(surv_phe) = c("FID", "IID", "status", "time")
  fwrite(surv_phe, file = paste0(output_path, phecode_id, "_pheSurv_useAttendTime.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

}
