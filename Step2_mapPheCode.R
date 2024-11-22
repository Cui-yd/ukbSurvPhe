### Declaration: This script is sourced from the project ukbb_pan_ancestry
### The sourced script is available at https://github.com/atgu/ukbb_pan_ancestry/blob/master/create_phecode_mapping.R. 
### It is licensed under the MIT license, and includes minor modifications by the current author.


suppressPackageStartupMessages(library(dplyr, quietly = T))
suppressPackageStartupMessages(library(stringr, quietly = T))
suppressPackageStartupMessages(library(purrr, quietly = T))
suppressPackageStartupMessages(library(data.table, quietly = T))


icd9key = fread("./results/UKB_PHENOME_ICD9_PHECODE_MAP.txt", colClasses = "character", data.table = F) 
icd10key = fread("./results/UKB_PHENOME_ICD10_PHECODE_MAP.txt", colClasses = "character", data.table = F) 

icd9_map = icd9key %>%
  group_by(phecode, sex, description, group) %>%
  summarize(
    icd9_code = str_c(str_remove(ICD9, "\\."), collapse = ',')
  ) %>%
  ungroup %>%
  arrange(phecode)

icd10_map = icd10key %>%
  group_by(phecode, sex, description, group) %>%
  summarize(
    icd10_code = str_c(str_remove(ICD10, "\\."), collapse = ',')
  ) %>%
  ungroup %>%
  arrange(phecode)

icd_all = merge(icd9_map, icd10_map, all = TRUE)

if (length(unique(icd_all$phecode)) != nrow(icd_all)) {
  stop("Inconsistent sex def.")
}

pheinfo = fread("./results/UKB_PHENOME_DESCRIPTION.txt", colClasses = "character", data.table = F)
icd_all = subset(icd_all, phecode %in% pheinfo$phecode) %>% arrange(phecode)

# source: https://github.com/umich-cphds/createUKBphenome/blob/master/scripts/function.expandPhecodes.r
expandPhecodes <- function(x,addIntegers=T){
  if(is.na(x) | x == "") return(NA)
  if(grepl("\\-",x)){
    # split range
    # character prefix
    i1 <- strsplit(x,"-")[[1]]
    
    # numeric length of digits before "."
    nprefix <- max(nchar(gsub("\\..+","",i1)))
    # numbers of digits
    ndigits <- max(c(nchar(gsub("^[0-9]+[\\.]{0,1}","",i1)),0))
    # add "." to length of formatted number if present
    addDot <- max(as.numeric(grepl("\\.",i1)))
    # create sequence of ICD codes
    seq1 <- seq(as.numeric(i1[1]),as.numeric(i1[2]),(1/10^ndigits))
    # format sequence to match intput
    seq1 <- formatC(seq1, format='f', digits=ndigits,width=nprefix+ndigits+addDot,flag=0)
    # add integers if within range
    if(addIntegers) seq1 <- unique(sort(c(seq1,gsub("\\..+","",seq1[which(round(as.numeric(seq1)) == as.numeric(seq1))]))))
    
    if(ndigits == 2){
      seq2 <- seq(as.numeric(i1[1]),as.numeric(i1[2]),(1/10^(ndigits-1)))
      seq2 <- formatC(seq2, format='f', digits=ndigits-1,width=nprefix+ndigits+addDot-1,flag=0)
      seq1 <- unique(sort(c(seq1,seq2)))
    }
    return(seq1)
  } else {
    return(x)
  }
}

icd_all$exclude_phecodes = map_chr(1:nrow(pheinfo), function(p) {
  phecode_remove <- ""
  phecode <- pheinfo$phecode[p]
  
  # collect phecodes to include from controls
  exclude_phecodes <- phecode
  if(pheinfo$phecode_exclude_range[p] != ""){
    phecode_remove <- unlist(strsplit(gsub(" ","",pheinfo$phecode_exclude_range[p]),",")[[1]])
    exclude_phecodes <- c(exclude_phecodes,unlist(sapply(phecode_remove,function(x) expandPhecodes(x,T),USE.NAMES=F)))
  }
  exclude_phecodes <- unique(exclude_phecodes[which(exclude_phecodes %in% pheinfo$phecode)])
  return(str_c(sort(exclude_phecodes), collapse=","))
})

write.table(icd_all, paste0("./results/UKB_PHECODE_v1.2b1_ICD_Mapping.txt"), quote = F, row.names = F, sep = "\t")

