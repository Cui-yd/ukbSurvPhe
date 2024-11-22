### Declaration: This script is sourced from the project createUKBphenome.
### Available at https://github.com/umich-cphds/createUKBphenome. 
### It is licensed under the GPL-3.0 license, and includes minor modifications by the current author.


# Input -----------------------------------------------------------------

options(stringsAsFactors=F)
suppressPackageStartupMessages(library("data.table",quietly = T))
suppressPackageStartupMessages(library("tidyr",quietly = T))
suppressPackageStartupMessages(library("parallel",quietly = T))
suppressPackageStartupMessages(library("intervals",quietly = T))
suppressPackageStartupMessages(library("XML",quietly = T))
suppressPackageStartupMessages(library("RCurl",quietly = T))
suppressPackageStartupMessages(library("rlist",quietly = T))
suppressPackageStartupMessages(library("bitops",quietly = T))

rm(list = ls())


setDTthreads(detectCores()/2)



### change this path by your own UKBB data
UKB_mainDat = "PATH/TO/ukb_main.csv"


## Some additional functions (also from createUKBphenome)
source("./scripts/function.getRanges.r")
source("./scripts/function.expandPhecodes.r")
source("./scripts/function.simpleCap.r")
source("./scripts/function.harmonizeICD9.r")
# source("./scripts/function.harmonizeICD10.r")


## List of people who withdrew
file.withdrawn <- sort(list.files("./data","w[0-9]+_.[0-9\\-]+.csv",full.name=T),decreasing=T)[1]
if(!file.exists(file.withdrawn)) {
  warning("Please add file with withdrawn samples 'w***_***.csv' to './data/' folder'!")
} else {
  withdrawn <- readLines(file.withdrawn)
}


## .fam file with genetic sex
file.fam <- sort(list.files("./data",".+\\.fam$",full.name=T),decreasing=T)[1]
if(!file.exists(file.fam)) {
  warning("Please add FAM file from genotype data '***.fam' to './data/' folder'!")
} else {
  fam <- fread(file.fam, select = c(1,5), header = FALSE)
  sex <- fam[,.(f.eid = V1, Sex = 2 - V5 ) ]
}


## this fold will save output files
if(!file.exists("./results")) dir.create("./results")   



## Prepare PheCode 
## read PheWAS map (downloaded from https://phewascatalog.org/phecodes; selected "export all" top right corner)
icd9map <- fread("./data/phecode_icd9_rolled.csv",colClasses="character")
## read PheWAS map (downloaded from https://phewascatalog.org/phecodes_icd10; selected "export all" top right corner)
icd10map <- fread("./data/phecode_icd10.csv",colClasses="character")


## collected phecode information by combining tables from https://github.com/PheWAS/ and the mapping tables
## Phecode information
load(file="./data/pheinfo.rda")
pheinfo <- data.table(pheinfo)
pheinfo <- pheinfo[,.(phecode,description,groupnum,group,color)]
## Sex rules: Both, Male, and Female
load(file="./data/sex_restriction.RData")
sex_restriction <- data.table(sex_restriction)
sex_restriction <- rbind(
  data.table('phecode'=sex_restriction[male_only == F & female_only == F,phecode],'sex'="Both"),
  data.table('phecode'=sex_restriction[male_only == T,phecode],'sex'="Male"),
  data.table('phecode'=sex_restriction[female_only == T,phecode],'sex'="Female")
)




# Mapping PheCode (ICD9 ICD10) -----------------------------------------------

## Exclusion / roll up rules
pheinfo2 <- unique(fread("./data/phecode_icd9_rolled.csv",colClasses="character",
                         select=c("PheCode","Excl. Phecodes","Excl. Phenotypes","Rollup","Leaf"),
                         col.names=c("phecode","phecode_exclude_range","phecode_exclude_phenotypes","rollup","leaf")))

## create one data.table with all criteria
pheinfo <- merge(pheinfo,pheinfo2,by="phecode")
pheinfo <- merge(pheinfo,sex_restriction,by="phecode")
pheinfo <- pheinfo[,c("phecode","description","sex","rollup","leaf","groupnum","group","color","phecode_exclude_range","phecode_exclude_phenotypes")]

## add phecode information that's missing (collected form earlier versions)
pheinfoOLD <- fread("./data/Phecode_Definitions_FullTable_Modified.txt",colClasses="character")
pheinfo <- rbind(pheinfo,pheinfoOLD[!phecode %in% pheinfo$phecode,])

## Phecode that should not be rolled up
norollup <- pheinfo$phecode[which(pheinfo$rollup == 0)]
## add manual no rollup rules:
# norollup <- c(norollup,fread("./data/no_rollup_extra.txt",colClasses="character")$phecode)



## Map ICD9 codes to PheCodes ####
message("Mapping ICD9 codes to PheCodes...")

## read UKB coding
ICD9codes <- fread("./data/icd9.tsv")
ICD9codes[!grepl("Block",coding),ICD9:=sapply(coding,harmonizeICD9)]
codeICD9 <- ICD9codes[,sort(unique(ICD9))]

mappedICD9Codes <- NULL
icd9map_new <- list()
for(i in 1:nrow(icd9map)){
	mapped9 <- grep(icd9map$ICD9[i],codeICD9)
	if(length(mapped9) > 0) {
		mappedICD9Codes <- unique(c(mappedICD9Codes,codeICD9[mapped9]))
		icd9map_new[[icd9map$ICD9[i]]] <- data.table(
			'phecode'=icd9map$PheCode[i],
			'ICD9'=codeICD9[mapped9])
	}
}
icd9key <- rbindlist(icd9map_new)
icd9unmapped <- codeICD9[!codeICD9 %in% mappedICD9Codes]

## roll up PheWAS codes
icd9key$added <- "original"
pcodes <- unique(c(icd9key$phecode,gsub("\\..+","",icd9key$phecode),gsub("(\\..).","\\1",icd9key$phecode)))
pcodes <- sort(pcodes)

for(p in 1:length(pcodes)){
    if(grepl("\\.",pcodes[p])) {
        iSub <- which(icd9key$phecode %in% pcodes[which(grepl(paste0("^",pcodes[p]),pcodes)
        				& nchar(pcodes) > nchar(pcodes[p]))]
        				& !icd9key$phecode %in% norollup)
    } else {
        iSub <- which(icd9key$phecode %in% pcodes[grep(paste0("^",pcodes[p],"\\."),pcodes)]
        				& !icd9key$phecode %in% norollup)
    }
    if(length(iSub) == 0) next
    iTop <- icd9key$ICD9[which(icd9key$phecode == pcodes[p])]
    addTop <- which(icd9key$ICD9 %in% unique(icd9key$ICD9[iSub]) & ! icd9key$ICD9 %in% iTop)
    if(length(addTop) == 0) next

    addKey <- icd9key[addTop,]
    addKey$phecode <- pcodes[p]
    addKey$added <- "rolled up PheWAS code"
    icd9key <- rbind(icd9key,addKey)
}

## Add ICD code description and phecode description
icd9key <- unique(icd9key)
print("Rollup of phewas codes (ICD9 code)")
print(table(icd9key$added))

icd9key <- merge(icd9key,ICD9codes[,.(ICD9,meaning,node_id,parent_id,selectable)],by="ICD9")
icd9key <- merge(icd9key,pheinfo,by="phecode")
icd9key <- icd9key[,c("ICD9","meaning","node_id","parent_id","selectable","phecode","description",
	"group","groupnum","added","sex","rollup","leaf")]

## remove any issues you might have found
# icd9map_remove <- fread("./data/remove_icd9_map.txt",colClasses="character")
# icd9key <- icd9key[!paste(ICD9,phecode,sep=":") %in% icd9map_remove[,paste(ICD9,phecode,sep=":")],]


## Map ICD10 codes to PheCodes ####
message("Mapping ICD10 codes to PheCodes")

## read UKB coding
ICD10codes <- fread("./data/icd10.tsv")
ICD10codes <- ICD10codes[!grepl("Block",coding),]
ICD10codes[,ICD10category:=gsub("([A-Z][0-9]{2}).+","\\1",coding)]
ICD10codes[,ICD10suffix:=gsub("[A-Z].+","",gsub("^[A-Z][0-9]{2}","",coding))]
ICD10codes[,ICD10:=paste0(ICD10category,ifelse(ICD10suffix == "","","."),ICD10suffix)]
ICD10codes <- ICD10codes[,c("ICD10category","ICD10suffix"):=NULL]

codeICD10 <- ICD10codes[,sort(unique(ICD10))]
mappedICD10Codes <- NULL
icd10map_new <- list()
for(i in 1:nrow(icd10map)){
	mapped10 <- grep(icd10map$ICD10[i],codeICD10)
	if(length(mapped10) > 0) {
		mappedICD10Codes <- unique(c(mappedICD10Codes,codeICD10[mapped10]))
		icd10map_new[[icd10map$ICD10[i]]] <- data.table(
			'phecode'=icd10map$PheCode[i],
			'ICD10'=codeICD10[mapped10])
	}
}
icd10key <- rbindlist(icd10map_new)
icd10unmapped <- codeICD10[!codeICD10 %in% mappedICD10Codes]

## roll up phewas codes
icd10key$added <- "original"
pcodes <- unique(c(icd10key$phecode,gsub("\\..+","",icd10key$phecode),gsub("(\\..).","\\1",icd10key$phecode)))
pcodes <- sort(pcodes)
for(p in 1:length(pcodes)){
    if(grepl("\\.",pcodes[p])) {
        iSub <- which(icd10key$phecode %in% pcodes[which(grepl(paste0("^",pcodes[p]),pcodes)
        				& nchar(pcodes) > nchar(pcodes[p]))]
        				& !icd10key$phecode %in% norollup)
    } else {
        iSub <- which(icd10key$phecode %in% pcodes[grep(paste0("^",pcodes[p],"\\."),pcodes)]
        				& !icd10key$phecode %in% norollup)
    }
    if(length(iSub) == 0) next
    iTop <- icd10key$ICD10[which(icd10key$phecode == pcodes[p])]
    addTop <- which(icd10key$ICD10 %in% unique(icd10key$ICD10[iSub]) & ! icd10key$ICD10 %in% iTop)
    if(length(addTop) == 0) next

    addKey <- icd10key[addTop,]
    addKey$phecode <- pcodes[p]
    addKey$added <- "rolled up PheWAS code"
    icd10key <- rbind(icd10key,addKey)
}
icd10key <- unique(icd10key)

print("Rollup of phewas codes (ICD10 code)")
print(table(icd10key$added))

icd10key <- merge(icd10key,ICD10codes[,.(ICD10,meaning,node_id,parent_id,selectable)],by="ICD10")
icd10key <- merge(icd10key,pheinfo,by="phecode")
icd10key <- icd10key[,c("ICD10","meaning","node_id","parent_id","selectable","phecode","description",
                        "group","groupnum","added","sex","rollup","leaf")]

## remove any issues you might have found
# icd10map_remove <- fread("./data/remove_icd10_map.txt",colClasses="character")
# icd10key <- icd10key[!paste(ICD10,phecode,sep=":") %in% icd10map_remove[,paste(ICD10,phecode,sep=":")],]

print(paste0("Mapping tables created and stored in ./results/UKB_PHENOME_ICD*_PHECODE_MAP",".txt"))
fwrite(icd9key, paste0("./results/UKB_PHENOME_ICD9_PHECODE_MAP",".txt"),sep="\t",quote=F)
fwrite(icd10key, paste0("./results/UKB_PHENOME_ICD10_PHECODE_MAP",".txt"),sep="\t",quote=F)
save.image(file = paste0("./results/UKB_MAPPING", ".RData"))



# Create UKBB Phenome -----------------------------------------------------

print("Collecting UKB fields and data")


## get all field ids using scripts/function.collectFields.r; the Field Listings store in file `./data/FieldListing.txt`
fields <- fread("./data/FieldListing.txt")
fieldids9 <- fields[grepl("ICD9",Description) & !grepl("addendum|date",Description,ignore.case=T),`Field ID`]
fieldids10 <- fields[grepl("ICD10",Description) & !grepl("addendum|date",Description,ignore.case=T),`Field ID`]


ukb_data_columns <- colnames(fread(UKB_mainDat, header = T, sep=",", nrows = 1))

## ICD9 columns in UKBB data
extractICD9columns <- grep(paste(fieldids9, sep = "", collapse = "|"),ukb_data_columns)
ICD9columnsname <- ukb_data_columns[c(1, extractICD9columns)]
ICD9Data <- fread(UKB_mainDat, header = T, select = ICD9columnsname)
ICD9 <- na.omit(gather(ICD9Data, "field", "ICD9", -eid))
rm(ICD9Data)
gc()

## ICD10 columns in UKBB data
extractICD10columns <- grep(paste(fieldids10, sep = "", collapse = "|"),ukb_data_columns)
ICD10columnsname <- ukb_data_columns[c(1, extractICD10columns)]
ICD10Data <- fread(UKB_mainDat, header = T, select = ICD10columnsname)
ICD10 <- na.omit(gather(ICD10Data, "field", "ICD10", -eid))
rm(ICD10Data)
gc()

## 31-0.0 sex; 22001-0.0 genetic sex
extractSexColumns <- which(ukb_data_columns %in% c("31-0.0","22001-0.0"))
SexColumnsname <- ukb_data_columns[c(1, extractSexColumns)]
SexDate <- fread(UKB_mainDat, header = T, select = SexColumnsname)
setnames(SexDate,c("31-0.0","22001-0.0"),c("Sex","GeneticSex"))
SexDate[,Sex:=ifelse(Sex != GeneticSex, NA, Sex)]
SexDate <- na.omit(SexDate)   # remove "Sex != GeneticSex"
females <- SexDate[Sex==0, eid]
males <- SexDate[Sex==1, eid]


## Remove field (variants names)
ICD9 <- as.data.table(ICD9)
ICD10 <- as.data.table(ICD10)
ICD9 <- unique(ICD9[,field:=NULL])
ICD10 <- unique(ICD10[,field:=NULL])


## Cleanup ICD9 codes; i.e. add dot and remove suffixes
ICD9[,ICD9har:=sapply(ICD9,harmonizeICD9)]
ICD9 <- unique(ICD9)
setnames(ICD9, c("ICD9", "ICD9har"), c("ICD9old", "ICD9"))


## Cleanup ICD10 codes; i.e. add dot and remove suffixes
ICD10[,ICD10category:=gsub("([A-Z][0-9]{2}).+","\\1",ICD10)]
ICD10[,ICD10suffix:=gsub("[A-Z].+","",gsub("^[A-Z][0-9]{2}","",ICD10))]
ICD10[,ICD10new:=paste0(ICD10category,ifelse(ICD10suffix == "","","."),ICD10suffix)]
ICD10[,c("ICD10category","ICD10suffix"):=NULL]
ICD10 <- unique(ICD10)
setnames(ICD10, c("ICD10", "ICD10new"), c("ICD10old", "ICD10"))


## Merge ICD9 with phewas map then remove ICD9 codes
phecode1 <- merge(ICD9,icd9key[phecode != "",.(phecode,ICD9)],by="ICD9",allow.cartesian=T)
phecode1[,ICD9:=NULL]
phecode1[,ICD9old:=NULL]

## Merge ICD10 with phewas map then remove ICD10 codes
phecode2 <- merge(ICD10,icd10key[phecode != "",.(phecode,ICD10)],by="ICD10",allow.cartesian=T)
phecode2[,ICD10:=NULL]
phecode2[,ICD10old:=NULL]


## combine ICD9 and ICD10-based phecodes
phenotypeData <- unique(rbind(phecode1,phecode2))
phecodes <- unique(phenotypeData$phecode)
phecodes <- phecodes[order(as.numeric(phecodes))]

## collect inclusions (at this point they are identical)
inclusions <- exclusions <- phenotypeData[,c("eid","phecode")]
inclusions <- exclusions <- split(inclusions$eid,inclusions$phecode)

## prepare summary table
pheinfo2 <- data.table(pheinfo[phecode %in% phecodes,
                               c("phecode","description","group","groupnum","color",
                                 "sex","phecode_exclude_range")],
                       'ncases'=0,'ncontrols'=0, 'case:control'=0)

## only keep samples with valid sex (see above)
sampleNames <- SexDate$eid

## create two empty data.frames (one with exclusion criteria applied, one without apply exclusion criteria)
## add X to phecodes
phenoOut <- phenoOut0 <- data.table('IID'=sampleNames,matrix(NA,ncol=length(phecodes),nrow=length(sampleNames)))
colnames(phenoOut) <- c("IID",paste0("X",phecodes))


## collect cases with potentially wrong sex
notFemale <- list()
notMale <- list()

## TRUE FALSE variable for sex
femaleTF <- sampleNames %in% females 
maleTF <- sampleNames %in% males 

message("Creating case control studies")
## create case control studies, one phecode at a time
pb = txtProgressBar(min = 0, max = nrow(pheinfo2), initial = 0, style = 3) 
for(p in 1:nrow(pheinfo2)){
  phecode_remove <- ""
  phecode <- pheinfo2$phecode[p]
  
  ## collect phecodes to include from controls
  exclude_phecodes <- phecode
  if(pheinfo2$phecode_exclude_range[p] != ""){
    phecode_remove <- unlist(strsplit(gsub(" ","",pheinfo2$phecode_exclude_range[p]),",")[[1]])
    exclude_phecodes <- c(exclude_phecodes,unlist(sapply(phecode_remove,function(x) expandPhecodes(x,T),USE.NAMES=F)))
  }
  exclude_phecodes <- unique(exclude_phecodes[which(exclude_phecodes %in% pheinfo2$phecode)])
  
  ## collected sample sets by inclusion and exclusion criteria
  pinclusions <- unlist(inclusions[phecode])
  pexclusions <- unique(unlist(exclusions[exclude_phecodes]))
  
  ## phecode based case control definitions (without sex filter)
  cases <- sampleNames %in% pinclusions
  controls <- !sampleNames %in% pexclusions
  
  ## case control status
  ccstatus <- ccstatus0 <- rep(NA,nrow(phenoOut))
  if(pheinfo2$sex[p] == "Both"){
    ## non-sex-specific traits
    ccstatus[controls] <- 0
    ccstatus[cases] <- ccstatus0[cases] <- 1
    ccstatus0[!cases] <- 0
  } else if (pheinfo2$sex[p] == "Female"){
    ## female-specific traits
    ccstatus[controls & femaleTF] <- 0
    ccstatus[cases & femaleTF] <- ccstatus0[cases & femaleTF] <- 1
    ccstatus0[!cases & femaleTF] <- 0
    checkFemales <- sampleNames[cases & !femaleTF]
    if(length(checkFemales) > 0) notFemale[[phecode]] <- data.table('IID'=checkFemales,phecode)
  } else if (pheinfo2$sex[p] == "Male"){
    ## male-specific traits
    ccstatus[controls & maleTF] <- 0
    ccstatus[cases & maleTF] <- ccstatus0[cases & maleTF] <- 1
    ccstatus0[!cases & maleTF] <- 0
    checkMales <- sampleNames[cases & !maleTF]
    if(length(checkMales) > 0) notMale[[phecode]] <- data.table('IID'=checkMales,phecode)
  }
  phenoOut[[paste0("X",phecode)]] <- ccstatus
  phenoOut0[[paste0("X",phecode)]] <- ccstatus0
  
  ## get sample sizes
  pheinfo2$ncontrols[p] <- length(which(ccstatus == 0))
  pheinfo2$ncases[p] <- length(which(ccstatus == 1))
  pheinfo2$`case:control`[p] <- paste0("1:", round(length(which(ccstatus == 0))/length(which(ccstatus == 1)), 2))
  setTxtProgressBar(pb,p)
}

message("Write output:")
## summary table
fwrite(pheinfo2,paste0("results/UKB_PHENOME_DESCRIPTION.txt"),sep="\t",row.names=F,col.names=T,quote=F)

## phenome
fwrite(phenoOut,paste0("results/UKB_PHENOME.txt"),sep="\t",row.names=F,col.names=T,quote=F)

## phenome without exclusion criteria (sex filter was applied)
fwrite(phenoOut0,paste0("results/UKB_PHENOME_NO_EXCLUSIONS.txt"),sep="\t",row.names=F,col.names=T,quote=F)

## Output individuals with diagnoses that don't match their sex
notFemale <- rbindlist(notFemale)
notFemale <- merge(notFemale,pheinfo[,.(phecode,description)],by="phecode")
notFemale[,Sex:="Male"]

notMale <- rbindlist(notMale)
notMale <- merge(notMale,pheinfo[,.(phecode,description)],by="phecode")
notMale[,Sex:="Female"]

sexIssues <- rbind(notFemale,notMale)
sexIssues <- sexIssues[order(IID,phecode),]
fwrite(sexIssues,paste0("./results/UKB_PHENOME_SEXMISMATCH_SEX-SPECIFIC_DISAGNOSES.txt"),
       sep="\t",row.names=F,col.names=T,quote=F)

## Output unmapped ICD9 codes
ICD9_U <- ICD9[ICD9 %in% icd9unmapped,]
ICD9_U <- data.table(table(ICD9_U$ICD9))
colnames(ICD9_U) <- c("ICD9","N")
ICD9_U <- merge(ICD9_U,ICD9codes,by="ICD9")[order(N,decreasing=T),]
fwrite(ICD9_U,paste0("./results/UKB_PHENOME_UNMAPPED_ICD9_CODES.txt"),
       sep="\t",row.names=F,col.names=T,quote=F)

## Output unmapped ICD10 codes
ICD10_U <- ICD10[ICD10 %in% icd10unmapped,]
ICD10_U <- data.table(table(ICD10_U$ICD10))
colnames(ICD10_U) <- c("ICD10","N")
ICD10_U <- merge(ICD10_U,ICD10codes,by="ICD10")[order(N,decreasing=T),]
fwrite(ICD10_U,paste0("./results/UKB_PHENOME_UNMAPPED_ICD10_CODES.txt"),
       sep="\t",row.names=F,col.names=T,quote=F)

print(paste0("Phenome created and stored in ./results/UKB_PHENOME*.txt"))


