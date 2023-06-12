#BMI and ECG parameters
#M Ardissino 
#Date: 16/9/2022

####Load packages####
library(data.table)
library(dplyr)
library(tidyr)
library(foreign)
library(tibble)
library(TwoSampleMR)
library(meta)
library(MendelianRandomization)
library(meta)
library(metafor)
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(phenoscanner)
library(beepr)
library(MRInstruments)
library(MRPRESSO)
library(R.utils)
library(doMC)
library(stringr)
library(nlmr)
library(PolyMR)
library(forestplot)
library(readxl)
library(openxlsx)
library(MVMR)


#### Make IV sets ####
# Meta-analysis of bmi in UK Biobank and GIANT data.max N = 806,834
# Beta: effect per 1-unit increase in BMI 
bmi <- as.data.frame(fread("https://portals.broadinstitute.org/collaboration/giant/images/1/14/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"))
head(bmi)
bmi$SNP<-gsub('.{4}$', '', bmi$SNP) #remove last 4 characters from SNP column (currently rs28788874:c:t rather than just rsID)
head(bmi)
setnames(bmi, old=c("SNP", "Tested_Allele","Other_Allele",  "BETA", "SE", "P", "CHR",  "POS", "Freq_Tested_Allele", "N"),
         new=c("SNP", "bmi_ea", "bmi_nea", "bmi_beta", "bmi_se", "bmi_p", "bmi_chr", "bmi_pos", "bmi_eaf", 'bmi_ss'))
bmi <- bmi[,c("SNP", "bmi_ea", "bmi_nea", "bmi_beta", "bmi_se", "bmi_p", "bmi_chr", "bmi_pos", "bmi_eaf", 'bmi_ss')]
bmi$bmi_chr<-as.numeric(bmi$bmi_chr)
bmi$bmi_ea<-toupper(bmi$bmi_ea)
bmi$bmi_p<-as.numeric(bmi$bmi_p)
bmi$bmi_eaf<-as.numeric(bmi$bmi_eaf)
write.table(bmi, '/Volumes/MADDY2/mrdata/bmi_mvmrformat.tsv', row.names = FALSE, col.names = TRUE)
bmi_iv<-as.data.frame(filter(bmi, bmi_p<5e-8))
rm(bmi)

# WHR 
whr <- as.data.frame(fread("https://portals.broadinstitute.org/collaboration/giant/images/6/6e/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"))
head(whr)
whr$SNP<-gsub('.{4}$', '', whr$SNP) #remove last 4 characters from SNP column (currently rs28788874:c:t rather than just rsID)
head(whr)
setnames(whr, old=c("SNP", "Tested_Allele","Other_Allele",  "BETA", "SE", "P", "CHR",  "POS", "Freq_Tested_Allele", "N"),
         new=c("SNP", "whr_ea", "whr_nea", "whr_beta", "whr_se", "whr_p", "whr_chr", "whr_pos", "whr_eaf", 'whr_ss'))
whr <- whr[,c("SNP", "whr_ea", "whr_nea", "whr_beta", "whr_se", "whr_p", "whr_chr", "whr_pos", "whr_eaf", 'whr_ss')]
whr$whr_chr<-as.numeric(whr$whr_chr)
whr$whr_ea<-toupper(whr$whr_ea)
whr$whr_p<-as.numeric(whr$whr_p)
whr$whr_eaf<-as.numeric(whr$whr_eaf)
write.table(whr, '/Volumes/MADDY2/mrdata/whr_mvmrformat.tsv', row.names = FALSE, col.names = TRUE)
whr_iv<-as.data.frame(filter(whr, whr_p<5e-8))
rm(whr)

# Height 
heig <- as.data.frame(fread("https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz"))
setnames(heig, old=c("SNP", "Tested_Allele","Other_Allele",  "BETA", "SE", "P", "CHR",  "POS", "Freq_Tested_Allele_in_HRS", "N"),
         new=c("SNP", "heig_ea", "heig_nea", "heig_beta", "heig_se", "heig_p", "heig_chr", "heig_pos", "heig_eaf", 'heig_ss'))
heig <- heig[,c("SNP", "heig_ea", "heig_nea", "heig_beta", "heig_se", "heig_p", "heig_chr", "heig_pos", "heig_eaf", 'heig_ss')]
heig$heig_chr<-as.numeric(heig$heig_chr)
heig$heig_ea<-toupper(heig$heig_ea)
heig$heig_p<-as.numeric(heig$heig_p)
heig$heig_eaf<-as.numeric(heig$heig_eaf)
write.table(heig, '/Volumes/MADDY2/mrdata/heig_mvmrformat.tsv', row.names = FALSE, col.names = TRUE)
heig_iv<-as.data.frame(filter(heig, heig_p<5e-8))
rm(heig)


# WBFM  WBLM and weight
# Load Neale reference set to merge RSIDs to Neale files
refgen<-fread("/Volumes/MADDY2/mrdata/ukbrefsmall.csv") # link = https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz
refgen$SNP<-refgen$rsid
refgen<-refgen[,c("variant", "SNP")]

# Whole body fat mass, link: https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23100_irnt.gwas.imputed_v3.both_sexes.tsv.bgz , n=354244, continuous
gunzip("~/Downloads/23100_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "/Volumes/MADDY2/mrdata/wbfm.tsv")
wbfm_iv <- as.data.frame(fread("/Volumes/MADDY2/mrdata/wbfm.tsv"))
head(wbfm_iv)
wbfm_iv<-merge(wbfm_iv, refgen, by=c("variant"), all.x=T, all.y=F) #merge rsid
wbfm_iv$variant<- sub(".+?:","",wbfm_iv$variant) #remove chr from variant column
wbfm_iv$variant<- sub(".+?:","",wbfm_iv$variant) #remove  pos from variant column
wbfm_iv$wbfm_ref<-sub(":.*","",wbfm_iv$variant) # make ref with the  allele after :
wbfm_iv$wbfm_ea<-sub(".*:","",wbfm_iv$variant) # make ea with the  allele before :
wbfm_iv$wbfm_maf <- ifelse(wbfm_iv$minor_allele == wbfm_iv$wbfm_ea, wbfm_iv$minor_AF, 1-wbfm_iv$minor_AF)
setnames(wbfm_iv, old=c("SNP", "wbfm_ea", "wbfm_ref", "beta", "se", "pval","wbfm_maf"),
         new=c("SNP", "wbfm_ea", "wbfm_nea", "wbfm_beta", "wbfm_se", "wbfm_p",  "wbfm_eaf"))
wbfm_iv<-wbfm_iv[,c("SNP", "wbfm_ea", "wbfm_nea", "wbfm_beta", "wbfm_se", "wbfm_p",  "wbfm_eaf")]
wbfm_iv$wbfm_ss <- 354244
write.table(wbfm_iv, '/Volumes/MADDY2/mrdata/wbfm_mvmrformat.tsv', row.names = FALSE, col.names = TRUE)
wbfm_iv<-subset(wbfm_iv, wbfm_p<5e-8)
rownames(wbfm_iv) <- c()
unlink("/Volumes/MADDY2/mrdata/wbfm.tsv")

# Whole body fat free (lean) mass, link: https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23101_irnt.gwas.imputed_v3.both_sexes.tsv.bgz , n=354808, continuous
gunzip("~/Downloads/23101_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "/Volumes/MADDY2/mrdata/wblm.tsv")
wblm_iv <- as.data.frame(fread("/Volumes/MADDY2/mrdata/wblm.tsv"))
wblm_iv<-merge(wblm_iv, refgen, by=c("variant"), all.x=T, all.y=F) #merge rsid
wblm_iv$variant<- sub(".+?:","",wblm_iv$variant) #remove chr from variant column
wblm_iv$variant<- sub(".+?:","",wblm_iv$variant) #remove  pos from variant column
wblm_iv$wblm_ref<-sub(":.*","",wblm_iv$variant) # make ref with the  allele after :
wblm_iv$wblm_ea<-sub(".*:","",wblm_iv$variant) # make ea with the  allele before :
wblm_iv$wblm_maf <- ifelse(wblm_iv$minor_allele == wblm_iv$wblm_ea, wblm_iv$minor_AF, 1-wblm_iv$minor_AF)
setnames(wblm_iv, old=c("SNP", "wblm_ea", "wblm_ref", "beta", "se", "pval","wblm_maf"),
         new=c("SNP", "wblm_ea", "wblm_nea", "wblm_beta", "wblm_se", "wblm_p",  "wblm_eaf"))
wblm_iv<-wblm_iv[,c("SNP", "wblm_ea", "wblm_nea", "wblm_beta", "wblm_se", "wblm_p",  "wblm_eaf")]
wblm_iv$wblm_ss <- 354808
write.table(wblm_iv, '/Volumes/MADDY2/mrdata/wblm_mvmrformat.tsv', row.names = FALSE, col.names = TRUE)
wblm_iv<-subset(wblm_iv, wblm_p<5e-8)
rownames(wblm_iv) <- c()
unlink("/Volumes/MADDY2/mrdata/wblm.tsv")

# Weight, link: https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23098_irnt.gwas.imputed_v3.both_sexes.tsv.bgz , n=360116, continuous
gunzip("~/Downloads/23098_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "/Volumes/MADDY2/mrdata/weigh.tsv")
weigh_iv <- as.data.frame(fread("/Volumes/MADDY2/mrdata/weigh.tsv"))
weigh_iv<-merge(weigh_iv, refgen, by=c("variant"), all.x=T, all.y=F) #merge rsid
weigh_iv$variant<- sub(".+?:","",weigh_iv$variant) #remove chr from variant column
weigh_iv$variant<- sub(".+?:","",weigh_iv$variant) #remove  pos from variant column
weigh_iv$weigh_ref<-sub(":.*","",weigh_iv$variant) # make ref with the  allele after :
weigh_iv$weigh_ea<-sub(".*:","",weigh_iv$variant) # make ea with the  allele before :
weigh_iv$weigh_maf <- ifelse(weigh_iv$minor_allele == weigh_iv$weigh_ea, weigh_iv$minor_AF, 1-weigh_iv$minor_AF)
setnames(weigh_iv, old=c("SNP", "weigh_ea", "weigh_ref", "beta", "se", "pval","weigh_maf"),
         new=c("SNP", "weigh_ea", "weigh_nea", "weigh_beta", "weigh_se", "weigh_p",  "weigh_eaf"))
weigh_iv<-weigh_iv[,c("SNP", "weigh_ea", "weigh_nea", "weigh_beta", "weigh_se", "weigh_p",  "weigh_eaf")]
weigh_iv$weigh_ss <- 360116
write.table(weigh_iv, '/Volumes/MADDY2/mrdata/weigh_mvmrformat.tsv', row.names = FALSE, col.names = TRUE)
weigh_iv<-subset(weigh_iv, weigh_p<5e-8)
rownames(weigh_iv) <- c()
unlink("/Volumes/MADDY2/mrdata/weigh.tsv")
rm(refgen)


genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls())))

# remove all iv sets with 1 or 0 snps only
to.rm <- unlist(eapply(.GlobalEnv, function(x) is.data.frame(x) && nrow(x)<2))
rm(list = names(to.rm)[to.rm], envir = .GlobalEnv)

##save all IVs as CSV files
rm(to.rm, genelist)
setwd("~/Desktop/bmi-ecg/ivs")
files <- mget(ls())
for (i in 1:length(files)){
  write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))
}
setwd("~/Desktop/bmi-ecg")
rm(list=ls())

#### Format outcome data ####
# Download file from: 
# # PR http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008042/
# # QRS: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008054/
# # QT: https://personal.broadinstitute.org/ryank/Nauffal_2022_QT_GWAS_SAIGE.zip 

print <- as.data.frame(fread("~/Downloads/pr_interval.combined.page.out"))
print$SNP<-print$rsid
setnames(print, old=c("rsid", "Effect-allele", "Other-allele", "Beta", "SE", "P-val", "Chr", "Position_hg19","Effect-allele-frequency", "Sample-size"),
         new=c("SNP", "print_ea", "print_nea", "print_beta", "print_se", "print_p", "chr", "pos","print_eaf", "print_ss"))
print<-print[, c("SNP", "print_ea", "print_nea", "print_beta", "print_se", "print_p", "chr", "pos","print_eaf", "print_ss")]
write.table(print, "out/print-woj.tsv", sep = "\t",quote = FALSE, col.names = TRUE,row.names = FALSE)
unlink("~/Downloads/pr_interval.combined.page.out")
rm(print)

qrsint <- as.data.frame(fread("~/Downloads/qrs_interval.combined.page.out"))
setnames(qrsint, old=c("rsid", "Effect-allele", "Other-allele", "Beta", "SE", "P-val", "Chr", "Position_hg19","Effect-allele-frequency", "Sample-size"),
         new=c("SNP", "qrsint_ea", "qrsint_nea", "qrsint_beta", "qrsint_se", "qrsint_p", "chr", "pos","qrsint_eaf", "qrsint_ss"))
qrsint<-qrsint[, c("SNP", "qrsint_ea", "qrsint_nea", "qrsint_beta", "qrsint_se", "qrsint_p", "chr", "pos","qrsint_eaf", "qrsint_ss")]
write.table(qrsint, "out/qrsint-woj.tsv", sep = "\t",quote = FALSE, col.names = TRUE,row.names = FALSE)
unlink("~/Downloads/qrs_interval.combined.page.out")
rm(qrsint)

pwd <- as.data.frame(fread("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004826/Pwave_duration_GWAS_ALL_maf0.01.txt.gz"))
setnames(pwd, old=c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "chr", "pos"),
         new=c("SNP", "pwd_ea", "pwd_nea", "pwd_beta", "pwd_se", "pwd_p", "chr", "pos"))
pwd<-pwd[, c("SNP", "pwd_ea", "pwd_nea", "pwd_beta", "pwd_se", "pwd_p", "chr", "pos")]
write.table(pwd, "out/pwd-chris.tsv", sep = "\t",quote = FALSE, col.names = TRUE,row.names = FALSE)
rm(pwd)
                       
# qtc nuffal
qtcint <- as.data.frame(fread("~/Downloads/Nauffal_2022_QT_GWAS_SAIGE/ukb_qtc_gwas_summary_statistics_final.txt"))
setnames(qtcint, old=c("ID_37", "Effect_allele", "Other_allele", "Beta", "SE", "P-value", "Chromosome", "Position_b37"),
         new=c("Markername", "qtcint_ea", "qtcint_nea", "qtcint_beta", "qtcint_se", "qtcint_p", "chr", "pos"))
qtcint<-qtcint[, c("Markername", "qtcint_ea", "qtcint_nea", "qtcint_beta", "qtcint_se", "qtcint_p", "chr", "pos")]
refgen<-fread("~/Desktop/MR_studies/MR_datasets/ukbrefsmall.csv")
refgen$SNP<-refgen$rsid
refgen$Markername <- str_c(refgen$chr, ":", refgen$pos)
refgen<-refgen[,c("Markername", "SNP")]
qtcint$Markername <- str_c(qtcint$chr, ":", qtcint$pos)
qtcint<-merge(qtcint, refgen, by=c("Markername"), all.x=T, all.y=F) #merge rsid
qtcint<-qtcint[, c("SNP", "qtcint_ea", "qtcint_nea", "qtcint_beta", "qtcint_se", "qtcint_p", "chr", "pos")]
write.table(qtcint, "out/qtcint-nuf.tsv", sep = "\t",quote = FALSE, col.names = TRUE,row.names = FALSE)
unlink("~/Downloads/Nauffal_2022_QT_GWAS_SAIGE/ukb_qtc_gwas_summary_statistics_final.txt")
rm(qtcint, refgen)

#####  OUTCOME 1 = printerval #### 
#load exposures -  all CSVs from IV folder 
setwd("~/Desktop/bmi-ecg/ivs")    #set wd to where i want to read IV sets from 
files = list.files(pattern="*.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets

instruments<- merge(bmi_iv, whr_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, heig_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wbfm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wblm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, weigh_iv[,-1], by="SNP", all=TRUE)

rm(files, data_list)
setwd("~/Desktop/bmi-ecg") #set wd back to project
#
#
# Extract outcome for all IVs - print
print <- as.data.frame(fread(("out/print-woj.tsv")))
print$print_ea<-toupper(print$print_ea)
print$print_nea<-toupper(print$print_nea)
colnames(print)

# merge with all exposures
instruments_print <- print[which(paste(print$SNP) %in% paste(instruments$SNP)),]
data <- merge(instruments, instruments_print, by=c("SNP"), all=FALSE)
data <- data[!duplicated(data$SNP),]
head(data)

# Separate 'data' file into exposure and outcome sets, and format respective dataframes for MR
outcomedata<-data[,c("SNP", "print_ea", "print_nea", "print_beta", "print_se", "print_p", "chr", "pos","print_eaf", "print_ss")]
outcomedat<-format_data(outcomedata,type="outcome",snp_col = "SNP",
                        beta_col = "print_beta",
                        se_col = "print_se",
                        eaf_col = "print_eaf",
                        effect_allele_col = "print_ea",
                        other_allele_col = "print_nea",
                        pos_col	="pos",
                        chr_col="chr",pval_col = "print_p", 
                        samplesize_col = "print_ss")

# BMI
exposuredata<-data[,c("SNP", "bmi_ea", "bmi_nea", "bmi_beta", "bmi_se", "bmi_p", "chr", "pos", "bmi_eaf", 'bmi_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "bmi_beta",
                         se_col = "bmi_se",
                         eaf_col = "bmi_eaf",
                         effect_allele_col = "bmi_ea",
                         other_allele_col = "bmi_nea",
                         pos_col	="pos",
                         chr_col="chr",pval_col="bmi_p", 
                         samplesize_col = "bmi_ss")
bmi_har<-harmonise_data(exposuredat, outcomedat, action = 2)
bmi_har<-bmi_har[which(bmi_har$mr_keep==TRUE),]

# WHR
exposuredata<-data[,c("SNP", "whr_ea", "whr_nea", "whr_beta", "whr_se", "whr_p", "chr", "pos", "whr_eaf", 'whr_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "whr_beta",
                         se_col = "whr_se",
                         eaf_col = "whr_eaf",
                         effect_allele_col = "whr_ea",
                         other_allele_col = "whr_nea",pos_col	="pos",
                         chr_col="chr",pval_col="whr_p", 
                         samplesize_col = "whr_ss")
whr_har<-harmonise_data(exposuredat, outcomedat, action = 2)
whr_har<-whr_har[which(whr_har$mr_keep==TRUE),]

# Heig
exposuredata<-data[,c("SNP", "heig_ea", "heig_nea", "heig_beta", "heig_se", "heig_p", "chr", "pos", "heig_eaf", 'heig_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "heig_beta",
                         se_col = "heig_se",
                         eaf_col = "heig_eaf",
                         effect_allele_col = "heig_ea",
                         other_allele_col = "heig_nea",pos_col	="pos",
                         chr_col="chr",pval_col="heig_p", 
                         samplesize_col = "heig_ss")
heig_har<-harmonise_data(exposuredat, outcomedat, action = 2)
heig_har<-heig_har[which(heig_har$mr_keep==TRUE),]


# wbfm
exposuredata<-data[,c("SNP", "wbfm_ea", "wbfm_nea", "wbfm_beta", "wbfm_se", "wbfm_p", "chr", "pos", "wbfm_eaf", 'wbfm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wbfm_beta",
                         se_col = "wbfm_se",
                         eaf_col = "wbfm_eaf",
                         effect_allele_col = "wbfm_ea",
                         other_allele_col = "wbfm_nea",
                         pval_col="wbfm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wbfm_ss")
wbfm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wbfm_har<-wbfm_har[which(wbfm_har$mr_keep==TRUE),]

# wblm
exposuredata<-data[,c("SNP", "wblm_ea", "wblm_nea", "wblm_beta", "wblm_se", "wblm_p", "chr", "pos", "wblm_eaf", 'wblm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wblm_beta",
                         se_col = "wblm_se",
                         eaf_col = "wblm_eaf",
                         effect_allele_col = "wblm_ea",
                         other_allele_col = "wblm_nea",
                         pval_col="wblm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wblm_ss")
wblm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wblm_har<-wblm_har[which(wblm_har$mr_keep==TRUE),]


# weigh
exposuredata<-data[,c("SNP", "weigh_ea", "weigh_nea", "weigh_beta", "weigh_se", "weigh_p", "chr", "pos", "weigh_eaf", 'weigh_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "weigh_beta",
                         se_col = "weigh_se",
                         eaf_col = "weigh_eaf",
                         effect_allele_col = "weigh_ea",
                         other_allele_col = "weigh_nea",
                         pval_col="weigh_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "weigh_ss")
weigh_har<-harmonise_data(exposuredat, outcomedat, action = 2)
weigh_har<-weigh_har[which(weigh_har$mr_keep==TRUE),]
rm(exposuredat,exposuredata,outcomedat,outcomedata)

# Clump
clp<-function(dat){
  dat <- get(dat, envir = .GlobalEnv)
  dat<-clump_data(dat, clump_r2=0.001, pop="EUR")
}
iv_list<- sapply(genelist, clp, simplify = FALSE)
invisible(lapply(names(iv_list), function(x) assign(x, iv_list[[x]],envir=.GlobalEnv)))

# Do MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  b<-0
  se<-0
  pval<-0
  n.SNP<-0
  
  if(nrow(dat)==1){
    dat$effect_allele_col<-"A"
    dat$other_allele_col<-"C"
    exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             pval_col="pval.exposure",
                             eaf_col = "eaf.exposure",
                             effect_allele_col = "effect_allele_col",
                             other_allele_col = "other_allele_col")
    
    
    outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            pval_col="pval.outcome",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele_col",
                            other_allele_col = "other_allele_col")
    
    harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
    output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
    b<-output[1,7]
    se<-output[1,8]
    pval<-output[1,9]
    n.SNP<-1
    
  }else{
    output<-TwoSampleMR::mr_ivw(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
    b<-unlist(output)[1]
    se<-unlist(output)[2]
    pval<-unlist(output)[3]
    n.SNP<-unlist(output)[4]
  }
  
  outputobj<-list(b,se,pval,n.SNP)
  return(unlist(outputobj))
}
mr_table2 <- t(data.frame((sapply(genelist,mrfunc2))))
# mr_table2$or<-exp(mr_table2$coef)
# mr_table2$lci<-exp((mr_table2$coef)-((mr_table2$se)*1.96))
# mr_table2$uci<-exp((mr_table2$coef)+((mr_table2$se)*1.96))
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/print_results2.csv')

# Do MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  b<-0
  se<-0
  pval<-0
  n.SNP<-0
  
  if(nrow(dat)==1){
    dat$effect_allele_col<-"A"
    dat$other_allele_col<-"C"
    exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             pval_col="pval.exposure",
                             eaf_col = "eaf.exposure",
                             effect_allele_col = "effect_allele_col",
                             other_allele_col = "other_allele_col")
    
    
    outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            pval_col="pval.outcome",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele_col",
                            other_allele_col = "other_allele_col")
    
    harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
    output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
    b<-output[1,7]
    se<-output[1,8]
    pval<-output[1,9]
    n.SNP<-1
    
  }else{
    output<-TwoSampleMR::mr_ivw_fe(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
    b<-unlist(output)[1]
    se<-unlist(output)[2]
    pval<-unlist(output)[3]
    n.SNP<-unlist(output)[4]
  }
  
  outputobj<-list(b,se,pval,n.SNP)
  return(unlist(outputobj))
}

mr_table2 <- t(data.frame((sapply(genelist,mrfunc2))))
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/print_results2_fixed.csv')

# Sentitivity analyses 
#wm 
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_weighted_median(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
mr_table2 <- as.data.frame((sapply(genelist,mrfunc2)))
mr_table2<-t(mr_table2)
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/print_results2_wm.csv')

#egger
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_egger_regression(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
a <- (sapply(genelist,mrfunc2))
write.csv(a, '~/Desktop/bmi-ecg/res/print_results2_egg.csv')

rm(list=ls())
#



#####  OUTCOME 2 = qrsinterval #### 
#load exposures -  all CSVs from IV folder 
setwd("~/Desktop/bmi-ecg/ivs")    #set wd to where i want to read IV sets from 
files = list.files(pattern="*.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
instruments<- merge(bmi_iv, whr_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, heig_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wbfm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wblm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, weigh_iv[,-1], by="SNP", all=TRUE)
rm(files, data_list)
setwd("~/Desktop/bmi-ecg") #set wd back to project
#
#
# Extract outcome for all IVs - qrsint
qrsint <- as.data.frame(fread(("out/qrsint-woj.tsv")))
qrsint$qrsint_ea<-toupper(qrsint$qrsint_ea)
qrsint$qrsint_nea<-toupper(qrsint$qrsint_nea)
colnames(qrsint)

# merge with all exposures
instruments_qrsint <- qrsint[which(paste(qrsint$SNP) %in% paste(instruments$SNP)),]
data <- merge(instruments, instruments_qrsint, by=c("SNP"), all=FALSE)
data <- data[!duplicated(data$SNP),]
head(data)

# Separate 'data' file into exposure and outcome sets, and format respective dataframes for MR
outcomedata<-data[,c("SNP", "qrsint_ea", "qrsint_nea", "qrsint_beta", "qrsint_se", "qrsint_p","qrsint_eaf", "qrsint_ss")]
outcomedat<-format_data(outcomedata,type="outcome",snp_col = "SNP",
                        beta_col = "qrsint_beta",
                        se_col = "qrsint_se",
                        eaf_col = "qrsint_eaf",
                        effect_allele_col = "qrsint_ea",
                        other_allele_col = "qrsint_nea",
                        pval_col = "qrsint_p", 
                        samplesize_col = "qrsint_ss")

# BMI
exposuredata<-data[,c("SNP", "bmi_ea", "bmi_nea", "bmi_beta", "bmi_se", "bmi_p", "bmi_chr", "bmi_pos", "bmi_eaf", 'bmi_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "bmi_beta",
                         se_col = "bmi_se",
                         eaf_col = "bmi_eaf",
                         effect_allele_col = "bmi_ea",
                         other_allele_col = "bmi_nea",pos_col	="pos",
                         chr_col="chr",pval_col="bmi_p", 
                         samplesize_col = "bmi_ss")
bmi_har<-harmonise_data(exposuredat, outcomedat, action = 2)
bmi_har<-bmi_har[which(bmi_har$mr_keep==TRUE),]

# WHR
exposuredata<-data[,c("SNP", "whr_ea", "whr_nea", "whr_beta", "whr_se", "whr_p", "whr_chr", "whr_pos", "whr_eaf", 'whr_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "whr_beta",
                         se_col = "whr_se",
                         eaf_col = "whr_eaf",
                         effect_allele_col = "whr_ea",
                         other_allele_col = "whr_nea",pos_col	="pos",
                         chr_col="chr",pval_col="whr_p", 
                         samplesize_col = "whr_ss")
whr_har<-harmonise_data(exposuredat, outcomedat, action = 2)
whr_har<-whr_har[which(whr_har$mr_keep==TRUE),]

# Heig
exposuredata<-data[,c("SNP", "heig_ea", "heig_nea", "heig_beta", "heig_se", "heig_p", "heig_chr", "heig_pos", "heig_eaf", 'heig_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "heig_beta",
                         se_col = "heig_se",
                         eaf_col = "heig_eaf",
                         effect_allele_col = "heig_ea",
                         other_allele_col = "heig_nea",pos_col	="pos",
                         chr_col="chr",pval_col="heig_p", 
                         samplesize_col = "heig_ss")
heig_har<-harmonise_data(exposuredat, outcomedat, action = 2)
heig_har<-heig_har[which(heig_har$mr_keep==TRUE),]

# wbfm
exposuredata<-data[,c("SNP", "wbfm_ea", "wbfm_nea", "wbfm_beta", "wbfm_se", "wbfm_p", "chr", "pos", "wbfm_eaf", 'wbfm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wbfm_beta",
                         se_col = "wbfm_se",
                         eaf_col = "wbfm_eaf",
                         effect_allele_col = "wbfm_ea",
                         other_allele_col = "wbfm_nea",
                         pval_col="wbfm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wbfm_ss")
wbfm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wbfm_har<-wbfm_har[which(wbfm_har$mr_keep==TRUE),]

# wblm
exposuredata<-data[,c("SNP", "wblm_ea", "wblm_nea", "wblm_beta", "wblm_se", "wblm_p", "chr", "pos", "wblm_eaf", 'wblm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wblm_beta",
                         se_col = "wblm_se",
                         eaf_col = "wblm_eaf",
                         effect_allele_col = "wblm_ea",
                         other_allele_col = "wblm_nea",
                         pval_col="wblm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wblm_ss")
wblm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wblm_har<-wblm_har[which(wblm_har$mr_keep==TRUE),]


# weigh
exposuredata<-data[,c("SNP", "weigh_ea", "weigh_nea", "weigh_beta", "weigh_se", "weigh_p", "chr", "pos", "weigh_eaf", 'weigh_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "weigh_beta",
                         se_col = "weigh_se",
                         eaf_col = "weigh_eaf",
                         effect_allele_col = "weigh_ea",
                         other_allele_col = "weigh_nea",
                         pval_col="weigh_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "weigh_ss")
weigh_har<-harmonise_data(exposuredat, outcomedat, action = 2)
weigh_har<-weigh_har[which(weigh_har$mr_keep==TRUE),]


rm(exposuredat,exposuredata,outcomedat,outcomedata)


# Delete rows with NAs frorm each IV list
genelist<- gsub("_iv","_har", genelist, fixed = TRUE) #make list with names of all harmonised sets

nas<-function(dat){
  dat <- get(dat, envir = .GlobalEnv)
  dat<-na.omit(dat)
}
iv_list<- sapply(genelist, nas, simplify = FALSE)
invisible(lapply(names(iv_list), function(x) assign(x,iv_list[[x]],envir=.GlobalEnv)))

# Clump
clp<-function(dat){
  dat <- get(dat, envir = .GlobalEnv)
  dat<-clump_data(dat, clump_r2=0.001, pop="EUR")
}
iv_list<- sapply(genelist, clp, simplify = FALSE)
invisible(lapply(names(iv_list), function(x) assign(x, iv_list[[x]],envir=.GlobalEnv)))

# Do MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  b<-0
  se<-0
  pval<-0
  n.SNP<-0
  
  if(nrow(dat)==1){
    dat$effect_allele_col<-"A"
    dat$other_allele_col<-"C"
    exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             pval_col="pval.exposure",
                             eaf_col = "eaf.exposure",
                             effect_allele_col = "effect_allele_col",
                             other_allele_col = "other_allele_col")
    
    
    outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            pval_col="pval.outcome",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele_col",
                            other_allele_col = "other_allele_col")
    
    harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
    output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
    b<-output[1,7]
    se<-output[1,8]
    pval<-output[1,9]
    n.SNP<-1
    
  }else{
    output<-TwoSampleMR::mr_ivw(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
    b<-unlist(output)[1]
    se<-unlist(output)[2]
    pval<-unlist(output)[3]
    n.SNP<-unlist(output)[4]
  }
  
  outputobj<-list(b,se,pval,n.SNP)
  return(unlist(outputobj))
}

mr_table2 <- t(data.frame((sapply(genelist,mrfunc2))))
# mr_table2$or<-exp(mr_table2$coef)
# mr_table2$lci<-exp((mr_table2$coef)-((mr_table2$se)*1.96))
# mr_table2$uci<-exp((mr_table2$coef)+((mr_table2$se)*1.96))
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/qrsint_results2.csv')

# Do MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  b<-0
  se<-0
  pval<-0
  n.SNP<-0
  
  if(nrow(dat)==1){
    dat$effect_allele_col<-"A"
    dat$other_allele_col<-"C"
    exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             pval_col="pval.exposure",
                             eaf_col = "eaf.exposure",
                             effect_allele_col = "effect_allele_col",
                             other_allele_col = "other_allele_col")
    
    
    outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            pval_col="pval.outcome",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele_col",
                            other_allele_col = "other_allele_col")
    
    harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
    output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
    b<-output[1,7]
    se<-output[1,8]
    pval<-output[1,9]
    n.SNP<-1
    
  }else{
    output<-TwoSampleMR::mr_ivw_fe(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
    b<-unlist(output)[1]
    se<-unlist(output)[2]
    pval<-unlist(output)[3]
    n.SNP<-unlist(output)[4]
  }
  
  outputobj<-list(b,se,pval,n.SNP)
  return(unlist(outputobj))
}

mr_table2 <- t(data.frame((sapply(genelist,mrfunc2))))
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/qrsint_results2_fixed.csv')

# Sentitivity analyses 
#wm 
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_weighted_median(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
mr_table2 <- as.data.frame((sapply(genelist,mrfunc2)))
mr_table2<-t(mr_table2)
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/qrsint_results2_wm.csv')

#egger
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_egger_regression(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
a <- (sapply(genelist,mrfunc2))
write.csv(a, '~/Desktop/bmi-ecg/res/qrsint_results2_egg.csv')

rm(list=ls())
# #


#####  OUTCOME 3 = pwd #### 
#load exposures -  all CSVs from IV folder 
setwd("~/Desktop/bmi-ecg/ivs")    #set wd to where i want to read IV sets from 
files = list.files(pattern="*.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
instruments<- merge(bmi_iv, whr_iv, by="SNP", all=TRUE)
instruments<- merge(instruments, heig_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wbfm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wblm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, weigh_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, lvm_iv, by="SNP", all=TRUE)
rm(files, data_list)
setwd("~/Desktop/bmi-ecg") #set wd back to project
#
#
# Extract outcome for all IVs - pwd
pwd <- as.data.frame(fread(("out/pwd-chris.tsv")))
pwd$pwd_ea<-toupper(pwd$pwd_ea)
pwd$pwd_nea<-toupper(pwd$pwd_nea)
colnames(pwd)

# merge with all exposures
instruments_pwd <- pwd[which(paste(pwd$SNP) %in% paste(instruments$SNP)),]
data <- merge(instruments, instruments_pwd, by=c("SNP"), all=FALSE)
data <- data[!duplicated(data$SNP),]
head(data)

# Separate 'data' file into exposure and outcome sets, and format respective dataframes for MR
outcomedata<-data[,c("SNP", "pwd_ea", "pwd_nea", "pwd_beta", "pwd_se", "pwd_p",  "bmi_eaf")]
outcomedat<-format_data(outcomedata,type="outcome",snp_col = "SNP",
                        beta_col = "pwd_beta",
                        se_col = "pwd_se",
                        effect_allele_col = "pwd_ea",
                        other_allele_col = "pwd_nea", pval_col = "pwd_p" )

# BMI
exposuredata<-data[,c("SNP", "bmi_ea", "bmi_nea", "bmi_beta", "bmi_se", "bmi_p", "bmi_chr", "bmi_pos", "bmi_eaf", 'bmi_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "bmi_beta",
                         se_col = "bmi_se",
                         eaf_col = "bmi_eaf",
                         effect_allele_col = "bmi_ea",
                         other_allele_col = "bmi_nea",pos_col	="pos",
                         chr_col="chr",pval_col="bmi_p", 
                         samplesize_col = "bmi_ss")
bmi_har<-harmonise_data(exposuredat, outcomedat, action = 2)
bmi_har<-bmi_har[which(bmi_har$mr_keep==TRUE),]

# WHR
exposuredata<-data[,c("SNP", "whr_ea", "whr_nea", "whr_beta", "whr_se", "whr_p", "whr_chr", "whr_pos", "whr_eaf", 'whr_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "whr_beta",
                         se_col = "whr_se",
                         eaf_col = "whr_eaf",
                         effect_allele_col = "whr_ea",
                         other_allele_col = "whr_nea",pos_col	="pos",
                         chr_col="chr",pval_col="whr_p", 
                         samplesize_col = "whr_ss")
whr_har<-harmonise_data(exposuredat, outcomedat, action = 2)
whr_har<-whr_har[which(whr_har$mr_keep==TRUE),]

# Heig
exposuredata<-data[,c("SNP", "heig_ea", "heig_nea", "heig_beta", "heig_se", "heig_p", "heig_chr", "heig_pos", "heig_eaf", 'heig_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "heig_beta",
                         se_col = "heig_se",
                         eaf_col = "heig_eaf",
                         effect_allele_col = "heig_ea",
                         other_allele_col = "heig_nea",pos_col	="pos",
                         chr_col="chr",pval_col="heig_p", 
                         samplesize_col = "heig_ss")
heig_har<-harmonise_data(exposuredat, outcomedat, action = 2)
heig_har<-heig_har[which(heig_har$mr_keep==TRUE),]

# wbfm
exposuredata<-data[,c("SNP", "wbfm_ea", "wbfm_nea", "wbfm_beta", "wbfm_se", "wbfm_p", "chr", "pos", "wbfm_eaf", 'wbfm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wbfm_beta",
                         se_col = "wbfm_se",
                         eaf_col = "wbfm_eaf",
                         effect_allele_col = "wbfm_ea",
                         other_allele_col = "wbfm_nea",
                         pval_col="wbfm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wbfm_ss")
wbfm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wbfm_har<-wbfm_har[which(wbfm_har$mr_keep==TRUE),]

# wblm
exposuredata<-data[,c("SNP", "wblm_ea", "wblm_nea", "wblm_beta", "wblm_se", "wblm_p", "chr", "pos", "wblm_eaf", 'wblm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wblm_beta",
                         se_col = "wblm_se",
                         eaf_col = "wblm_eaf",
                         effect_allele_col = "wblm_ea",
                         other_allele_col = "wblm_nea",
                         pval_col="wblm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wblm_ss")
wblm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wblm_har<-wblm_har[which(wblm_har$mr_keep==TRUE),]


# weigh
exposuredata<-data[,c("SNP", "weigh_ea", "weigh_nea", "weigh_beta", "weigh_se", "weigh_p", "chr", "pos", "weigh_eaf", 'weigh_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "weigh_beta",
                         se_col = "weigh_se",
                         eaf_col = "weigh_eaf",
                         effect_allele_col = "weigh_ea",
                         other_allele_col = "weigh_nea",
                         pval_col="weigh_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "weigh_ss")
weigh_har<-harmonise_data(exposuredat, outcomedat, action = 2)
weigh_har<-weigh_har[which(weigh_har$mr_keep==TRUE),]

rm(exposuredat,exposuredata,outcomedat,outcomedata)

genelist<- gsub("_iv","_har", genelist, fixed = TRUE) #make list with names of all harmonised sets


# Clump
clp<-function(dat){
  dat <- get(dat, envir = .GlobalEnv)
  dat<-clump_data(dat, clump_r2=0.001, pop="EUR")
}
iv_list<- sapply(genelist, clp, simplify = FALSE)
invisible(lapply(names(iv_list), function(x) assign(x, iv_list[[x]],envir=.GlobalEnv)))

# Do MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  b<-0
  se<-0
  pval<-0
  n.SNP<-0
  
  if(nrow(dat)==1){
    dat$effect_allele_col<-"A"
    dat$other_allele_col<-"C"
    exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             pval_col="pval.exposure",
                             eaf_col = "eaf.exposure",
                             effect_allele_col = "effect_allele_col",
                             other_allele_col = "other_allele_col")
    
    
    outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            pval_col="pval.outcome",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele_col",
                            other_allele_col = "other_allele_col")
    
    harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
    output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
    b<-output[1,7]
    se<-output[1,8]
    pval<-output[1,9]
    n.SNP<-1
    
  }else{
    output<-TwoSampleMR::mr_ivw(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
    b<-unlist(output)[1]
    se<-unlist(output)[2]
    pval<-unlist(output)[3]
    n.SNP<-unlist(output)[4]
  }
  
  outputobj<-list(b,se,pval,n.SNP)
  return(unlist(outputobj))
}

mr_table2 <- t(data.frame((sapply(genelist,mrfunc2))))
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/pwd_results2.csv')


# Do MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  b<-0
  se<-0
  pval<-0
  n.SNP<-0
  
  if(nrow(dat)==1){
    dat$effect_allele_col<-"A"
    dat$other_allele_col<-"C"
    exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             pval_col="pval.exposure",
                             eaf_col = "eaf.exposure",
                             effect_allele_col = "effect_allele_col",
                             other_allele_col = "other_allele_col")
    
    
    outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            pval_col="pval.outcome",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele_col",
                            other_allele_col = "other_allele_col")
    
    harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
    output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
    b<-output[1,7]
    se<-output[1,8]
    pval<-output[1,9]
    n.SNP<-1
    
  }else{
    output<-TwoSampleMR::mr_ivw_fe(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
    b<-unlist(output)[1]
    se<-unlist(output)[2]
    pval<-unlist(output)[3]
    n.SNP<-unlist(output)[4]
  }
  
  outputobj<-list(b,se,pval,n.SNP)
  return(unlist(outputobj))
}

mr_table2 <- t(data.frame((sapply(genelist,mrfunc2))))
# mr_table2$or<-exp(mr_table2$coef)
# mr_table2$lci<-exp((mr_table2$coef)-((mr_table2$se)*1.96))
# mr_table2$uci<-exp((mr_table2$coef)+((mr_table2$se)*1.96))
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/pwd_results2_fixed.csv')


# Sentitivity analyses 
#wm 
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_weighted_median(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
mr_table2 <- as.data.frame((sapply(genelist,mrfunc2)))
mr_table2<-t(mr_table2)
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/pwd_results2_wm.csv')

#egger
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_egger_regression(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
a <- (sapply(genelist,mrfunc2))
write.csv(a, '~/Desktop/bmi-ecg/res/pwd_results2_egg.csv')

rm(list=ls())
# #



#####  OUTCOME 4 = qtcinterval #### 
#load exposures -  all CSVs from IV folder 
setwd("~/Desktop/bmi-ecg/ivs")    #set wd to where i want to read IV sets from 
files = list.files(pattern="*.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
instruments<- merge(bmi_iv, whr_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, heig_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wbfm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wblm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, weigh_iv[,-1], by="SNP", all=TRUE)
rm(files, data_list)
setwd("~/Desktop/bmi-ecg") #set wd back to project
#
#
# Extract outcome for all IVs - qtcint
qtcint <- as.data.frame(fread(("out/qtcint-nuf.tsv")))
qtcint$qtcint_ea<-toupper(qtcint$qtcint_ea)
qtcint$qtcint_nea<-toupper(qtcint$qtcint_nea)
colnames(qtcint)

# merge with all exposures
instruments_qtcint <- qtcint[which(paste(qtcint$SNP) %in% paste(instruments$SNP)),]
data <- merge(instruments, instruments_qtcint, by=c("SNP"), all=FALSE)
data <- data[!duplicated(data$SNP),]
head(data)

# Separate 'data' file into exposure and outcome sets, and format respective dataframes for MR
outcomedata<-data[,c("SNP", "qtcint_ea", "qtcint_nea", "qtcint_beta", "qtcint_se", "qtcint_p", "chr", "pos")]
outcomedat<-format_data(outcomedata,type="outcome",snp_col = "SNP",
                        beta_col = "qtcint_beta",
                        se_col = "qtcint_se",
                        effect_allele_col = "qtcint_ea",
                        other_allele_col = "qtcint_nea",
                        pos_col	="pos",
                        chr_col="chr",pval_col = "qtcint_p")

# BMI
exposuredata<-data[,c("SNP", "bmi_ea", "bmi_nea", "bmi_beta", "bmi_se", "bmi_p", "chr", "pos", "bmi_eaf", 'bmi_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "bmi_beta",
                         se_col = "bmi_se",
                         eaf_col = "bmi_eaf",
                         effect_allele_col = "bmi_ea",
                         other_allele_col = "bmi_nea",pos_col	="pos",
                         chr_col="chr",pval_col="bmi_p", 
                         samplesize_col = "bmi_ss")
bmi_har<-harmonise_data(exposuredat, outcomedat, action = 2)
bmi_har<-bmi_har[which(bmi_har$mr_keep==TRUE),]

# WHR
exposuredata<-data[,c("SNP", "whr_ea", "whr_nea", "whr_beta", "whr_se", "whr_p", "chr", "pos", "whr_eaf", 'whr_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "whr_beta",
                         se_col = "whr_se",
                         eaf_col = "whr_eaf",
                         effect_allele_col = "whr_ea",
                         other_allele_col = "whr_nea",pos_col	="pos",
                         chr_col="chr",pval_col="whr_p", 
                         samplesize_col = "whr_ss")
whr_har<-harmonise_data(exposuredat, outcomedat, action = 2)
whr_har<-whr_har[which(whr_har$mr_keep==TRUE),]

# Heig
exposuredata<-data[,c("SNP", "heig_ea", "heig_nea", "heig_beta", "heig_se", "heig_p", "chr", "pos", "heig_eaf", 'heig_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "heig_beta",
                         se_col = "heig_se",
                         eaf_col = "heig_eaf",
                         effect_allele_col = "heig_ea",
                         other_allele_col = "heig_nea",pos_col	="pos",
                         chr_col="chr",pval_col="heig_p", 
                         samplesize_col = "heig_ss")
heig_har<-harmonise_data(exposuredat, outcomedat, action = 2)
heig_har<-heig_har[which(heig_har$mr_keep==TRUE),]


# wbfm
exposuredata<-data[,c("SNP", "wbfm_ea", "wbfm_nea", "wbfm_beta", "wbfm_se", "wbfm_p", "chr", "pos", "wbfm_eaf", 'wbfm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wbfm_beta",
                         se_col = "wbfm_se",
                         eaf_col = "wbfm_eaf",
                         effect_allele_col = "wbfm_ea",
                         other_allele_col = "wbfm_nea",
                         pval_col="wbfm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wbfm_ss")
wbfm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wbfm_har<-wbfm_har[which(wbfm_har$mr_keep==TRUE),]

# wblm
exposuredata<-data[,c("SNP", "wblm_ea", "wblm_nea", "wblm_beta", "wblm_se", "wblm_p", "chr", "pos", "wblm_eaf", 'wblm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wblm_beta",
                         se_col = "wblm_se",
                         eaf_col = "wblm_eaf",
                         effect_allele_col = "wblm_ea",
                         other_allele_col = "wblm_nea",
                         pval_col="wblm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wblm_ss")
wblm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wblm_har<-wblm_har[which(wblm_har$mr_keep==TRUE),]


# weigh
exposuredata<-data[,c("SNP", "weigh_ea", "weigh_nea", "weigh_beta", "weigh_se", "weigh_p", "chr", "pos", "weigh_eaf", 'weigh_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "weigh_beta",
                         se_col = "weigh_se",
                         eaf_col = "weigh_eaf",
                         effect_allele_col = "weigh_ea",
                         other_allele_col = "weigh_nea",
                         pval_col="weigh_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "weigh_ss")
weigh_har<-harmonise_data(exposuredat, outcomedat, action = 2)
weigh_har<-weigh_har[which(weigh_har$mr_keep==TRUE),]


rm(exposuredat,exposuredata,outcomedat,outcomedata)


# Clump

genelist<- gsub("_iv","_har", genelist, fixed = TRUE) #make list with names of all harmonised sets


clp<-function(dat){
  dat <- get(dat, envir = .GlobalEnv)
  dat<-clump_data(dat, clump_r2=0.001, pop="EUR")
}
iv_list<- sapply(genelist, clp, simplify = FALSE)
invisible(lapply(names(iv_list), function(x) assign(x, iv_list[[x]],envir=.GlobalEnv)))

# Do MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  b<-0
  se<-0
  pval<-0
  n.SNP<-0
  
  if(nrow(dat)==1){
    dat$effect_allele_col<-"A"
    dat$other_allele_col<-"C"
    exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             pval_col="pval.exposure",
                             eaf_col = "eaf.exposure",
                             effect_allele_col = "effect_allele_col",
                             other_allele_col = "other_allele_col")
    
    
    outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            pval_col="pval.outcome",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele_col",
                            other_allele_col = "other_allele_col")
    
    harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
    output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
    b<-output[1,7]
    se<-output[1,8]
    pval<-output[1,9]
    n.SNP<-1
    
  }else{
    output<-TwoSampleMR::mr_ivw(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
    b<-unlist(output)[1]
    se<-unlist(output)[2]
    pval<-unlist(output)[3]
    n.SNP<-unlist(output)[4]
  }
  
  outputobj<-list(b,se,pval,n.SNP)
  return(unlist(outputobj))
}

mr_table2 <- t(data.frame((sapply(genelist,mrfunc2))))
# mr_table2$or<-exp(mr_table2$coef)
# mr_table2$lci<-exp((mr_table2$coef)-((mr_table2$se)*1.96))
# mr_table2$uci<-exp((mr_table2$coef)+((mr_table2$se)*1.96))
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/qtcint_results2.csv')


# Do MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  b<-0
  se<-0
  pval<-0
  n.SNP<-0
  
  if(nrow(dat)==1){
    dat$effect_allele_col<-"A"
    dat$other_allele_col<-"C"
    exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             pval_col="pval.exposure",
                             eaf_col = "eaf.exposure",
                             effect_allele_col = "effect_allele_col",
                             other_allele_col = "other_allele_col")
    
    
    outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            pval_col="pval.outcome",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele_col",
                            other_allele_col = "other_allele_col")
    
    harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
    output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
    b<-output[1,7]
    se<-output[1,8]
    pval<-output[1,9]
    n.SNP<-1
    
  }else{
    output<-TwoSampleMR::mr_ivw_fe(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
    b<-unlist(output)[1]
    se<-unlist(output)[2]
    pval<-unlist(output)[3]
    n.SNP<-unlist(output)[4]
  }
  
  outputobj<-list(b,se,pval,n.SNP)
  return(unlist(outputobj))
}

mr_table2 <- t(data.frame((sapply(genelist,mrfunc2))))
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/qtcint_results2_fixed.csv')

# Sentitivity analyses 
#wm 
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_weighted_median(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
mr_table2 <- as.data.frame((sapply(genelist,mrfunc2)))
mr_table2<-t(mr_table2)
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/qtcint_results2_wm.csv')

#egger
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_egger_regression(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
a <- (sapply(genelist,mrfunc2))
write.csv(a, '~/Desktop/bmi-ecg/res/qtcint_results2_egg.csv')

rm(list=ls())
#


#### MVMR PWD ####
# Load exposures and mediators
bmi <- as.data.frame(fread("/Volumes/MADDY2/mrdata/bmi_mvmrformat.tsv"))
heig <- as.data.frame(fread("/Volumes/MADDY2/mrdata/heig_mvmrformat.tsv"))
wblm <- as.data.frame(fread("/Volumes/MADDY2/mrdata/wblm_mvmrformat.tsv"))
wbfm <- as.data.frame(fread("/Volumes/MADDY2/mrdata/wbfm_mvmrformat.tsv"))
whr <- as.data.frame(fread("/Volumes/MADDY2/mrdata/whr_mvmrformat.tsv"))
weigh <- as.data.frame(fread("/Volumes/MADDY2/mrdata/weigh_mvmrformat.tsv"))

# load outcome
pwd<-as.data.frame(fread("/Volumes/MADDY2/mrdata/pwd-chris.tsv"))


bmi$bmi_p<-as.numeric(bmi$bmi_p)
whr$whr_p<-as.numeric(whr$whr_p)
wblm$wblm_p<-as.numeric(wblm$wblm_p)
wbfm$wbfm_p<-as.numeric(wbfm$wbfm_p)
heig$heig_p<-as.numeric(heig$heig_p)
weigh$weigh_p<-as.numeric(weigh$weigh_p)

#Select instruments
bmi_i <- subset(bmi, bmi_p<5e-8)
heig_i <- subset(heig, heig_p<5e-8)
weigh_i <- subset(weigh, weigh_p<5e-8)
whr_i <- subset(whr, whr_p<5e-8)
wblm_i <- subset(wblm, wblm_p<5e-8)
wbfm_i <- subset(wbfm, wbfm_p<5e-8)

##Merge instruments
instruments <- merge(bmi_i,  whr_i, by="SNP", all=TRUE)
instruments <- merge(instruments, heig_i, by="SNP", all=TRUE)
instruments <- merge(instruments, wblm_i, by="SNP", all=TRUE)
instruments <- merge(instruments, wbfm_i, by="SNP", all=TRUE)
instruments <- merge(instruments, weigh_i, by="SNP", all=TRUE)
dim(instruments)
instruments <- instruments[!duplicated(instruments$SNP),]
dim(instruments)

#Merge association estimates
instruments_whr <- whr[which(paste(whr$SNP) %in% paste(instruments$SNP)),]
instruments_bmi <- bmi[which(paste(bmi$SNP) %in% paste(instruments$SNP)),]
instruments_weigh <- weigh[which(paste(weigh$SNP) %in% paste(instruments$SNP)),]
instruments_heig <- heig[which(paste(heig$SNP) %in% paste(instruments$SNP)),]
instruments_wblm <- wblm[which(paste(wblm$SNP) %in% paste(instruments$SNP)),]
instruments_wbfm <- wbfm[which(paste(wbfm$SNP) %in% paste(instruments$SNP)),]
instruments_pwd <- pwd[which(paste(pwd$SNP) %in% paste(instruments$SNP)),]

data <- merge(instruments_bmi, instruments_whr, by="SNP", all=TRUE)
data <- merge(data, instruments_heig, by="SNP", all=TRUE)
data <- merge(data, instruments_weigh, by="SNP", all=TRUE)
data <- merge(data, instruments_wbfm, by="SNP", all=TRUE)
data <- merge(data, instruments_wblm, by="SNP", all=TRUE)
data <- merge(data, instruments_pwd, by="SNP", all=TRUE)
dim(data)
data <- data[!duplicated(data$SNP),]
dim(data)
head(data)

databmi<-data[,c("SNP","bmi_ea","bmi_nea","bmi_eaf","bmi_beta","bmi_se","bmi_p")]
datawhr<-data[,c("SNP","whr_ea","whr_nea","whr_eaf","whr_beta","whr_se","whr_p")]
datawblm<-data[,c("SNP","wblm_ea","wblm_nea","wblm_eaf","wblm_beta","wblm_se","wblm_p")]
datawbfm<-data[,c("SNP","wbfm_ea","wbfm_nea","wbfm_eaf","wbfm_beta","wbfm_se","wbfm_p")]
dataweigh<-data[,c("SNP","weigh_ea","weigh_nea","weigh_eaf","weigh_beta","weigh_se","weigh_p")]
dataheig<-data[,c("SNP","heig_ea","heig_nea","heig_eaf","heig_beta","heig_se","heig_p")]
datapwd<-data[,c("SNP","pwd_ea","pwd_nea","pwd_beta","pwd_se","pwd_p")]
databmi$phenotype="bmi"
datawbfm$phenotype="wbfm"
datawblm$phenotype="wblm"
dataheig$phenotype="heig"
dataweigh$phenotype="weigh"
datawhr$phenotype="whr"
datapwd$phenotype="pwd"


colnames(databmi)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(datawbfm)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(datawblm)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(dataheig)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(dataweigh)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(datawhr)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(datapwd)<-c("SNP","effect_allele","other_allele","beta","se","pval","phenotype")

setwd("~/Desktop/bmi-ecg/mvmr/tempdata/pwd")
write.csv(databmi,  "databmi_beforeharmonizing.csv",row.names = FALSE)
write.csv(datawbfm,  "datawbfm_beforeharmonizing.csv",row.names = FALSE)
write.csv(datawblm,  "datawblm_beforeharmonizing.csv",row.names = FALSE)
write.csv(dataheig,  "dataheig_beforeharmonizing.csv",row.names = FALSE)
write.csv(dataweigh,  "dataweigh_beforeharmonizing.csv",row.names = FALSE)
write.csv(datawhr,  "datawhr_beforeharmonizing.csv",row.names = FALSE)
write.csv(datapwd,  "datapwd_beforeharmonizing.csv",row.names = FALSE)

rm(list=ls())

#MVMR
setwd("~/Desktop/bmi-ecg/mvmr/tempdata/pwd")
wb <- createWorkbook()

# bmi-wbfm on pwd
exposuredat=mv_extract_exposures_local(filenames_exposure=c("databmi_beforeharmonizing.csv","datawbfm_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("datapwd_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata<-mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results<-mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "bmi-wbfm on pwd")
writeData(wb, "bmi-wbfm on pwd", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres1 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres1
pres1 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres1


# bmi-wblm on pwd
exposuredat = mv_extract_exposures_local(filenames_exposure=c("databmi_beforeharmonizing.csv","datawblm_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("datapwd_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata <- mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results <- mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "bmi-wblm on pwd")
writeData(wb, "bmi-wblm on pwd", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres2 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres2

pres2 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres2

fst <- rbind(sres1, sres2)

# bmi-heig on pwd
exposuredat = mv_extract_exposures_local(filenames_exposure=c("databmi_beforeharmonizing.csv","dataheig_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("datapwd_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata <- mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results <- mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "bmi-heig on pwd")
writeData(wb, "bmi-heig on pwd", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres2 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres2

pres2 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres2

fst <- rbind(fst, sres2)

# wbfm-wblm on pwd
exposuredat = mv_extract_exposures_local(filenames_exposure=c("datawbfm_beforeharmonizing.csv","datawblm_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("datapwd_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata <- mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results <- mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "wbfm-wblm on pwd")
writeData(wb, "wbfm-wblm on pwd", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres2 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres2

pres2 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres2

fst <- rbind(fst, sres2)

# wblm-heig on pwd
exposuredat = mv_extract_exposures_local(filenames_exposure=c("datawblm_beforeharmonizing.csv","dataheig_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("datapwd_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata <- mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results <- mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "wblm-heig on pwd")
writeData(wb, "wblm-heig on pwd", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres2 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres2

pres2 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres2

fst <- rbind(fst, sres2)


colnames(fst)<- c("exposure1", "exposure2")
rownames(fst)<- c("bmixwbfm",  "bmixwblm", "bmixheig",  "wblmxheig")
write.csv(fst, "~/Desktop/bmi-ecg/mvmr/res/fstatpwd.csv")
saveWorkbook(wb, file = "~/Desktop/bmi-ecg/mvmr/res/mediation_pwd.xlsx", overwrite = TRUE)


# whr-wbfm-wblm-heig on pwd
setwd("~/Desktop/bmi-ecg/mvmr/tempdata/pwd")
wb <- createWorkbook()
exposuredat=mv_extract_exposures_local(filenames_exposure=c("datawhr_beforeharmonizing.csv","datawbfm_beforeharmonizing.csv",
                                                            "datawblm_beforeharmonizing.csv","dataheig_beforeharmonizing.csv"),
                                       sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("datapwd_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata<-mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results<-mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "whr-wbfm-wblm-heig on pwd")
writeData(wb, "whr-wbfm-wblm-heig on pwd", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres1 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres1
colnames(sres1)<- c("whr", "wbfm", 'wblm', 'heig')
write.csv(sres1, "~/Desktop/bmi-ecg/mvmr/res/fstatpwd_all3.csv")
saveWorkbook(wb, file = "~/Desktop/bmi-ecg/mvmr/res/mediation_all3_pwd.xlsx", overwrite = TRUE)


# whr-wbfm-wblm-heig on pwd
setwd("~/Desktop/bmi-ecg/mvmr/tempdata/pwd")
wb <- createWorkbook()
exposuredat=mv_extract_exposures_local(filenames_exposure=c("datawhr_beforeharmonizing.csv","datawbfm_beforeharmonizing.csv",
                                                            "datawblm_beforeharmonizing.csv","dataheig_beforeharmonizing.csv",
                                                            "databmi_beforeharmonizing.csv"),
                                       sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("datapwd_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata<-mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results<-mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "whr-wbfm-wblm-heig-bmi on pwd")
writeData(wb, "whr-wbfm-wblm-heig-bmi on pwd", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres1 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres1
colnames(sres1)<- c("whr", "wbfm", 'wblm', 'heig', 'bmi')
write.csv(sres1, "~/Desktop/bmi-ecg/mvmr/res/fstatpwd_all2.csv")
saveWorkbook(wb, file = "~/Desktop/bmi-ecg/mvmr/res/mediation_all2_pwd.xlsx", overwrite = TRUE)

rm(list=ls())



#### MVMR QTC ####
# Load exposures and mediators
bmi <- as.data.frame(fread("/Volumes/MADDY2/mrdata/bmi_mvmrformat.tsv"))
heig <- as.data.frame(fread("/Volumes/MADDY2/mrdata/heig_mvmrformat.tsv"))
wblm <- as.data.frame(fread("/Volumes/MADDY2/mrdata/wblm_mvmrformat.tsv"))
wbfm <- as.data.frame(fread("/Volumes/MADDY2/mrdata/wbfm_mvmrformat.tsv"))
whr <- as.data.frame(fread("/Volumes/MADDY2/mrdata/whr_mvmrformat.tsv"))
weigh <- as.data.frame(fread("/Volumes/MADDY2/mrdata/weigh_mvmrformat.tsv"))

# load outcome
qtc<-as.data.frame(fread("/Volumes/MADDY2/mrdata/qtcint-nuf.tsv"))
head(qtc)
qtc <- qtc[complete.cases(qtc$SNP), ] 


bmi$bmi_p<-as.numeric(bmi$bmi_p)
whr$whr_p<-as.numeric(whr$whr_p)
wblm$wblm_p<-as.numeric(wblm$wblm_p)
wbfm$wbfm_p<-as.numeric(wbfm$wbfm_p)
heig$heig_p<-as.numeric(heig$heig_p)
weigh$weigh_p<-as.numeric(weigh$weigh_p)

#Select instruments
bmi_i <- subset(bmi, bmi_p<5e-8)
heig_i <- subset(heig, heig_p<5e-8)
weigh_i <- subset(weigh, weigh_p<5e-8)
whr_i <- subset(whr, whr_p<5e-8)
wblm_i <- subset(wblm, wblm_p<5e-8)
wbfm_i <- subset(wbfm, wbfm_p<5e-8)

##Merge instruments
instruments <- merge(bmi_i,  whr_i, by="SNP", all=TRUE)
instruments <- merge(instruments, heig_i, by="SNP", all=TRUE)
instruments <- merge(instruments, wblm_i, by="SNP", all=TRUE)
instruments <- merge(instruments, wbfm_i, by="SNP", all=TRUE)
instruments <- merge(instruments, weigh_i, by="SNP", all=TRUE)
dim(instruments)
instruments <- instruments[!duplicated(instruments$SNP),]
dim(instruments)

#Merge association estimates
instruments_whr <- whr[which(paste(whr$SNP) %in% paste(instruments$SNP)),]
instruments_bmi <- bmi[which(paste(bmi$SNP) %in% paste(instruments$SNP)),]
instruments_weigh <- weigh[which(paste(weigh$SNP) %in% paste(instruments$SNP)),]
instruments_heig <- heig[which(paste(heig$SNP) %in% paste(instruments$SNP)),]
instruments_wblm <- wblm[which(paste(wblm$SNP) %in% paste(instruments$SNP)),]
instruments_wbfm <- wbfm[which(paste(wbfm$SNP) %in% paste(instruments$SNP)),]
instruments_qtc <- qtc[which(paste(qtc$SNP) %in% paste(instruments$SNP)),]

data <- merge(instruments_bmi, instruments_whr, by="SNP", all=TRUE)
data <- merge(data, instruments_heig, by="SNP", all=TRUE)
data <- merge(data, instruments_weigh, by="SNP", all=TRUE)
data <- merge(data, instruments_wbfm, by="SNP", all=TRUE)
data <- merge(data, instruments_wblm, by="SNP", all=TRUE)
data <- merge(data, instruments_qtc, by="SNP", all=FALSE)
dim(data)
data <- data[!duplicated(data$SNP),]
dim(data)
head(data)
rm(whr, bmi, wblm, wbfm, weigh, heig)

databmi<-data[,c("SNP","bmi_ea","bmi_nea","bmi_eaf","bmi_beta","bmi_se","bmi_p")]
datawhr<-data[,c("SNP","whr_ea","whr_nea","whr_eaf","whr_beta","whr_se","whr_p")]
datawblm<-data[,c("SNP","wblm_ea","wblm_nea","wblm_eaf","wblm_beta","wblm_se","wblm_p")]
datawbfm<-data[,c("SNP","wbfm_ea","wbfm_nea","wbfm_eaf","wbfm_beta","wbfm_se","wbfm_p")]
dataweigh<-data[,c("SNP","weigh_ea","weigh_nea","weigh_eaf","weigh_beta","weigh_se","weigh_p")]
dataheig<-data[,c("SNP","heig_ea","heig_nea","heig_eaf","heig_beta","heig_se","heig_p")]
dataqtc<-data[,c("SNP","qtcint_ea","qtcint_nea","qtcint_beta","qtcint_se","qtcint_p")]
databmi$phenotype="bmi"
datawbfm$phenotype="wbfm"
datawblm$phenotype="wblm"
dataheig$phenotype="heig"
dataweigh$phenotype="weigh"
datawhr$phenotype="whr"
dataqtc$phenotype="qtc"


colnames(databmi)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(datawbfm)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(datawblm)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(dataheig)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(dataweigh)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(datawhr)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval","phenotype")
colnames(dataqtc)<-c("SNP","effect_allele","other_allele","beta","se","pval","phenotype")

setwd("~/Desktop/bmi-ecg/mvmr/tempdata/qtc")
write.csv(databmi,  "databmi_beforeharmonizing.csv",row.names = FALSE)
write.csv(datawbfm,  "datawbfm_beforeharmonizing.csv",row.names = FALSE)
write.csv(datawblm,  "datawblm_beforeharmonizing.csv",row.names = FALSE)
write.csv(dataheig,  "dataheig_beforeharmonizing.csv",row.names = FALSE)
write.csv(dataweigh,  "dataweigh_beforeharmonizing.csv",row.names = FALSE)
write.csv(datawhr,  "datawhr_beforeharmonizing.csv",row.names = FALSE)
write.csv(dataqtc,  "dataqtc_beforeharmonizing.csv",row.names = FALSE)

rm(list=ls())

#MVMR
setwd("~/Desktop/bmi-ecg/mvmr/tempdata/qtc")
library(openxlsx)
wb <- createWorkbook()
library(MVMR)

# bmi-wbfm on qtc
exposuredat=mv_extract_exposures_local(filenames_exposure=c("databmi_beforeharmonizing.csv","datawbfm_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("dataqtc_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata<-mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results<-mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "bmi-wbfm on qtc")
writeData(wb, "bmi-wbfm on qtc", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres1 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres1
pres1 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres1


# bmi-wblm on qtc
exposuredat = mv_extract_exposures_local(filenames_exposure=c("databmi_beforeharmonizing.csv","datawblm_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("dataqtc_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata <- mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results <- mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "bmi-wblm on qtc")
writeData(wb, "bmi-wblm on qtc", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres2 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres2

pres2 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres2

fst <- rbind(sres1, sres2)

# bmi-heig on qtc
exposuredat = mv_extract_exposures_local(filenames_exposure=c("databmi_beforeharmonizing.csv","dataheig_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("dataqtc_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata <- mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results <- mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "bmi-heig on qtc")
writeData(wb, "bmi-heig on qtc", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres2 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres2

pres2 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres2

fst <- rbind(fst, sres2)

# wblm-wbfm on qtc
exposuredat = mv_extract_exposures_local(filenames_exposure=c("datawblm_beforeharmonizing.csv","datawbfm_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("dataqtc_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata <- mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results <- mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "wblm-wbfm on qtc")
writeData(wb, "wblm-wbfm on qtc", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres2 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres2

pres2 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres2

fst <- rbind(fst, sres2)

# wblm-heig on qtc
exposuredat = mv_extract_exposures_local(filenames_exposure=c("datawblm_beforeharmonizing.csv","dataheig_beforeharmonizing.csv"),sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("dataqtc_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata <- mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results <- mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "wblm-heig on qtc")
writeData(wb, "wblm-heig on qtc", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres2 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres2

pres2 <- pleiotropy_mvmr(r_input = F.stat, gencov = 0)
pres2

fst <- rbind(fst, sres2)


colnames(fst)<- c("exposure1", "exposure2")
rownames(fst)<- c("bmixwbfm",  "bmixwblm", "bmixheig",  "wblmxheig")
write.csv(fst, "~/Desktop/bmi-ecg/mvmr/res/fstatqtc.csv")

saveWorkbook(wb, file = "~/Desktop/bmi-ecg/mvmr/res/mediation_qtc.xlsx", overwrite = TRUE)


# bmi-wbfm-wblm-heig on qtc
setwd("~/Desktop/bmi-ecg/mvmr/tempdata/qtc")
wb <- createWorkbook()
exposuredat=mv_extract_exposures_local(filenames_exposure=c("databmi_beforeharmonizing.csv","datawbfm_beforeharmonizing.csv",
                                                            "datawblm_beforeharmonizing.csv","dataheig_beforeharmonizing.csv"),
                                       sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("dataqtc_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata<-mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results<-mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "bmi-wbfm-wblm-heig on qtc")
writeData(wb, "bmi-wbfm-wblm-heig on qtc", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres1 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres1
colnames(sres1)<- c("bmi", "wbfm", 'wblm', 'heig')
write.csv(sres1, "~/Desktop/bmi-ecg/mvmr/res/fstatqtc_all.csv")
saveWorkbook(wb, file = "~/Desktop/bmi-ecg/mvmr/res/mediation_all_qtc.xlsx", overwrite = TRUE)


# whr-wbfm-wblm-heig on qtc
setwd("~/Desktop/bmi-ecg/mvmr/tempdata/qtc")
wb <- createWorkbook()
exposuredat=mv_extract_exposures_local(filenames_exposure=c("datawhr_beforeharmonizing.csv","datawbfm_beforeharmonizing.csv",
                                                            "datawblm_beforeharmonizing.csv","dataheig_beforeharmonizing.csv"),
                                       sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("dataqtc_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata<-mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results<-mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "whr-wbfm-wblm-heig on qtc")
writeData(wb, "whr-wbfm-wblm-heig on qtc", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres1 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres1
colnames(sres1)<- c("whr", "wbfm", 'wblm', 'heig')
write.csv(sres1, "~/Desktop/bmi-ecg/mvmr/res/fstatqtc_all3.csv")
saveWorkbook(wb, file = "~/Desktop/bmi-ecg/mvmr/res/mediation_all3_qtc.xlsx", overwrite = TRUE)


# whr-wbfm-wblm-heig-bmi on qtc
setwd("~/Desktop/bmi-ecg/mvmr/tempdata/qtc")
wb <- createWorkbook()
exposuredat=mv_extract_exposures_local(filenames_exposure=c("datawhr_beforeharmonizing.csv","datawbfm_beforeharmonizing.csv",
                                                            "datawblm_beforeharmonizing.csv","dataheig_beforeharmonizing.csv",
                                                            "databmi_beforeharmonizing.csv"),
                                       sep = ",",phenotype_col="phenotype")
outcomedat = read_outcome_data("dataqtc_beforeharmonizing.csv",sep=",",phenotype_col = "phenotype")
mvdata<-mv_harmonise_data(exposuredat, outcomedat, harmonise_strictness = 2)
results<-mv_multiple(mvdata, pval_threshold = 5e-08)
addWorksheet(wb, "whr-wbfm-wblm-heig-bmi on qtc")
writeData(wb, "whr-wbfm-wblm-heig-bmi on qtc", results$result, startRow = 1, startCol = 1)
# conditional f-stats
F.stat <- format_mvmr(BXGs = mvdata[["exposure_beta"]], #beta of exp and mediator
                      BYG = mvdata[["outcome_beta"]],#beta outcome
                      seBXGs = mvdata[["exposure_se"]], 
                      seBYG = mvdata[["outcome_se"]],
                      RSID = rownames(mvdata[["exposure_beta"]]))
sres1 <- strength_mvmr(F.stat,  0) #as non-overlapping samples for exposure, therefore cov=0
sres1
colnames(sres1)<- c("whr", "wbfm", 'wblm', 'heig', 'bmi',)
write.csv(sres1, "~/Desktop/bmi-ecg/mvmr/res/fstatqtc_all2.csv")
saveWorkbook(wb, file = "~/Desktop/bmi-ecg/mvmr/res/mediation_all2_qtc.xlsx", overwrite = TRUE)

rm(list=ls())

rm(list=ls())




#### Phenoscanner ####
library(phenoscanner)
setwd("~/Desktop/bmi-ecg/ivs")    #set wd to where i want to read IV sets from 

# BMI 
bmi_iv <- as.data.frame(fread('bmi_iv.csv'))
bmi_iv<-clump_data(bmi_iv, clump_r2=0.001, pop="EUR")

# phenoscan for each set
bmisnps <- as.character(bmi_iv$SNP)
# write.csv(bmisnps, "~/Desktop/bmi-ecg/phenoscanner/snp-bmi.csv")
psbmi <- phenoscanner(snpquery=bmisnps[1:99], pvalue = 5e-8)
psbmi1 <- phenoscanner(snpquery=bmisnps[100:199], pvalue = 5e-8)
psbmi2 <- phenoscanner(snpquery=bmisnps[200:299], pvalue = 5e-8)
psbmi3 <- phenoscanner(snpquery=bmisnps[300:399], pvalue = 5e-8)
psbmi4 <- phenoscanner(snpquery=bmisnps[400:499], pvalue = 5e-8)
psbmi5 <- phenoscanner(snpquery=bmisnps[500:516], pvalue = 5e-8)
psbmi <- rbind(psbmi$results, psbmi1$results)
psbmi <- rbind(psbmi, psbmi2$results)
psbmi <- rbind(psbmi, psbmi3$results)
psbmi <- rbind(psbmi, psbmi4$results)
psbmi <- rbind(psbmi, psbmi5$results)

write.csv(psbmi, "~/Desktop/birthw-cmr/phenoscanner/psbmi.csv")

#####  OUTCOME 5 = pwtforce #### 
#load exposures -  all CSVs from IV folder 
setwd("~/Desktop/bmi-ecg/ivs")    #set wd to where i want to read IV sets from 
files = list.files(pattern="*.csv")   #make list of all csv names
data_list = lapply(files, read.csv, header = TRUE)    #apply reead csv to all names
names(data_list) <- gsub(".csv","",
                         list.files(full.names = FALSE),
                         fixed = TRUE)   #make names for GE - withoout the .csv
invisible(lapply(names(data_list), function(x) assign(x,data_list[[x]],envir=.GlobalEnv)))  #pull imported data out from df list and assign names
genelist<-ls(Filter(function(x) is(x, "data.frame"), mget(ls()))) #make list with names of all IV sets
instruments<- merge(bmi_iv, whr_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, heig_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wbfm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, wblm_iv[,-1], by="SNP", all=TRUE)
instruments<- merge(instruments, weigh_iv[,-1], by="SNP", all=TRUE)
rm(files, data_list)
setwd("~/Desktop/bmi-ecg") #set wd back to project
#
#
# Extract outcome for all IVs - pwtforce
pwtforce <- as.data.frame(fread(("out/pwtforce-chris.tsv")))
pwtforce$pwtforce_ea<-toupper(pwtforce$pwtforce_ea)
pwtforce$pwtforce_nea<-toupper(pwtforce$pwtforce_nea)
colnames(pwtforce)

# merge with all exposures
instruments_pwtforce <- pwtforce[which(paste(pwtforce$SNP) %in% paste(instruments$SNP)),]
data <- merge(instruments, instruments_pwtforce, by=c("SNP"), all=FALSE)
data <- data[!duplicated(data$SNP),]
head(data)

# Separate 'data' file into exposure and outcome sets, and format respective dataframes for MR
outcomedata<-data[,c("SNP", "pwtforce_ea", "pwtforce_nea", "pwtforce_beta", "pwtforce_se", "pwtforce_p", "bmi_chr", "bmi_pos", "bmi_eaf")]
outcomedat<-format_data(outcomedata,type="outcome",snp_col = "SNP",
                        beta_col = "pwtforce_beta",
                        se_col = "pwtforce_se",
                        eaf_col = "bmi_eaf",
                        effect_allele_col = "pwtforce_ea",
                        other_allele_col = "pwtforce_nea",
                        pos_col	="bmi_pos",
                        chr_col="bmi_chr",pval_col = "pwtforce_p" )

# BMI
exposuredata<-data[,c("SNP", "bmi_ea", "bmi_nea", "bmi_beta", "bmi_se", "bmi_p", "bmi_chr", "bmi_pos", "bmi_eaf", 'bmi_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "bmi_beta",
                         se_col = "bmi_se",
                         eaf_col = "bmi_eaf",
                         effect_allele_col = "bmi_ea",
                         other_allele_col = "bmi_nea",pos_col	="pos",
                         chr_col="chr",pval_col="bmi_p", 
                         samplesize_col = "bmi_ss")
bmi_har<-harmonise_data(exposuredat, outcomedat, action = 2)
bmi_har<-bmi_har[which(bmi_har$mr_keep==TRUE),]

# WHR
exposuredata<-data[,c("SNP", "whr_ea", "whr_nea", "whr_beta", "whr_se", "whr_p", "whr_chr", "whr_pos", "whr_eaf", 'whr_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "whr_beta",
                         se_col = "whr_se",
                         eaf_col = "whr_eaf",
                         effect_allele_col = "whr_ea",
                         other_allele_col = "whr_nea",pos_col	="pos",
                         chr_col="chr",pval_col="whr_p", 
                         samplesize_col = "whr_ss")
whr_har<-harmonise_data(exposuredat, outcomedat, action = 2)
whr_har<-whr_har[which(whr_har$mr_keep==TRUE),]

# Heig
exposuredata<-data[,c("SNP", "heig_ea", "heig_nea", "heig_beta", "heig_se", "heig_p", "heig_chr", "heig_pos", "heig_eaf", 'heig_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "heig_beta",
                         se_col = "heig_se",
                         eaf_col = "heig_eaf",
                         effect_allele_col = "heig_ea",
                         other_allele_col = "heig_nea",pos_col	="pos",
                         chr_col="chr",pval_col="heig_p", 
                         samplesize_col = "heig_ss")
heig_har<-harmonise_data(exposuredat, outcomedat, action = 2)
heig_har<-heig_har[which(heig_har$mr_keep==TRUE),]

# wbfm
exposuredata<-data[,c("SNP", "wbfm_ea", "wbfm_nea", "wbfm_beta", "wbfm_se", "wbfm_p", "chr", "pos", "wbfm_eaf", 'wbfm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wbfm_beta",
                         se_col = "wbfm_se",
                         eaf_col = "wbfm_eaf",
                         effect_allele_col = "wbfm_ea",
                         other_allele_col = "wbfm_nea",
                         pval_col="wbfm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wbfm_ss")
wbfm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wbfm_har<-wbfm_har[which(wbfm_har$mr_keep==TRUE),]

# wblm
exposuredata<-data[,c("SNP", "wblm_ea", "wblm_nea", "wblm_beta", "wblm_se", "wblm_p", "chr", "pos", "wblm_eaf", 'wblm_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "wblm_beta",
                         se_col = "wblm_se",
                         eaf_col = "wblm_eaf",
                         effect_allele_col = "wblm_ea",
                         other_allele_col = "wblm_nea",
                         pval_col="wblm_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "wblm_ss")
wblm_har<-harmonise_data(exposuredat, outcomedat, action = 2)
wblm_har<-wblm_har[which(wblm_har$mr_keep==TRUE),]


# weigh
exposuredata<-data[,c("SNP", "weigh_ea", "weigh_nea", "weigh_beta", "weigh_se", "weigh_p", "chr", "pos", "weigh_eaf", 'weigh_ss')]
exposuredat<-format_data(exposuredata,type="exposure",snp_col = "SNP",
                         beta_col = "weigh_beta",
                         se_col = "weigh_se",
                         eaf_col = "weigh_eaf",
                         effect_allele_col = "weigh_ea",
                         other_allele_col = "weigh_nea",
                         pval_col="weigh_p", pos_col	="pos",
                         chr_col="chr",
                         samplesize_col = "weigh_ss")
weigh_har<-harmonise_data(exposuredat, outcomedat, action = 2)
weigh_har<-weigh_har[which(weigh_har$mr_keep==TRUE),]

rm(exposuredat,exposuredata,outcomedat,outcomedata)

genelist<- gsub("_iv","_har", genelist, fixed = TRUE) #make list with names of all harmonised sets


# Clump
clp<-function(dat){
  dat <- get(dat, envir = .GlobalEnv)
  dat<-clump_data(dat, clump_r2=0.001, pop="EUR")
}
iv_list<- sapply(genelist, clp, simplify = FALSE)
invisible(lapply(names(iv_list), function(x) assign(x, iv_list[[x]],envir=.GlobalEnv)))

# Do MR
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  b<-0
  se<-0
  pval<-0
  n.SNP<-0
  
  if(nrow(dat)==1){
    dat$effect_allele_col<-"A"
    dat$other_allele_col<-"C"
    exposuredat<-format_data(dat,type="exposure",snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             pval_col="pval.exposure",
                             eaf_col = "eaf.exposure",
                             effect_allele_col = "effect_allele_col",
                             other_allele_col = "other_allele_col")
    
    
    outcomedat<-format_data(dat,type="outcome",snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            pval_col="pval.outcome",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele_col",
                            other_allele_col = "other_allele_col")
    
    harmonizeddata<-harmonise_data(exposuredat, outcomedat, action = 2)
    output<-mr_singlesnp(harmonizeddata,single_method="mr_wald_ratio")
    b<-output[1,7]
    se<-output[1,8]
    pval<-output[1,9]
    n.SNP<-1
    
  }else{
    output<-TwoSampleMR::mr_ivw(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
    b<-unlist(output)[1]
    se<-unlist(output)[2]
    pval<-unlist(output)[3]
    n.SNP<-unlist(output)[4]
  }
  
  outputobj<-list(b,se,pval,n.SNP)
  return(unlist(outputobj))
}

mr_table2 <- t(data.frame((sapply(genelist,mrfunc2))))
# mr_table2$or<-exp(mr_table2$coef)
# mr_table2$lci<-exp((mr_table2$coef)-((mr_table2$se)*1.96))
# mr_table2$uci<-exp((mr_table2$coef)+((mr_table2$se)*1.96))
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/pwtforce_results2.csv')

# Sentitivity analyses 
#wm 
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_weighted_median(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
mr_table2 <- as.data.frame((sapply(genelist,mrfunc2)))
mr_table2<-t(mr_table2)
write.csv(mr_table2, '~/Desktop/bmi-ecg/res/pwtforce_results2_wm.csv')

#egger
mrfunc2<-function(dat) {
  dat <- get(dat, envir = .GlobalEnv)
  mr_egger_regression(b_exp = dat$beta.exposure, se_exp = dat$se.exposure, b_out = dat$beta.outcome, se_out = dat$se.outcome)
}
a <- (sapply(genelist,mrfunc2))
write.csv(a, '~/Desktop/bmi-ecg/res/pwtforce_results2_egg.csv')

# 
rm(list=ls())
# #




