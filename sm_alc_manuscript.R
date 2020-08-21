setwd("/Volumes/015/working/data/hnc_dpw_si/manuscript")
install.packages("devtools")
library(devtools)
devtools::install_github("MRCIEU/TwoSampleMR@0.4.26")
library(TwoSampleMR)
devtools::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
devtools::install_github("WSpiller/MVMR")
library(MVMR)
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
install.packages("MendelianRandomization")
library(MendelianRandomization)

#Table 1 
#Univariable MR of smoking initiation with HNC 
exposure_dat <-read_exposure_data("SI_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- read_outcome_data("smandalc_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
SI_SNPs <- as.character(dat$SNP[dat$mr_keep=="TRUE"])
write.csv(SI_SNPs, "SI_SNPs_final.csv", quote=F)
dat$exposure <- "smoking initiation"
dat$outcome <- "oral/oropharyngeal cancer"
dat <- dat[dat$mr_keep=="TRUE",]
#176 SNPs - Supplementary Table 1 
write.csv(dat, "SI_dat.csv")

mr_results <- mr(dat)
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "si_hnc.txt")

#MR-PRESSO 
dat <- dat[dat$mr_keep=="TRUE",]
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", data=dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 10000,  SignifThreshold = 0.05)

#Univariable MR for SI w/o UK Biobank and 23andMe
exposure_dat <- read_exposure_data("smandalc_si_noukb.txt", sep="\t")
exposure_dat <- exposure_dat[exposure_dat$SNP %in% SI_SNPs,]
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Head and neck cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "si_hnc_woukb.txt")

dat <- dat[dat$mr_keep=="TRUE",]
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", data=dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 10000,  SignifThreshold = 0.05)

#Univariable MR for SI in Biobank only 
exposure_dat <-read_exposure_data("si_ukb.csv", sep=",")
exposure_dat$beta.exposure <- exposure_dat$beta.exposure / ((exposure_dat$ncase.exposure/exposure_dat$samplesize.exposure)*(1-(exposure_dat$ncase.exposure/exposure_dat$samplesize.exposure)))
exposure_dat$se.exposure <- exposure_dat$se.exposure / ((exposure_dat$ncase.exposure/exposure_dat$samplesize.exposure)*(1-(exposure_dat$ncase.exposure/exposure_dat$samplesize.exposure)))
outcome_dat <- read_outcome_data("smandalc_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Head and neck cancer"
dat <- dat[dat$mr_keep=="TRUE",]
mr_results <- mr(dat)
mr_results$b = mr_results$b / ((dat$ncase.exposure/(dat$ncase.exposure+dat$ncontrol.exposure))*(1-(dat$ncase.exposure/(dat$ncase.exposure+dat$ncontrol.exposure))))[1] 
mr_results$se = mr_results$se / ((dat$ncase.exposure/(dat$ncase.exposure+dat$ncontrol.exposure))*(1-(dat$ncase.exposure/(dat$ncase.exposure+dat$ncontrol.exposure))))[1] 
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "si_hnc_biobank.txt")

dat <- dat[dat$mr_keep=="TRUE",]
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", data=dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 10000,  SignifThreshold = 0.05)

#Univariable MR with CSI and HNC 
exposure_dat <-read_exposure_data("CSI_exposure.csv", sep=",")
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- read_outcome_data("smandalc_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
CSI_SNPs <- as.character(dat$SNP[dat$mr_keep=="TRUE"])
write.csv(CSI_SNPs, "CSI_SNPs.csv", quote=F)
dat$exposure <- "comprehensive smoking index"
dat$outcome <- "oral/oropharyngeal cancer"
dat <- dat[dat$mr_keep=="TRUE",]
#108 SNPs - Supplementary Table 1 
write.csv(dat, "CSI_dat.csv")
dat$samplesize.exposure <- 462690
mr_results <- mr(dat)
#to get the outcome per standard deviation increase then multiply beta and se by 0.6940093
mr_results$b <- mr_results$b*0.6940093
mr_results$se <- mr_results$se*0.6940093
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "csi_hnc.txt")

dat <- dat[dat$mr_keep=="TRUE",]
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", data=dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 10000,  SignifThreshold = 0.05)

#Table 2
#Univariable MR with drinks per week and HNC (units from paper)
exposure_dat <-read_exposure_data("DPW_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- read_outcome_data("smandalc_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
DPW_SNPs <- as.character(dat$SNP[dat$mr_keep=="TRUE"])
write.csv(DPW_SNPs, "DPW_SNPs_final.csv")
SI_DPW_SNPs <- c(SI_SNPs, DPW_SNPs)
write.csv(SI_DPW_SNPs, "SI_DPW_SNPs_final.csv")
dat <- dat[dat$mr_keep=="TRUE",]
dat$outcome <- "Head and neck cancer"
#60 SNPs - Supplementary Table 1 
write.csv(dat, "DPW_dat.csv")
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_hnc.txt")

dat <- dat[dat$mr_keep=="TRUE",]
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", data=dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 10000,  SignifThreshold = 0.05)

#Univariable MR for drinks per week and HNC (excluding ADH1B)
exposure_dat <-read_exposure_data("DPW_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
exposure_dat <- exposure_dat[exposure_dat$SNP != "rs1229984",]
outcome_dat <- read_outcome_data("smandalc_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[dat$mr_keep=="TRUE",]
dat$outcome <- "Head and neck cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_hnc_noadh1b.txt")

#Univariable MR for DPW w/o UK Biobank and 23andMe
exposure_dat <- read_exposure_data("smandalc_dpw_noukb.txt", sep="\t")
exposure_dat <- exposure_dat[exposure_dat$SNP %in% DPW_SNPs,]
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Head and neck cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_hnc_woukb.txt")

dat <- dat[dat$mr_keep=="TRUE",]
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", data=dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 10000,  SignifThreshold = 0.05)

#Univariable MR for DPW in Biobank only 
exposure_dat <- read_exposure_data("smandalc_dpw_ukb.txt", sep="\t")
exposure_dat <- exposure_dat[exposure_dat$SNP %in% DPW_SNPs,]
outcome_dat <- read_outcome_data("smandalc_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[dat$mr_keep=="TRUE",]
dat$outcome <- "Head and neck cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_hnc_biobank.txt")

dat <- dat[dat$mr_keep=="TRUE",]
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", data=dat, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 10000,  SignifThreshold = 0.05)

#Univariable MR for DPW in UK Biobank only (excluding ADH1B)
exposure_dat <- read_exposure_data("smandalc_dpw_ukb.txt", sep="\t")
exposure_dat <- exposure_dat[exposure_dat$SNP %in% DPW_SNPs,]
exposure_dat <- exposure_dat[exposure_dat$SNP != "rs1229984",]
outcome_dat <- read_outcome_data("smandalc_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[dat$mr_keep=="TRUE",]
dat$outcome <- "Head and neck cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_hnc_biobank_noadh1b.txt")

#Figure 2 
#Univariable MR of smoking initiation with OC and OPC 
exposure_dat <-read_exposure_data("SI_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- read_outcome_data("smandalc_oc.txt", sep="\t")
outcome_dat <- outcome_dat[outcome_dat$SNP %in% exposure_dat$SNP,]
dat <- harmonise_data(exposure_dat, outcome_dat)
SI_SNPs <- as.character(dat$SNP[dat$mr_keep=="TRUE"])
dat$exposure <- "smoking initiation"
dat$outcome <- "oral cancer"
dat <- dat[dat$mr_keep=="TRUE",]
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "si_oc.txt")

exposure_dat <-read_exposure_data("SI_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- read_outcome_data("smandalc_opc.txt", sep="\t")
outcome_dat <- outcome_dat[outcome_dat$SNP %in% exposure_dat$SNP,]
dat <- harmonise_data(exposure_dat, outcome_dat)
SI_SNPs <- as.character(dat$SNP[dat$mr_keep=="TRUE"])
dat$exposure <- "smoking initiation"
dat$outcome <- "oropharyngeal cancer"
dat <- dat[dat$mr_keep=="TRUE",]
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "si_opc.txt")

#Univariable MR of CSI with OC and OPC 
exposure_dat <-read_exposure_data("CSI_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- read_outcome_data("smandalc_oc.txt", sep="\t")
outcome_dat <- outcome_dat[outcome_dat$SNP %in% exposure_dat$SNP,]
dat <- harmonise_data(exposure_dat, outcome_dat)
CSI_SNPs <- as.character(dat$SNP[dat$mr_keep=="TRUE"])
dat$exposure <- "comprehensive smoking index"
dat$outcome <- "oral cancer"
dat <- dat[dat$mr_keep=="TRUE",]
mr_results <- mr(dat)
#to get the outcome per standard deviation increase then multiply beta and se by 0.6940093
mr_results$b <- mr_results$b*0.6940093
mr_results$se <- mr_results$se*0.6940093
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "csi_oc.txt")

exposure_dat <-read_exposure_data("CSI_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- read_outcome_data("smandalc_opc.txt", sep="\t")
outcome_dat <- outcome_dat[outcome_dat$SNP %in% exposure_dat$SNP,]
dat <- harmonise_data(exposure_dat, outcome_dat)
CSI_SNPs <- as.character(dat$SNP[dat$mr_keep=="TRUE"])
dat$exposure <- "comprehensive smoking index"
dat$outcome <- "oropharyngeal cancer"
dat <- dat[dat$mr_keep=="TRUE",]
mr_results <- mr(dat)
#to get the outcome per standard deviation increase then multiply beta and se by 0.6940093
mr_results$b <- mr_results$b*0.6940093
mr_results$se <- mr_results$se*0.6940093
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "csi_opc.txt")

#Univariable MR of alcohol intake with OC and OPC 
exposure_dat <-read_exposure_data("DPW_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- read_outcome_data("smandalc_oc.txt", sep="\t")
outcome_dat <- outcome_dat[outcome_dat$SNP %in% exposure_dat$SNP,]
dat <- harmonise_data(exposure_dat, outcome_dat)
DPW_SNPs <- as.character(dat$SNP[dat$mr_keep=="TRUE"])
dat$exposure <- "drinks per week"
dat$outcome <- "oral cancer"
dat <- dat[dat$mr_keep=="TRUE",]
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_oc2.txt")

exposure_dat <-read_exposure_data("DPW_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- read_outcome_data("smandalc_opc.txt", sep="\t")
outcome_dat <- outcome_dat[outcome_dat$SNP %in% exposure_dat$SNP,]
dat <- harmonise_data(exposure_dat, outcome_dat)
DPW_SNPs <- as.character(dat$SNP[dat$mr_keep=="TRUE"])
dat$exposure <- "drinks per week"
dat$outcome <- "oropharyngeal cancer"
dat <- dat[dat$mr_keep=="TRUE",]
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_opc2.txt")

#Univariable MR of alcohol intake with OC and OPC (without ADH1B)
exposure_dat <-read_exposure_data("DPW_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
exposure_dat <- exposure_dat[exposure_dat$SNP != "rs1229984",]
outcome_dat <- read_outcome_data("smandalc_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[dat$mr_keep=="TRUE",]
dat$outcome <- "oral cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_oc_noadh1b.txt")

exposure_dat <-read_exposure_data("DPW_exposure.csv", sep=",", samplesize_col = "N")
exposure_dat <- clump_data(exposure_dat)
exposure_dat <- exposure_dat[exposure_dat$SNP != "rs1229984",]
outcome_dat <- read_outcome_data("smandalc_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[dat$mr_keep=="TRUE",]
dat$outcome <- "oropharyngeal cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_opc_noadh1b.txt")

#Univariable MR for DPW with OC and OPC in UK Biobank only (excluding ADH1B)
exposure_dat <- read_exposure_data("smandalc_dpw_ukb.txt", sep="\t")
exposure_dat <- exposure_dat[exposure_dat$SNP %in% DPW_SNPs,]
exposure_dat <- exposure_dat[exposure_dat$SNP != "rs1229984",]
outcome_dat <- read_outcome_data("smandalc_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[dat$mr_keep=="TRUE",]
dat$outcome <- "oral cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_hnc_biobank_oc_noadh1b.txt")

exposure_dat <- read_exposure_data("smandalc_dpw_ukb.txt", sep="\t")
exposure_dat <- exposure_dat[exposure_dat$SNP %in% DPW_SNPs,]
exposure_dat <- exposure_dat[exposure_dat$SNP != "rs1229984",]
outcome_dat <- read_outcome_data("smandalc_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[dat$mr_keep=="TRUE",]
dat$outcome <- "oropharyngeal cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
write.table(or_results, "dpw_hnc_biobank_opc_noadh1b.txt")

#Table 3 
#MVMR with CSI in UKB and DPW in GSCAN (-UKB)
XGs <- read.csv("MVMR_dpwandcsi.csv")
YG <- read.csv("MVMR_dpwandcsi_hnc.csv")

#restrict to 108 SNPs for CSI and 60 SNPs for DPW 
CSI_DPW_SNPs <- read.csv("CSI_DPW_SNPs.csv")
XGs <- XGs[XGs$SNP %in% CSI_DPW_SNPs$SNP,]
YG <- YG[YG$SNP %in% CSI_DPW_SNPs$SNP,]

library(tidyr)
XGs_betas <- XGs[,c(1:2,4)]
XGs_betas <- spread(XGs_betas, Exposure, xg)
XGs_se <- XGs[,c(1,3:4)]
XGs_se <- spread(XGs_se, Exposure, xgse)

#Remove NAs 
XGs_betas <- na.omit(XGs_betas)
XGs_se <- na.omit(XGs_se)

YG <- YG[YG$SNP %in% XGs_betas$SNP,]
XGs_betas <- XGs_betas[XGs_betas$SNP %in% YG$SNP,]
XGs_se <- XGs_se[XGs_se$SNP %in% YG$SNP,]

XGs_betas <- XGs_betas[order(XGs_betas$SNP),]
XGs_se <- XGs_se[order(XGs_se$SNP),]
YG <- YG[order(YG$SNP),]

mvmr <- format_mvmr(XGs_betas[,c(2:3)], YG$yg, XGs_se[,c(2:3)], YG$ygse, XGs_betas$SNP)
mvmr_res <- mvmr(mvmr, 0, 1)
save(mvmr_res, file="mvmr_res_dpwandcsi_final.Rdata")

#MVMR-Egger 
mr_mvivw <- mr_mvivw(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI), bxse = cbind(XGs_se$DPW, XGs_se$CSI), by = YG$yg, YG$ygse))
save(mr_mvivw, file="mr_mvivw_dpwandcsi_final.Rdata")

mr_mvegger <- mr_mvegger(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI), bxse = cbind(XGs_se$DPW, XGs_se$CSI), by = YG$yg, YG$ygse), orientate = 1)
save(mr_mvegger, file="mr_mvegger_dpwandcsi_final.Rdata")

mvmr_results_DPW <- c(exp(mr_mvivw$Estimate[1]), exp(mr_mvivw$CILower[1]), exp(mr_mvivw$CIUpper[1]), mr_mvivw$Pvalue[1], exp(mr_mvegger$Estimate[1]), exp(mr_mvegger$CILower.Est[1]), exp(mr_mvegger$CIUpper.Est[1]), mr_mvegger$Pvalue.Est[1])
mvmr_results_DPW <- c(format(mvmr_results_DPW, scientific=F))
names(mvmr_results_DPW) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_CSI <- c(exp(mr_mvivw$Estimate[2]*0.6940093), exp(mr_mvivw$CILower[2]*0.6940093), exp(mr_mvivw$CIUpper[2]*0.6940093), mr_mvivw$Pvalue[2], exp(mr_mvegger$Estimate[2]*0.6940093), exp(mr_mvegger$CILower.Est[2]*0.6940093), exp(mr_mvegger$CIUpper.Est[2]*0.6940093), mr_mvegger$Pvalue.Est[2])
mvmr_results_CSI <- c(format(mvmr_results_CSI, scientific=F))
names(mvmr_results_CSI) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results <- cbind(mvmr_results_DPW, mvmr_results_CSI)  
#to get the outcome per standard deviation increase then multiply beta and se by 0.6940093
write.csv(mvmr_results, "mvmr_results.csv")    

mvmregger_results_CSI <- c(mr_mvegger$Intercept[1]*0.6940093, mr_mvegger$StdError.Int[1]*0.6940093, mr_mvegger$CILower.Int[1]*0.6940093, mr_mvegger$CIUpper.Int[1]*0.6940093, mr_mvegger$Pvalue.Int[1])
names(mvmregger_results_CSI) <- c("intercept", "stderror.int", "cilower.int", "ciupper.int", "pvalue.int")

#Instrument strength 
strength_mvmr <- strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr <- pleiotropy_mvmr(mvmr, gencov=0)

#without ADH1B
XGs_betas <- XGs_betas[XGs_betas$SNP != "rs1229984",]
XGs_se <- XGs_se[XGs_se$SNP != "rs1229984",]

YG <- YG[YG$SNP %in% XGs_betas$SNP,]
XGs_betas <- XGs_betas[XGs_betas$SNP %in% YG$SNP,]
XGs_se <- XGs_se[XGs_se$SNP %in% YG$SNP,]

mvmr <- format_mvmr(XGs_betas[,c(2:3)], YG$yg, XGs_se[,c(2:3)], YG$ygse, XGs_betas$SNP)
mvmr_res <- mvmr(mvmr, 0, 1)

save(mvmr_res, file="mvmr_res_dpwandcsi_noadh1b_final.Rdata")

#MVMR-Egger 
mr_mvivw <- mr_mvivw(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI), bxse = cbind(XGs_se$DPW, XGs_se$CSI), by = YG$yg, YG$ygse))
save(mr_mvivw, file="mr_mvivw_dpwandcsi_noadh1b_final.Rdata")

mr_mvegger <- mr_mvegger(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI), bxse = cbind(XGs_se$DPW, XGs_se$CSI), by = YG$yg, YG$ygse), orientate = 1)
save(mr_mvegger, file="mr_mvegger_dpwandcsi_noadh1b_final.Rdata")

mvmr_results_DPW <- c(exp(mr_mvivw$Estimate[1]), exp(mr_mvivw$CILower[1]), exp(mr_mvivw$CIUpper[1]), mr_mvivw$Pvalue[1], exp(mr_mvegger$Estimate[1]), exp(mr_mvegger$CILower.Est[1]), exp(mr_mvegger$CIUpper.Est[1]), mr_mvegger$Pvalue.Est[1])
mvmr_results_DPW <- c(format(mvmr_results_DPW, scientific=F))
names(mvmr_results_DPW) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_CSI <- c(exp(mr_mvivw$Estimate[2]*0.6940093), exp(mr_mvivw$CILower[2]*0.6940093), exp(mr_mvivw$CIUpper[2]*0.6940093), mr_mvivw$Pvalue[2], exp(mr_mvegger$Estimate[2]*0.6940093), exp(mr_mvegger$CILower.Est[2]*0.6940093), exp(mr_mvegger$CIUpper.Est[2]*0.6940093), mr_mvegger$Pvalue.Est[2])
mvmr_results_CSI <- c(format(mvmr_results_CSI, scientific=F))
names(mvmr_results_CSI) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_noadh1b <- cbind(mvmr_results_DPW, mvmr_results_CSI)                 
write.csv(mvmr_results_noadh1b, "mvmr_results_noadh1b.csv")    

mvmregger_results_CSI <- c(mr_mvegger$Intercept[1]*0.6940093, mr_mvegger$StdError.Int[1]*0.6940093, mr_mvegger$CILower.Int[1]*0.6940093, mr_mvegger$CIUpper.Int[1]*0.6940093, mr_mvegger$Pvalue.Int[1])
names(mvmregger_results_CSI) <- c("intercept", "stderror.int", "cilower.int", "ciupper.int", "pvalue.int")

#Instrument strength 
strength_mvmr <- strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr <- pleiotropy_mvmr(mvmr, gencov=0)

#Figure 3 
#MVMR with CSI in UKB and DPW in GSCAN (-UKB) for OPC 
XGs <- read.csv("MVMR_dpwandcsi.csv")
YG <- read.csv("MVMR_dpwandcsi_opc.csv")

#restrict to 108 SNPs for CSI and 60 SNPs for DPW 
CSI_DPW_SNPs <- read.csv("CSI_DPW_SNPs.csv")
XGs <- XGs[XGs$SNP %in% CSI_DPW_SNPs$SNP,]
YG <- YG[YG$SNP %in% CSI_DPW_SNPs$SNP,]

library(tidyr)
XGs_betas <- XGs[,c(1:2,4)]
XGs_betas <- spread(XGs_betas, Exposure, xg)
XGs_se <- XGs[,c(1,3:4)]
XGs_se <- spread(XGs_se, Exposure, xgse)

#Remove NAs 
XGs_betas <- na.omit(XGs_betas)
XGs_se <- na.omit(XGs_se)

YG <- YG[YG$SNP %in% XGs_betas$SNP,]
XGs_betas <- XGs_betas[XGs_betas$SNP %in% YG$SNP,]
XGs_se <- XGs_se[XGs_se$SNP %in% YG$SNP,]

XGs_betas <- XGs_betas[order(XGs_betas$SNP),]
XGs_se <- XGs_se[order(XGs_se$SNP),]
YG <- YG[order(YG$SNP),]

#without ADH1B
XGs_betas <- XGs_betas[XGs_betas$SNP != "rs1229984",]
XGs_se <- XGs_se[XGs_se$SNP != "rs1229984",]

YG <- YG[YG$SNP %in% XGs_betas$SNP,]
XGs_betas <- XGs_betas[XGs_betas$SNP %in% YG$SNP,]
XGs_se <- XGs_se[XGs_se$SNP %in% YG$SNP,]

mvmr <- format_mvmr(XGs_betas[,c(2:3)], YG$yg, XGs_se[,c(2:3)], YG$ygse, XGs_betas$SNP)
mvmr_res <- mvmr(mvmr, 0, 1)

save(mvmr_res, file="mvmr_res_dpwandcsi_noadh1b_opc_final.Rdata")

#MVMR-Egger 
mr_mvivw <- mr_mvivw(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI), bxse = cbind(XGs_se$DPW, XGs_se$CSI), by = YG$yg, YG$ygse))
save(mr_mvivw, file="mr_mvivw_dpwandcsi_noadh1b_opc_final.Rdata")

mr_mvegger <- mr_mvegger(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI), bxse = cbind(XGs_se$DPW, XGs_se$CSI), by = YG$yg, YG$ygse), orientate = 1)
save(mr_mvegger, file="mr_mvegger_dpwandcsi_noadh1b_opc_final.Rdata")

mvmr_results_DPW <- c(exp(mr_mvivw$Estimate[1]), exp(mr_mvivw$CILower[1]), exp(mr_mvivw$CIUpper[1]), mr_mvivw$Pvalue[1], exp(mr_mvegger$Estimate[1]), exp(mr_mvegger$CILower.Est[1]), exp(mr_mvegger$CIUpper.Est[1]), mr_mvegger$Pvalue.Est[1])
mvmr_results_DPW <- c(format(mvmr_results_DPW, scientific=F))
names(mvmr_results_DPW) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_CSI <- c(exp(mr_mvivw$Estimate[2]*0.6940093), exp(mr_mvivw$CILower[2]*0.6940093), exp(mr_mvivw$CIUpper[2]*0.6940093), mr_mvivw$Pvalue[2], exp(mr_mvegger$Estimate[2]*0.6940093), exp(mr_mvegger$CILower.Est[2]*0.6940093), exp(mr_mvegger$CIUpper.Est[2]*0.6940093), mr_mvegger$Pvalue.Est[2])
mvmr_results_CSI <- c(format(mvmr_results_CSI, scientific=F))
names(mvmr_results_CSI) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_opc <- cbind(mvmr_results_DPW, mvmr_results_CSI)                 
write.csv(mvmr_results_opc, "mvmr_results_opc.csv")    

mvmregger_results_CSI <- c(mr_mvegger$Intercept[1]*0.6940093, mr_mvegger$StdError.Int[1]*0.6940093, mr_mvegger$CILower.Int[1]*0.6940093, mr_mvegger$CIUpper.Int[1]*0.6940093, mr_mvegger$Pvalue.Int[1])
names(mvmregger_results_CSI) <- c("intercept", "stderror.int", "cilower.int", "ciupper.int", "pvalue.int")

#Instrument strength 
strength_mvmr <- strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr <- pleiotropy_mvmr(mvmr, gencov=0)

#MVMR with CSI in UKB and DPW in GSCAN (-UKB) for OC 
XGs <- read.csv("MVMR_dpwandcsi.csv")
YG <- read.csv("MVMR_dpwandcsi_oc.csv")

#restrict to 108 SNPs for CSI and 60 SNPs for DPW 
CSI_DPW_SNPs <- read.csv("CSI_DPW_SNPs.csv")
XGs <- XGs[XGs$SNP %in% CSI_DPW_SNPs$SNP,]
YG <- YG[YG$SNP %in% CSI_DPW_SNPs$SNP,]

library(tidyr)
XGs_betas <- XGs[,c(1:2,4)]
XGs_betas <- spread(XGs_betas, Exposure, xg)
XGs_se <- XGs[,c(1,3:4)]
XGs_se <- spread(XGs_se, Exposure, xgse)

#Remove NAs 
XGs_betas <- na.omit(XGs_betas)
XGs_se <- na.omit(XGs_se)

YG <- YG[YG$SNP %in% XGs_betas$SNP,]
XGs_betas <- XGs_betas[XGs_betas$SNP %in% YG$SNP,]
XGs_se <- XGs_se[XGs_se$SNP %in% YG$SNP,]

XGs_betas <- XGs_betas[order(XGs_betas$SNP),]
XGs_se <- XGs_se[order(XGs_se$SNP),]
YG <- YG[order(YG$SNP),]

#without ADH1B
XGs_betas <- XGs_betas[XGs_betas$SNP != "rs1229984",]
XGs_se <- XGs_se[XGs_se$SNP != "rs1229984",]

YG <- YG[YG$SNP %in% XGs_betas$SNP,]
XGs_betas <- XGs_betas[XGs_betas$SNP %in% YG$SNP,]
XGs_se <- XGs_se[XGs_se$SNP %in% YG$SNP,]

mvmr <- format_mvmr(XGs_betas[,c(2:3)], YG$yg, XGs_se[,c(2:3)], YG$ygse, XGs_betas$SNP)
mvmr_res <- mvmr(mvmr, 0, 1)

save(mvmr_res, file="mvmr_res_dpwandcsi_noadh1b_oc_final.Rdata")

#MVMR-Egger 
mr_mvivw <- mr_mvivw(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI), bxse = cbind(XGs_se$DPW, XGs_se$CSI), by = YG$yg, YG$ygse))
save(mr_mvivw, file="mr_mvivw_dpwandcsi_noadh1b_oc_final.Rdata")

mr_mvegger <- mr_mvegger(mr_mvinput(bx = cbind(XGs_betas$DPW, XGs_betas$CSI), bxse = cbind(XGs_se$DPW, XGs_se$CSI), by = YG$yg, YG$ygse), orientate = 1)
save(mr_mvegger, file="mr_mvegger_dpwandcsi_noadh1b_oc_final.Rdata")

mvmr_results_DPW <- c(exp(mr_mvivw$Estimate[1]), exp(mr_mvivw$CILower[1]), exp(mr_mvivw$CIUpper[1]), mr_mvivw$Pvalue[1], exp(mr_mvegger$Estimate[1]), exp(mr_mvegger$CILower.Est[1]), exp(mr_mvegger$CIUpper.Est[1]), mr_mvegger$Pvalue.Est[1])
mvmr_results_DPW <- c(format(mvmr_results_DPW, scientific=F))
names(mvmr_results_DPW) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_CSI <- c(exp(mr_mvivw$Estimate[2]*0.6940093), exp(mr_mvivw$CILower[2]*0.6940093), exp(mr_mvivw$CIUpper[2]*0.6940093), mr_mvivw$Pvalue[2], exp(mr_mvegger$Estimate[2]*0.6940093), exp(mr_mvegger$CILower.Est[2]*0.6940093), exp(mr_mvegger$CIUpper.Est[2]*0.6940093), mr_mvegger$Pvalue.Est[2])
mvmr_results_CSI <- c(format(mvmr_results_CSI, scientific=F))
names(mvmr_results_CSI) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_oc <- cbind(mvmr_results_DPW, mvmr_results_CSI)                 
write.csv(mvmr_results_oc, "mvmr_results_oc.csv")    

mvmregger_results_CSI <- c(mr_mvegger$Intercept[1]*0.6940093, mr_mvegger$StdError.Int[1]*0.6940093, mr_mvegger$CILower.Int[1]*0.6940093, mr_mvegger$CIUpper.Int[1]*0.6940093, mr_mvegger$Pvalue.Int[1])
names(mvmregger_results_CSI) <- c("intercept", "stderror.int", "cilower.int", "ciupper.int", "pvalue.int")

#Instrument strength 
strength_mvmr <- strength_mvmr(mvmr, gencov=0)
pleiotropy_mvmr <- pleiotropy_mvmr(mvmr, gencov=0)

#Figure 4 - risk taking and number of sexual partners 
exposure_dat <- read_exposure_data("risktolerance.txt", sep="\t")
outcome_dat <- read_outcome_data("risk_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "risk_hnc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
write.table(or_results, "risk_hnc_results.txt")

outcome_dat <- read_outcome_data("risk_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "risk_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
write.table(or_results, "risk_oc_results.txt")

outcome_dat <- read_outcome_data("risk_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "risk_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
write.table(or_results, "risk_opc_results.txt")

exposure_dat <- read_exposure_data("nsp_exposure.txt", sep="\t")
outcome_dat <- read_outcome_data("nsp_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "nsp_hnc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
write.table(or_results, "nsp_hnc_results.txt")

outcome_dat <- read_outcome_data("nsp_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "nsp_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
write.table(or_results, "nsp_oc_results.txt")

outcome_dat <- read_outcome_data("nsp_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "nsp_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
write.table(or_results, "nsp_opc_results.txt")

#Supplementary Table 2
######################################################################################
# Regression dilution I2 GX
######################################################################################
# I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#calculate Isq wieghted and unweighted
I2<-c()

#Rename required columns
BetaXG <- dat$beta.exposure[dat$mr_keep=="TRUE"]
BetaYG <- dat$beta.outcome[dat$mr_keep=="TRUE"]
seBetaXG <- dat$se.exposure[dat$mr_keep=="TRUE"]
seBetaYG <- dat$se.outcome[dat$mr_keep=="TRUE"]
BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
I2<-cbind(mF, Isq_unweighted, Isq_weighted)
colnames(I2) <- c("mF", "Isq_unweighted", "Isq_weighted")
I2

#Run the simex correction for each phenotype
#install.packages("simex")
#load package
library(simex) 

#create empty dataframe to store output
simexegger<-c()

#run simex 
BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

#Supplementary Table 3
#Rucker framework 
mrrucker <- mr_rucker(dat)
#Rucker 1 - Start by performing IVW under fixed effects and calculate Q 
#Rucker 2 - If heterogeneity, switch to a random effects IVW 
#Rucker 3 - Calculate Q' from MR-Egger regression and calculate difference Q-Q'. If difference then switch to Egger
#Rucker 4 - If Q' reveals heterogeneity then switch to random effects MR-Egger 
mrrucker[[1]]$rucker
mrrucker[[1]]$intercept
mrrucker[[1]]$Q
mrrucker[[1]]$selected

#Supplementary Table 4
pleio <- mr_pleiotropy_test(dat)
#to get the outcome per standard deviation for CSI increase then multiply beta and se by 0.6940093
pleio$egger_intercept <- pleio$egger_intercept*0.6940093
pleio$se <- pleio$se*0.6940093
