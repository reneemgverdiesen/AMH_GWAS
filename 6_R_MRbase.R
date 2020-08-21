##########################
## MR analyses AMH GWAS ##
##                      ##
##    Renée Verdiesen   ##
##########################


# Install and load needed packages
# install.packages('data.table')
library(data.table)
# install.packages('devtools')
library(devtools)
# install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
# install.packages('ggplot2')
library(ggplot2)


# This script contains the following steps:
# 1. Prepare summary statistics of AMH GWAS meta-analysis for Mendelian Randomization analyses
# 2. Perform MR analyses for:
#    2.1 Breast cancer
#    2.2 PCOS



# Step 1: Prepare summary statistics AMH GWAS (own GWAS data)

# Load AMH summary stats for four genome-wide significant loci into R
amh_dat <- read_exposure_data(
  filename = "/path/to/AMH_GWAS/Meta-analysis output/AMH_MRBaseformat_main.txt",
  sep = '\t',
  snp_col = 'SNP',
  beta_col = 'beta',
  se_col = 'se',
  effect_allele_col = 'effect_allele',
  phenotype_col = '',
  units_col = '',
  other_allele_col = 'other_allele',
  eaf_col = 'eaf',
  samplesize_col = '',
  ncase_col = '',
  ncontrol_col = '',
  gene_col = '',
  pval_col = 'pval'
) # 46 SNPs with p-value < 5e-8 = correct

# Perform clumping of these 46 AMH GWAS SNPs - use same thresholds as we did for FUMA and DEPICT analyses
amh_dat <- clump_data(amh_dat,
                      clump_kb = 500,
                      clump_r2 = 0.1) # I checked if default setting have different results -- this is not the case

amh_dat # correct, these are the 4 identified lead SNPs


# Change summary statistics to AMH increasing effects --> flip SNP rs10417628, rs13009019 and rs11683493
amh_dat_flip <- amh_dat
amh_dat_flip$effect_allele.exposure2 <- amh_dat_flip$effect_allele.exposure
amh_dat_flip$other_allele.exposure2 <- amh_dat_flip$other_allele.exposure
# Betas
amh_dat_flip$beta.exposure[amh_dat_flip$SNP == "rs10417628"] <- abs(amh_dat_flip$beta.exposure)[amh_dat_flip$SNP == "rs10417628"]
amh_dat_flip$beta.exposure[amh_dat_flip$SNP == "rs13009019"] <- abs(amh_dat_flip$beta.exposure)[amh_dat_flip$SNP == "rs13009019"]
amh_dat_flip$beta.exposure[amh_dat_flip$SNP == "rs11683493"] <- abs(amh_dat_flip$beta.exposure)[amh_dat_flip$SNP == "rs11683493"]
# Effect alleles
amh_dat_flip$effect_allele.exposure[amh_dat_flip$SNP == "rs10417628"] <- amh_dat_flip$other_allele.exposure2[amh_dat_flip$SNP == "rs10417628"]
amh_dat_flip$effect_allele.exposure[amh_dat_flip$SNP == "rs13009019"] <- amh_dat_flip$other_allele.exposure2[amh_dat_flip$SNP == "rs13009019"]
amh_dat_flip$effect_allele.exposure[amh_dat_flip$SNP == "rs11683493"] <- amh_dat_flip$other_allele.exposure2[amh_dat_flip$SNP == "rs11683493"]
# Other alleles
amh_dat_flip$other_allele.exposure[amh_dat_flip$SNP == "rs10417628"] <- amh_dat_flip$effect_allele.exposure2[amh_dat_flip$SNP == "rs10417628"]
amh_dat_flip$other_allele.exposure[amh_dat_flip$SNP == "rs13009019"] <- amh_dat_flip$effect_allele.exposure2[amh_dat_flip$SNP == "rs13009019"]
amh_dat_flip$other_allele.exposure[amh_dat_flip$SNP == "rs11683493"] <- amh_dat_flip$effect_allele.exposure2[amh_dat_flip$SNP == "rs11683493"]
# Change eaf
amh_dat_flip$eaf.exposure[amh_dat_flip$SNP == "rs10417628"] <- (1-amh_dat_flip$eaf.exposure)[amh_dat_flip$SNP == "rs10417628"]
amh_dat_flip$eaf.exposure[amh_dat_flip$SNP == "rs13009019"] <- (1-amh_dat_flip$eaf.exposure)[amh_dat_flip$SNP == "rs13009019"]
amh_dat_flip$eaf.exposure[amh_dat_flip$SNP == "rs11683493"] <- (1-amh_dat_flip$eaf.exposure)[amh_dat_flip$SNP == "rs11683493"]

amh_dat
amh_dat_flip # correct

amh_dat_flip <- amh_dat_flip[, -c(13:14)] 
# Conlcusion: flipping to all positive betas does not change results -- but we use flipped data because interpretation is more intuitive

# Remove original data from R environment
rm(amh_dat)



# Step 2.1 MR analysis AMH and breast cancer

# We used breast cancer summary statistics corresponding to the manuscript of K Michailidou et al. (2017)
# These data were downloaded on 29.11.2019 from http://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-combined-summary-result/gwas-summary-results-breast-cancer-risk-2017
# Extract only AMH SNPs from breast cancer outcome file
bcac_dat <- read_outcome_data(
  filename = "/path/to/AMH GWAS/C Data/MR_analyses/Summarystatistics_outcomes/BCAC_MR.txt.gz",
  sep = "\t",
  phenotype_col = "",
  snp_col = "snp",
  snps = c("rs10417628", "rs13009019", "rs11683493", "rs16991615"),
  beta_col = "bcac_onco_icogs_gwas_beta",
  se_col = "bcac_onco_icogs_gwas_se",
  eaf_col = "bcac_onco_icogs_gwas_eaf_controls",
  effect_allele_col = "a1",
  other_allele_col = "a0",
  pval_col = "bcac_onco_icogs_gwas_P1df",
  units_col = "",
  ncase_col = "",
  ncontrol_col = "",
  samplesize_col = ""
)


# Harmonize data
dat_amh_bc <- harmonise_data(amh_dat_flip, 
                             bcac_dat, 
                             action = 2)


# Perform MR analysis AMH and breast cancer
results <- mr(dat_amh_bc, 
              method_list "mr_ivw") # mr_ivw --> random effects IVW
generate_odds_ratios(results) 
clean_results = generate_odds_ratios(results)[, c(5, 7:14)]
clean_results[2:9] <- round(clean_results[, c(2:9)], digits = 2)

# Make plot
p1 <- mr_scatter_plot(results, dat_amh_bc)
p1[[1]]
# Save plot
# ggsave(filename = "AMHGWAS_MR_BC_scatterplot.png",
#       plot = p1[[1]],
#       path = "/path/to/AMH GWAS/D Output/1 Article_report/MR_analyses/Breast cancer/",
#       dpi = 1200,
#       limitsize = F)


# Perform single SNP analysis
singlesnp <- mr_singlesnp(dat_amh_bc)
generate_odds_ratios(singlesnp)
clean_singlesnp <- generate_odds_ratios(singlesnp)[, c(6:14)]
clean_singlesnp[2:9] <- round(clean_singlesnp[, c(2:9)], digits = 2)
clean_singlesnp
# Make forrest plot single snp analysis
p2 <- mr_forest_plot(singlesnp, exponentiate = T)
p2[[1]]
# Save plot
# ggsave(filename = "AMHGWAS_MR_BC_singlesnp_forrestplot.png",
#       plot = p2[[1]],
#       path = "/path/to/AMH GWAS/D Output/1 Article_report/MR_analyses/Breast cancer/",
#       dpi = 1200,
#       limitsize = F)


# Perform leave one out analysis
# IVW method
 leaveoneout_ivw <- mr_leaveoneout(dat_amh_bc, method = mr_ivw)
 generate_odds_ratios(leaveoneout_ivw)
 p3 <- mr_leaveoneout_plot(leaveoneout_ivw)
 p3[[1]]
# Save plot
# ggsave(filename = "AMHGWAS_MR_BC_leaveoneout_ivw.png",
#       plot = p3[[1]],
#       path = "/path/to/AMH GWAS/D Output/1 Article_report/MR_analyses/Breast cancer/",
#       dpi = 1200,
#       limitsize = F)


# Assess heterogeneity
mr_heterogeneity(dat_amh_bc) 


# Save results
# write.table(clean_results,
#            file = "/path/to/AMH GWAS/D Output/1 Article_report/MR_analyses/Breast cancer/AMHGWAS_final_Estimates.csv",
#            row.names = F,
#            quote = F,
#            sep = ";")


# Remove breast cancer objects from environment
bc <- c(ls())
bc <- bc[-1]
rm(list = bc, bc)


# Restart R environment to empty working memory
.rs.restartR()



# Step 2.2 MR analysis AMH and PCOS
# We used PCOS summary statistics corresponding to the manuscript of F Day et al. (2018), 
# These data were downloaded on 02-12-2019 from https://www.repository.cam.ac.uk/handle/1810/283491
# Extract only AMH SNPs from PCOS outcome file
pcos_dat <- read_outcome_data(
  filename = "/path/to/AMH GWAS/C Data/MR_analyses/Summarystatistics_outcomes/PCOS/PCOS_summary_data_19092018.txt",
  sep = "auto",
  phenotype_col = "",
  snp_col = "MarkerName",
  snps = c("20:5948227", "19:2251817", "2:145670572", "2:174259325"),
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "EAF",
  effect_allele_col = "Effect_allele",
  other_allele_col = "Other_allele",
  pval_col = "Pvalue",
  units_col = "",
  ncase_col = "",
  ncontrol_col = "",
  samplesize_col = "TotalSampleSize"
)

pcos_dat

# Make variable with rsids for pcos_dat
pcos_dat$SNP_rs[pcos_dat$SNP=="2:145670572"] <- "rs13009019"
pcos_dat$SNP_rs[pcos_dat$SNP=="20:5948227"] <- "rs16991615"
pcos_dat$SNP_rs[pcos_dat$SNP=="2:174259325"] <- "rs11683493"
pcos_dat$SNP_rs[pcos_dat$SNP=="19:2251817"] <- "rs10417628"

# Remove variable SNP
pcos_dat$SNP <- NULL

# Rename SNP_rs
colnames(pcos_dat)[13] <- "SNP"


# Harmonize data
dat_amh_pcos <- harmonise_data(amh_dat_flip, 
                             pcos_dat, 
                             action = 2)

# Perform MR analysis AMH - PCOS
results <- mr(dat_amh_pcos, 
              method_list "mr_ivw")
generate_odds_ratios(results)
clean_results = generate_odds_ratios(results)[, c(5, 7:14)]
clean_results[2:9] <- round(clean_results[, c(2:9)], digits = 2)

# Make plot
p1 <- mr_scatter_plot(results, dat_amh_pcos)
p1[[1]]
# Save plot
# ggsave(filename = "AMHGWAS_MR_PCOS_scatterplot.png",
#       plot = p1[[1]],
#       path = "/path/to/AMH GWAS/D Output/1 Article_report/MR_analyses/PCOS/",
#       dpi = 1200,
#       limitsize = F)


# Perform single SNP analysis
singlesnp <- mr_singlesnp(dat_amh_pcos)
generate_odds_ratios(singlesnp)
clean_singlesnp <- generate_odds_ratios(singlesnp)[, c(6:14)]
clean_singlesnp[2:9] <- round(clean_singlesnp[, c(2:9)], digits = 2)
# Make forrest plot singl snp analysis
p2 <- mr_forest_plot(singlesnp, exponentiate = T)
p2[[1]]
# Save plot
# ggsave(filename = "AMHGWAS_MR_PCOS_singlesnp_forrestplot.png",
#       plot = p2[[1]],
#       path = "/path/to/AMH GWAS/D Output/1 Article_report/MR_analyses/PCOS/",
#       dpi = 1200,
#       limitsize = F)


# Perform leave one out analysis
# IVW method
leaveoneout_ivw <- mr_leaveoneout(dat_amh_pcos, method = mr_ivw)
generate_odds_ratios(leaveoneout_ivw)
p3 <- mr_leaveoneout_plot(leaveoneout_ivw)
p3[[1]]
# Save plot
# ggsave(filename = "AMHGWAS_MR_PCOS_leaveoneout_ivw.png",
#       plot = p3[[1]],
#       path = "/path/to/AMH GWAS/D Output/1 Article_report/MR_analyses/PCOS/",
#       dpi = 1200,
#       limitsize = F)


# Assess heterogeneity
mr_heterogeneity(dat_amh_pcos)


# Save results
# write.table(clean_results,
#            file = "/path/to/AMH GWAS/D Output/1 Article_report/MR_analyses/PCOS/AMHGWAS_final_Estimates.csv",
#            row.names = F,
#            quote = F,
#            sep = ";")



# Remove all objects from R environment
rm(list=ls())



### THE END ###
