options(scipen = 999)
setwd("C:/Users/riocx/Documents/Masters new/Thesis/TSV/rsidmap/")
library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)
library(MVMR)
library(dplyr)
library(ggplot2)
library(gt)
data_path <- 'C:/Users/riocx/Documents/Masters new/Thesis/TSV/rsidmap/'
data_path1 <- 'C:/Users/riocx/Documents/Masters new/Thesis/TSV/rsidmap/implicated/'
perio <- vroom(paste0(data_path, "AcutePerioJiang_Adj.txt")) 
#PerioAC
#ChronicPerioJiang
#AcutePerioJiang
#AcutePerioUKBB
exposure <- vroom(paste0(data_path, "sorted_vitDdeficiency_ebi_header.txt"))#poultry.txt, meat_tophits_MAFadj_ukb.txt, poultry_tophitsMAF_ukb.txt
#Outcome data prep
perio <- na.omit(perio)
perio <- format_data(perio, type = "outcome",
                     snp_col = "SNP",
                     beta_col = "beta_EUR",
                     se_col = "se_EUR",
                     effect_allele_col = "ref",
                     other_allele_col = "alt",
                     eaf_col = "af_cases_EUR",
                     pval_col = "neglog10_pval_EUR") %>%
  mutate(outcome = 'Acute Perio')
#perio <- subset(perio, select = c("ref", "alt", "af_controls_EUR", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "rsid"))
#perio <- subset(perio, select = -exposure) #to remove coloumn from dataframe

#perio <- perio %>% filter(eaf.outcome > 0.99 | eaf.outcome < 0.01)
vroom_write(perio, 'AcutePerioUKBB_Adj.txt')
#Exposure data prep
exposure <- na.omit(exposure)
exposure <- format_data(exposure, type = "exposure",
                        snp_col = "variant_id",
                        beta_col = "beta",
                        se_col = "standard_error",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele",
                        eaf_col = "effect_allele_frequency",
                        pval_col = "p_value") %>%
  mutate(exposure = 'vitDdeficiency')
vroom_write(exposure, ('vitDdeficiency_ebi_Adj.txt'))
#snp_col = "SNP",
#beta_col = "beta_EUR",
#se_col = "se_EUR",
#effect_allele_col = "ref",
#other_allele_col = "alt",
#eaf_col = "af_cases_EUR",
#pval_col = "neglog10_pval_EUR"
#variant_id, base_pair_location, effect_allele, other_allele, n, effect_allele_frequency, beta, standard_error, p_value
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
tophits <- exposure %>% 
  filter(pval.exposure <= 5e-05) %>% #0.0000005
  clump_data(.,
             clump_kb = 10000,
             clump_r2 = 0.001,
             clump_p1 = 5e-05,
             clump_p2 = 1,
             pop = "EUR")
vroom_write(tophits, ('exposure_tophits.txt'))
#common_values <- intersect(poultry$SNP, perio$SNP)
#head(common_values)# Display common values
#alt way to extract_outcome_data:
#outcome_dat <- perio %>% filter(SNP %in% poultry$SNP) #no need for this?
#Harmonize data
Harmonized <- harmonise_data(tophits, perio) #check why action 2, it may be necessary to set action = 3
#check phenoscanner to ensure SNPs are not pleiotropic check page 4 fasting glucose levels and periodontitis
#library(phenoscanner)
#res <- phenoscanner(snpquery = test)
#head(res$results)
#res$snps
#MR analysis (with five models, alpha significance : 0.05)
mr(Harmonized, method_list = c("mr_egger_regression", "mr_ivw", "mr_ivw_radial",
                               "mr_ivw_mre", "mr_weighted_median"))
#Causal Effects Tests (in-depth analysis)
b.1test = mr_heterogeneity(Harmonized)
mr_heterogeneity(Harmonized)
#Pleiotropy_Test (H_0: There is no Pleiotropy alpha=0.05)
b.2test = mr_pleiotropy_test(Harmonized)
mr_pleiotropy_test(Harmonized)
#Leave-one-out Sensitivity_Test
b.3test = mr_leaveoneout(Harmonized) #delete the first snp from the bottom that comes up from the graph after checking pleiotropy from phenoscanner
mr_leaveoneout_plot(b.3test)
#MR Scatter plot #check the graph label
sp1 = mr(Harmonized)
mr_scatter_plot(sp1, Harmonized)
#MR Forest plot
mr_forest_plot(mr_singlesnp(Harmonized))
#Funnel plot (symmetry ? or Not?, symmetric plot represents a valid causal result.)
mr_funnel_plot(mr_singlesnp(Harmonized))

#Identify top variants by estimating effect sizes
