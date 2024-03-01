### Neale Lab LDSR Heritability for chosen vs. omitted traits

library(tidyverse)

# Load heritability results
# original file available at 'https://nealelab.github.io/UKBB_ldsc/downloads.html' as described in Data Availability section of manuscript
h2 <- read_tsv('heritability/ukb31063_h2_all.02Oct2019.tsv.gz') 
h2 <- filter(h2, sex == "both_sexes")

# Load phenotype lists (all & curated)
# original phenotype manifest file available from Neale Lab v3 release. See Data Availability section, or, 'https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=227859291'
allpheno <- read_tsv('metadata/phenotype selection and review/phenotypes.both_sexes.tsv.bgz')
chosenpheno <- read_csv('metadata/UKBBv3_AllTraits2Category_FINAL_v2.csv')
chosenphenovec <- chosenpheno$`Phenotype Code`
omitphenovec <- filter(allpheno, !phenotype %in% chosenphenovec)
omitphenovec <- omitphenovec$phenotype


# Annotate heriability df with h2 values for chosen phenotypes. 
h2$included <- ifelse(h2$phenotype %in% chosenphenovec, "include", "exclude")
h2chosen <- filter(h2, phenotype %in% chosenphenovec)
h2omit <- filter(h2, !phenotype %in% chosenphenovec)
n_chosen <- nrow(h2chosen)
n_omit <- nrow(h2omit)
nh_chosen <- filter(h2chosen, h2_p<0.05 & h2_observed>0) %>% nrow()
nh_omit <- filter(h2omit, h2_p<0.05 & h2_observed>0) %>% nrow()
nh_chosen/n_chosen
nh_omit/n_omit
# the % of heritable phenotypes is not greater in chosen vs omitted


# For ICD10-derived "clinical outcome" traits
h2_heritable <- filter(h2, h2_p<0.05 & h2_observed>0)
h2_heritable_icd10 <- filter(h2_heritable, source == "icd10")

p1.icd <- ggplot(h2_heritable_icd10, aes(x=-log10(h2_p), y=(h2_observed), color=included)) +
  geom_point()
p1.icd 

p2.icd <- ggplot(h2_heritable_icd10, aes(x=h2_observed, y=included)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Observed Heritability (h2)") +
  ylab("Inclusion Status of Trait")
p2.icd # The average heritability of chosen ICD10-derived traits is higher than that of excluded ICD10 traits
ggsave(filename="/Users/abbyrich/Desktop/CellDeath/figures/2023Feb06/HeritabilityOfIncludedVsExcludedTraits.tiff",
       plot=p5, device="tiff", units="in", height=3, width=5)

h2_heritable_icd10_chosen <- filter(h2_heritable_icd10, included == "include")
h2_heritable_icd10_omit <- filter(h2_heritable_icd10, included == "exclude")
t.test(h2_heritable_icd10_chosen$h2_observed, h2_heritable_icd10_omit$h2_observed, alternative="greater")
# 1-tailed t test p=0.06243


# Additional comparisons if of interest -------

# For all traits, both clincal outcomes and continuous 
h2_heritable <- filter(h2, h2_p<0.05 & h2_observed>0)
p1.all <- ggplot(h2_heritable, aes(x=-log10(h2_p), y=(h2_observed), color=included)) +
  geom_point()
p1.all  
p2.all <- ggplot(h2_heritable, aes(x=h2_observed, y=included)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Observed Heritability (h2)") +
  ylab("Inclusion Status of Trait")
p2.all # the significance level & heritability scores of all chosen traits vs omitted traits not more sig/larger than omitted
h2_heritable_all_chosen <- filter(h2_heritable, included == "include")
h2_heritable_all_omit <- filter(h2_heritable, included == "exclude")
t.test(h2_heritable_all_chosen$h2_observed, h2_heritable_all_omit$h2_observed, alternative="greater")
# 1-tailed t test p-value=1

# For only EHR-derived traits (both ICD10 and FinnGen-derived)
h2_heritable_ehr <- filter(h2_heritable, source %in% c("icd10","finngen"))
p1.ehr <- ggplot(h2_heritable_ehr, aes(x=-log10(h2_p), y=(h2_observed), color=included)) +
  geom_point()
p1.ehr 
p2.ehr <- ggplot(h2_heritable_ehr, aes(x=h2_observed, y=included)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Observed Heritability (h2)") +
  ylab("Inclusion Status of Trait")
p2.ehr
h2_heritable_ehr_chosen <- filter(h2_heritable_ehr, included == "include")
h2_heritable_ehr_omit <- filter(h2_heritable_ehr, included == "exclude")
t.test(h2_heritable_ehr_chosen$h2_observed, h2_heritable_ehr_omit$h2_observed, alternative="greater")
# 1-tailed t test p=0.9035

