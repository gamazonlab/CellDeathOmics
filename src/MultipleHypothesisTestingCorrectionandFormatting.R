### Multiple Testing Correction (FDR)
#   Performed separately for both ICD10/Finngen traits and continuous lab traits

### UKBB PheWAS Multiple Hypothesis Testing Correction -------

# Load packages 
library(tidyverse)
library(data.table)

# Load all results where p<0.05 across health traits
df <- read_csv('results/UKBBv3_SPrediXcanResults_PCDGenes_FINAL.csv')

# Load cell death gene list of interest dataframe
pws <- read_csv('metadata/ProgrammedCellDeathGeneList_Pathways.csv')
genes <- pws$`Gene Name` 

# Load phenotypes of interest dataframe
pheno <- read_csv('metadata/UKBBv3_AllTraits2Category_FINAL_v2_plusCaseControl.csv')
pheno_ehr <- filter(pheno, Derivation == "ICD10/FinnGen")
pheno_cont <- filter(pheno, Derivation == "Lab/Continuous")

# FDR Correction, separately on each dataset -------
# Separate the formatted data
df_ehr_f <- filter(df, Derivation == "ICD10/FinnGen")
df_cont_f <- filter(df, Derivation == "Lab/Continuous")

# Determine count of number of tests performed
JTItissues_cdgenes <- read_csv('metadata/GenesTestedForEachJTITissue_PCDonly.csv')
JTItissues_cdgenects <- count(JTItissues_cdgenes, gene_name) %>% arrange(desc(n)) # 5 genes not tested: ACSL4, NOX1/2/3, XIAP
ntests_ehr <- (sum(JTItissues_cdgenects$n))*(nrow(pheno_ehr))
ntests_cont <- (sum(JTItissues_cdgenects$n))*(nrow(pheno_cont))

# FDR calculations for all genes and phenotypes of interest
df_ehr_f$FDR_all <- p.adjust(df_ehr_f$pvalue, method="BH", n=ntests_ehr)
df_ehr_f <- df_ehr_f %>% arrange(FDR_all) # arrange by FDR
df_cont_f$FDR_all <- p.adjust(df_cont_f$pvalue, method="BH", n=ntests_cont)
df_cont_f <- df_cont_f %>% arrange(FDR_all) # arrange by FDR

# Harmonize and format dataframe for saving
df_fdr <- rbind(df_ehr_f, df_cont_f)

write_csv(df_fdr, "results/UKBBv3_SPrediXcanResults_PCDGenes_FINAL.csv")
# Note: FDR calculations were made in the same fashion for the ferroptosis/cuproptosis/parthanatos analysis


### BioVU Validation Multiple Hypothesis Testing Correction -------

# Load data
df <- read_csv('results/BioVU_SPrediXcanResults_PCDGenes.csv')

# Annotate with corresponding trait name
pheno <- read_csv('metadata/Dennis2021_LabGWAS_TraitList_v2_Annotated.csv')
pheno <- select(pheno, Phecode, `Phenotype Description`) %>% drop_na()
df <- merge(df, pheno, by.x="trait", by.y="Phecode")

# Determine count of number of tests performed
ntests_validation <- nrow(df)

# FDR calculations for all genes and phenotypes of interest
df$FDR_all <- p.adjust(df$pvalue, method="BH", n=ntests_validation)
df <- df %>% arrange(FDR_all) # arrange by FDR

# Harmonize and format dataframe for saving
df$tissue <- str_replace_all(df$tissue, "_", " ")

# Save dataframe
write_csv(df, 'results/BioVU_SPrediXcanResults_PCDGenes.csv')


### XIAP Results Multiple Hypothesis Testing Correction (Bonferroni) -------

# Load data
df <- read_csv('results/UKBB_SPrediXcanResults_XIAP.csv')

# Determine count of number of tests performed
ntests_ehr = 2014 # 2014 tests (19 tissues * 106 phenotypes)
bonferroni = 0.05/1629

# Subset for results that meet significance
dfsig <- filter(df, pvalue < bonferroni)
