### Multiple Hypothesis Testing Correction (FDR) of all results; Formatting phenotype descriptions
#   Performed separately for both ICD10/Finngen traits and continuous lab traits
 
# Load packages 
library(tidyverse)
library(data.table)
 
# Load all results where p<0.05 across health traits
df <- read_csv('data/UKBBv3_AllTraits_FINAL_AllAssociations.csv.gz')
 
# Load cell death gene list of interest dataframe
pws <- read_csv('metadata/ProgrammedCellDeathGeneList_Pathways.csv')
genes <- pws$`Gene Name` 
 
# Load phenotypes of interest dataframe
pheno <- read_csv('metadata/UKBBv3_AllTraits2Category_FINAL_v2_plusCaseControl.csv')
 
# Subset dataframe for genes and phenotypes of interest
df <- filter(df, gene_name %in% genes, `Phenotype Code` %in% phenos)
df <- merge(df, pheno, by=c("Phenotype Code", "Derivation"), all.x=T)
df <- merge(df, pws, by.x="gene_name", by.y="Gene Name", all.x=T)
df <- select(df, gene_name, Pathway, SubPathway, `Phenotype Description`, Category, zscore, 
             effect_size, pvalue, Tissue, `Phenotype Code`, gene, 
             `Protein Name`,var_g, pred_perf_r2, pred_perf_pval, pred_perf_qval, 
             n_snps_used, n_snps_in_cov, n_snps_in_model, Derivation, File)
df$Tissue <- str_remove(df$Tissue, "JTI_") %>% str_replace_all("_", " ")
 
 
# FDR Correction, separately on each dataset -------
# Separate the formatted data
df_ehr_f <- filter(df, Derivation == "ICD10/FinnGen")
df_cont_f <- filter(df, Derivation == "Lab/Continuous")
 
# Determine count of number of tests performed
JTItissues_genes <- read_csv('metadata/GenesTestedForEachJTITissue.csv')
JTItissues_cdgenes <- filter(JTItissues_genes, gene_name %in% genes)
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
df_fdr <- df_fdr[,c(1:8,22,9:21)]
 
write_csv(df_fdr, "data/UKBBv3_SPrediXcanResults_PCDGenes_FINAL.csv")
