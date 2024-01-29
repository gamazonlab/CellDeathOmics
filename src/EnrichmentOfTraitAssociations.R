### Enrichment of Cell Death Gene/Trait Associations Across Pathways
 
library(tidyverse)
 
# Load both ICD10/FinnGen traits and Lab/continuous traits
phenos <- read_csv('metadata/UKBBv3_AllTraits2Category_FINAL_v2_plusCaseControl.csv')
 
# Load genes of interest & their pathway annotations
pws <- read_csv('metadata/ProgrammedCellDeathGeneList_Pathways.csv')
  
# Load associations for all genes, tissues, and traits where p<0.05
df_all <- read_csv('data/UKBBv3_AllTraits_FINAL_AllAssociations.csv.gz')
 
 
# FDR calculations for each trait across all tissues -------
# The purpose of this is to identify the p-value threshold at which an FDR<0.05 holds for 
# multiple testing correction done using hte 
# Load information on genes tested in each JTI tissue model
genepertissue <- read_csv('metadata/GenesTestedForEachJTITissue.csv')
genepertissue_cd <- filter(genepertissue, gene_name %in% pws$`Gene Name`)
 
traits <- as.vector(phenos$`Phenotype Code`, mode="list")
df_cd <- filter(df_all, gene_name %in% pws$`Gene Name`)
df_trait <- df_all[1,]
df_trait$FDR_trait <- NA
df_trait <- df_trait[NA,]
 
FDRforTrait <- function(X) {
  dftemp <- filter(df_cd, `Phenotype Code` %in% X)
  ntests <- nrow(genepertissue_cd)
  dftemp$FDR_trait <- p.adjust(dftemp$pvalue, method="BH", n=ntests)
  df_trait <- rbind(df_trait, dftemp)
  df_trait[-1,]
}
df_trait_fdr <- lapply(traits, FDRforTrait)
df_trait_fdr <- do.call("rbind", df_trait_fdr)
df_trait_fdr_sig <- filter(df_trait_fdr, FDR_trait<0.05)
###### p-value cutoff corresponds to about 0.02
rm(df_trait, df_trait_fdr, df_trait_fdr_sig, FDRforTrait, traits)
 
# Fisher's Exact test for enrichment of PCD genes across all traits -------
  
# In below function, "df" represents all "significant" associations (p<0.02): 
df <- filter(df_all, pvalue<0.02)
 
TestEnrichmentOnePheno <- function(X) {
  # Significant and PCD gene
  sig_goi <- filter(df, `Phenotype Code` %in% X) %>% filter(gene_name %in% pws$`Gene Name`) %>% nrow()
  # Significant and not PCD gene
  sig_notgoi <- filter(df, `Phenotype Code` %in% X) %>% filter(!gene_name %in% pws$`Gene Name`) %>% nrow()
  # Not significant and PCD gene
  notsig_goi <- nrow(genepertissue_cd) - sig_goi # subtract all significant PCD associations from all PCD associations
  # Not significant and not PCD gene
  notsig_notgoi <- nrow(genepertissue) - nrow(genepertissue_cd) - sig_notgoi # subtract all significant non PCD associations from all non PCD associations
  # Matrix structure
  # (sig & goi) (not sig and goi) (sig and not goi) (not sig and not goi)
  mat <- matrix(c(sig_goi, notsig_goi, sig_notgoi, notsig_notgoi), nrow=2, ncol=2, )
  fisher <- fisher.test(mat, alternative = "greater")
  pval <- fisher[["p.value"]]
  res <- data.frame(`Phenotype Code`=X, pvalue=pval)
  return(res)
}
 
# Run function then combine, annotate, and save results
reslist <- lapply(traits, TestEnrichmentOnePheno)
pcd_enrich <- do.call(rbind, reslist)
pcd_enrich <- merge(pcd_enrich, phenos, by.x="Phenotype.Code", by.y="Phenotype Code", all.x=T)
pcd_enrich <- arrange(pcd_enrich, pvalue)
pcd_enrich$GeneGrouping <- "All PCD Genes"
pcd_enrich$FDR <- p.adjust(pcd_enrich$pvalue, method="BH")
write_csv(pcd_enrich, 'results/TraitEnrichment_FishersExact_PCDGenes.csv')
