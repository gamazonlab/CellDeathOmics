### Enrichment of Cell Death Gene/Trait Associations Across Pathways
 
library(tidyverse)
 
# Load both ICD10/FinnGen traits and Lab/continuous traits
phenos <- read_csv('data/UKBBv3_AllTraits2Category_FINAL_v2_plusCaseControl.csv')
 
# Load genes of interest & their pathway annotations
pws <- read_csv('data/ProgrammedCellDeathGeneList_Pathways.csv')
  
# Load associations for all genes, tissues, and traits where p<0.05
df_all <- read_csv('data/UKBBv3_AllTraits_FINAL_AllAssociations.csv.gz')
 
 
# FDR calculations for each trait across all tissues -------
# The purpose of this is to identify the p-value threshold at which an FDR<0.05 holds for 
# multiple testing correction done using hte 
# Load information on genes tested in each JTI tissue model
genepertissue <- read_csv('data/GenesTestedForEachJTITissue.csv')
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
 
 
# Fisher's Exact tests for enrichment of apoptosis genes across all traits -------
apogenes <- filter(pws, Pathway %in% "APOPTOSIS")
apogenes <- apogenes$`Gene Name`
TestEnrichmentOnePheno_Apo <- function(X) {
  genepertissue_goi <- filter(genepertissue, gene_name %in% apogenes)
  sig_goi <- filter(df, `Phenotype Code` %in% X) %>% filter(gene_name %in% apogenes) %>% nrow()
  sig_notgoi <- filter(df, `Phenotype Code` %in% X) %>% filter(!gene_name %in% apogenes) %>% nrow()
  notsig_goi <- nrow(genepertissue_goi) - sig_goi # subtract all significant PCD associations from all PCD associations
  notsig_notgoi <- nrow(genepertissue) - nrow(genepertissue_goi) - sig_notgoi # subtract all significant non PCD associations from all non PCD associations
  mat <- matrix(c(sig_goi, notsig_goi, sig_notgoi, notsig_notgoi), nrow=2, ncol=2, )
  fisher <- fisher.test(mat, alternative = "greater")
  pval <- fisher[["p.value"]]
  res <- data.frame(`Phenotype Code`=X, pvalue=pval)
  return(res)
}
reslist2 <- lapply(traits, TestEnrichmentOnePheno_Apo)
apo_enrich <- do.call(rbind, reslist2)
apo_enrich <- merge(apo_enrich, phenos, by.x="Phenotype.Code", by.y="Phenotype Code", all.x=T)
apo_enrich <- arrange(apo_enrich, pvalue)
apo_enrich$GeneGrouping <- "Apoptosis Genes"
 
 
# Fisher's Exact tests for enrichment of necroptosis genes across all traits -------
necgenes <- filter(pws, Pathway %in% "NECROPTOSIS")
necgenes <- necgenes$`Gene Name`
TestEnrichmentOnePheno_Nec <- function(X) {
  genepertissue_goi <- filter(genepertissue, gene_name %in% necgenes)
  sig_goi <- filter(df, `Phenotype Code` %in% X) %>% filter(gene_name %in% necgenes) %>% nrow()
  sig_notgoi <- filter(df, `Phenotype Code` %in% X) %>% filter(!gene_name %in% necgenes) %>% nrow()
  notsig_goi <- nrow(genepertissue_goi) - sig_goi # subtract all significant PCD associations from all PCD associations
  notsig_notgoi <- nrow(genepertissue) - nrow(genepertissue_goi) - sig_notgoi # subtract all significant non PCD associations from all non PCD associations
  mat <- matrix(c(sig_goi, notsig_goi, sig_notgoi, notsig_notgoi), nrow=2, ncol=2, )
  fisher <- fisher.test(mat, alternative = "greater")
  pval <- fisher[["p.value"]]
  res <- data.frame(`Phenotype Code`=X, pvalue=pval)
  return(res)
}
reslist3 <- lapply(traits, TestEnrichmentOnePheno_Nec)
nec_enrich <- do.call(rbind, reslist3)
nec_enrich <- merge(nec_enrich, phenos, by.x="Phenotype.Code", by.y="Phenotype Code", all.x=T)
nec_enrich <- arrange(nec_enrich, pvalue)
nec_enrich$GeneGrouping <- "Necroptosis Genes"
 
 
# Fisher's Exact tests for enrichment of pyroptosis genes across all traits -------
pyrogenes <- filter(pws, Pathway %in% "PYROPTOSIS")
pyrogenes <- pyrogenes$`Gene Name`
TestEnrichmentOnePheno_Pyro <- function(X) {
  genepertissue_goi <- filter(genepertissue, gene_name %in% pyrogenes)
  sig_goi <- filter(df, `Phenotype Code` %in% X) %>% filter(gene_name %in% pyrogenes) %>% nrow()
  sig_notgoi <- filter(df, `Phenotype Code` %in% X) %>% filter(!gene_name %in% pyrogenes) %>% nrow()
  notsig_goi <- nrow(genepertissue_goi) - sig_goi # subtract all significant PCD associations from all PCD associations
  notsig_notgoi <- nrow(genepertissue) - nrow(genepertissue_goi) - sig_notgoi # subtract all significant non PCD associations from all non PCD associations
  mat <- matrix(c(sig_goi, notsig_goi, sig_notgoi, notsig_notgoi), nrow=2, ncol=2, )
  fisher <- fisher.test(mat, alternative = "greater")
  pval <- fisher[["p.value"]]
  res <- data.frame(`Phenotype Code`=X, pvalue=pval)
  return(res)
}
reslist4 <- lapply(traits, TestEnrichmentOnePheno_Pyro)
pyro_enrich <- do.call(rbind, reslist4)
pyro_enrich <- merge(pyro_enrich, phenos, by.x="Phenotype.Code", by.y="Phenotype Code", all.x=T)
pyro_enrich <- arrange(pyro_enrich, pvalue)
pyro_enrich$GeneGrouping <- "Pyroptosis Genes"
 
 
# Fisher's exact test for enrichment for categories in gene sets
pathways <- unique(pws$Pathway)
categorypathwaycombos <- expand.grid(phenos$Category, pathways, stringsAsFactors = FALSE) %>% unique()
combolist <- t(categorypathwaycombos) %>% as.data.frame %>% as.list.data.frame(categorypathwaycombos)
 
TestEnrichmentCategoryOfPhenos <- function(X) {
  xCategory = X[1]
  xPathway = X[2]
  genes <- filter(pws, Pathway %in% xPathway)
  genes <- genes$`Gene Name`
  pheno <- filter(phenos, Category %in% xCategory)
  pheno <- pheno$`Phenotype Code`
  genepertissue_goi <- filter(genepertissue, gene_name %in% genes)
  sig_goi <- filter(df, `Phenotype Code` %in% pheno) %>% filter(gene_name %in% genes) %>% nrow()
  sig_notgoi <- filter(df, `Phenotype Code` %in% pheno) %>% filter(!gene_name %in% genes) %>% nrow()
  notsig_goi <- (nrow(genepertissue_goi))*(length(pheno)) - sig_goi # subtract all significant PCD associations from all PCD associations
  notsig_notgoi <- (nrow(genepertissue))*(length(pheno)) - (nrow(genepertissue_goi))*(length(pheno)) - sig_notgoi # subtract all significant non PCD associations from all non PCD associations
  mat <- matrix(c(sig_goi, notsig_goi, sig_notgoi, notsig_notgoi), nrow=2, ncol=2, )
  fisher <- fisher.test(mat)
  pval <- fisher[["p.value"]]
  res <- data.frame(Category=xCategory, Pathway=xPathway, pvalue=pval)
  return(res)
}
reslist5 <- lapply(combolist, TestEnrichmentCategoryOfPhenos)
cat_enrich <- do.call(rbind, reslist5)
cat_enrich <- arrange(cat_enrich, pvalue)
cat_enrich$FDR <- p.adjust(cat_enrich_FDR$pvalue, method = "BH")
 
 
# Fisher's exact test for enrichment for categories for all PCD genes
categories <- unique(phenos$Category)
categorylist <- as.vector(categories, mode="list")
 
TestEnrichmentCategoryOfPhenos <- function(X) {
  pheno <- filter(phenos, Category %in% X)
  pheno <- pheno$`Phenotype Code`
  genepertissue_goi <- filter(genepertissue, gene_name %in% pws$`Gene Name`)
  sig_goi <- filter(df, `Phenotype Code` %in% pheno) %>% filter(gene_name %in% pws$`Gene Name`) %>% nrow()
  sig_notgoi <- filter(df, `Phenotype Code` %in% pheno) %>% filter(!gene_name %in% pws$`Gene Name`) %>% nrow()
  notsig_goi <- (nrow(genepertissue_goi))*(length(pheno)) - sig_goi # subtract all significant PCD associations from all PCD associations
  notsig_notgoi <- (nrow(genepertissue))*(length(pheno)) - (nrow(genepertissue_goi))*(length(pheno)) - sig_notgoi # subtract all significant non PCD associations from all non PCD associations
  mat <- matrix(c(sig_goi, notsig_goi, sig_notgoi, notsig_notgoi), nrow=2, ncol=2, )
  fisher <- fisher.test(mat)
  pval <- fisher[["p.value"]]
  res <- data.frame(Category=X, pvalue=pval)
  return(res)
}
reslist6 <- lapply(categorylist, TestEnrichmentCategoryOfPhenos)
cat_enrich_all <- do.call(rbind, reslist6)
cat_enrich_all <- arrange(cat_enrich_all, pvalue)
cat_enrich_all$FDR <- p.adjust(cat_enrich_all$pvalue, method = "BH")
cat_enrich_all$Pathway <- "ALL"
 
# Create enrichment result dataframe for trait categories
cat_enrich_both <- rbind(cat_enrich_all, cat_enrich)
cat_enrich_both2 <- cat_enrich_both %>% select(-pvalue) %>% pivot_wider(names_from=Pathway, values_from = c(FDR))
 
 
# Create, annotate, and save large enrichment result dataframe
all_enrich <- rbind(apo_enrich, nec_enrich, pyro_enrich)
all_enrich2 <- pivot_wider(all_enrich, names_from=GeneGrouping, values_from=pvalue)
 
# For individual groups only
all_enrich$FDR <- p.adjust(all_enrich$pvalue, method = "BH")
all_enrich <- arrange(all_enrich, FDR)
 
write_csv(all_enrich, 'data/TraitEnrichment_FishersExact_AllGeneGroupings.csv')
