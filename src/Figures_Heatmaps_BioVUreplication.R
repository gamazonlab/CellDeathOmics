# BioVU Blood Trait Validation - Overlap with UKBB results

library(tidyverse)

# Load data and subset for significant (FDR) or nominally significant (p-value) results
uk <- read_csv('results/UKBBv3_SPrediXcanResults_PCDGenes_FINAL.csv')
uk <- filter(uk, Derivation == "Lab/Continuous", FDR_all < 0.01)
vu <- read_csv('results/BioVU_SPrediXcanResults_PCDGenes.csv')
vu <- filter(vu, Derivation == "Lab/Continuous", pvalue < 0.05)


# Identify shared associations
uk$Direction <- ifelse(uk$zscore>0, "positive", "negative")
vu$Direction <- ifelse(vu$zscore>0, "positive", "negative")

uk2 <- select(uk, gene_name, `Phenotype Description`, Direction) %>% unique()
vu2 <- select(vu, gene_name, `Phenotype Description`, Direction) %>% unique()

harm <- inner_join(uk2, vu2, by = c("gene_name", "Phenotype Description", "Direction"))

# Filter the dataframes for these
uk_h <- merge(harm, uk, by=c("gene_name", "Phenotype Description", "Direction"), all.x=T)
vu_h <- merge(harm, vu, by=c("gene_name", "Phenotype Description", "Direction"), all.x=T)


# Heatmap of replicated associations -------
library(reshape2)

genes <- read_csv('metadata/ProgrammedCellDeathGeneList_Pathways.csv')
pheno <- read_csv('metadata/UKBBv3_AllTraits2Category_FINAL_v2_plusCaseControl.csv')

vu_h <- merge(vu_h, pheno, all.x=T)
vu_h <- merge(vu_h, genes, by.x="gene_name", by.y="Gene Name")

vu_h$SubPathway <- factor(vu_h$SubPathway, 
                            levels = c("PROSURVIVAL", "PROAPOPTOTIC BH3-ONLY", 
                                       "PROAPOPTOTIC EFFECTOR", "INITIATOR CASPASE", 
                                       "APOPTOSOME", 
                                       "EXECUTIONER CASPASE", "CASPASE ANTAGONIST",
                                       "NECROPTOSIS", "PYROPTOSIS"))

df.cont.med.all <-  vu_h %>% 
  group_by(gene_name, `Phenotype Description`) %>% 
  mutate(`Median Beta` = median(effect_size)) %>% 
  group_by(gene_name, `Phenotype Description`) %>% slice_min(FDR_all) %>%
  select(gene_name, `Phenotype Description`, Category, Pathway, SubPathway, `Median Beta`, pvalue, FDR_all) %>% 
  unique()
min.median.beta <- min(df.cont.med.all$`Median Beta`)
max.median.beta <- max(df.cont.med.all$`Median Beta`)
sig.cont.traits <- df.cont.med.all %>% filter(FDR_all < 0.05) %>% select(`Phenotype Description`) %>% unique()
sig.cont.traits <- sig.cont.traits$`Phenotype Description`

p.cont.trait.all.2 <- df.cont.med.all %>% 
  filter(`Phenotype Description` %in% sig.cont.traits) %>%
  select(gene_name, `Phenotype Description`, `Median Beta`, SubPathway, Category, FDR_all)
p.cont.trait.all.2b = subset(p.cont.trait.all.2, FDR_all<0.05)

p.cont.trait.all.heatplot <- ggplot(df.cont.med.all, aes(x = gene_name, y = `Phenotype Description`, group="FDR01")) +
  geom_tile(aes(fill=`Median Beta`), color="black", size=0.2) +
  geom_text(data = p.cont.trait.all.2b, aes(x = gene_name, y = `Phenotype Description`, label = "*", group="FDR01", size = 2, vjust = 0.5, hjust = 0.5)) +
  xlab("Gene Name") + ylab("Phenotype Description") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_bw() +
  theme(strip.text.x = element_blank(), strip.text.y=element_text(angle=0), axis.text.x=element_text(angle=90, hjust=1, face = "italic")) +
  facet_grid(Category ~ SubPathway, scales = "free", space = "free")
p.cont.trait.all.heatplot

ggsave('figures/Fig6_Heatplot_CBCs.tiff', 
       plot = p.cont.trait.all.heatplot, device = "tiff", 
       dpi = 200, width = 10, height = 4, units= "in", limitsize = TRUE)
