### Manhattan Plots For UKBB Cell Death PheWAS Results
#   For both continuous traits and clinical diagnoses
#   Updated on 2024JAN17

# Load packages
library(tidyverse)
library(ggrepel)
 
# Load data from original autosomal analysis of UKBB
unannotateddf <- read_csv('results/UKBBv3_SPrediXcanResults_PCDGenes_FINAL.csv')
 
# Load phenotype mapping, color, position information for graphing
phenos <- read_csv('metadata/categorycolors2phenotypes.csv')
 
# Annotate data with phenotype information for graphing
df <- merge(unannotateddf, phenos, by=c("Phenotype Code", "Phenotype Description", "Derivation"))
colnames(df)[colnames(df) == "Category.y"] <- "Category"
df_ehr <- filter(df, Derivation=="ICD10/FinnGen")
df_ehr <- df_ehr %>% group_by(gene_name, `Phenotype Description`) %>% 
  mutate(mostsigGenePheno = min(pvalue)) 
df_ehr <- df_ehr %>% group_by(gene_name, Category) %>% 
  mutate(mostsigGeneCat = min(pvalue))
df_ehr <- df_ehr %>% group_by(Category) %>% 
  mutate(mostsigCat = min(pvalue)) 
df_ehr <- df_ehr %>% group_by(Category) %>% mutate(xlabelpoint = median(Position)) %>%
  arrange(groupnum)
 
 
# Clinical trait Manhattan plot -------
  
# Top gene-trait associations per category given FDR < 0.25 cutoff
cutoff = 7.7e-5
labelcutoff=7.7e-5
df_ehr$LabelCat <- ifelse(df_ehr$pvalue %in% mostsigCat & df_ehr$pvalue < labelcutoff, paste0(df_ehr$gene_name, " - ", df_ehr$`Phenotype Description`), NA)

xBreaks <- unique(df_ehr$xlabelpoint)
xBreaks[15] <- 488 # requires custom nudging
xBreaks[16] <- 499
colors4graphing <- unique(df_ehr$color)
colors <- df_ehr %>% select(groupnum, color, Category) %>% 
  filter(color %in% colors4graphing) %>% unique() %>% arrange(groupnum)
catColors <- colors$color
names(catColors) = colors$Category
xLabels <- as.ordered(colors$Category)
 
df_ehr$Category <- as.ordered(df_ehr$Category)

# Manhattan plot
p_ehr_fdr <- ggplot(df_ehr, aes(x=Position, y = -log10(pvalue), color = Category, label=LabelCat)) +
  geom_point(size = 1.5) +
  geom_label_repel(size=3, min.segment.length = 0.5, max.overlaps = Inf, nudge_y = 5, force=40, force_pull=0.2, box.padding = unit(0.35, "lines"), fill = "white") +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust = 1),
        legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(values=catColors) +
  scale_x_continuous(labels = xLabels, breaks = xBreaks) +
  geom_hline(yintercept=-log10(cutoff), color = "red")+
  labs(y = "-log10(p-value)")
p_ehr_fdr
ggsave('Fig4A_Manhattan_EHR.tiff', plot=p_ehr_fdr, units="in", width=8.5, height=5)


# For Continuous Traits -------
df_lab <- filter(df, Derivation=="Lab/Continuous")
bloodtraits2omit <- c("High light scatter reticulocyte percentage", 
                        "Haematocrit percentage", "Lymphocyte percentage", 
                        "Monocyte percentage", "Neutrophil percentage",
                        "Reticulocyte percentage", "Mean sphered cell volume")

df_lab <- filter(df_lab, !`Phenotype Description` %in% bloodtraits2omit)
df_lab$absZ <- abs(df_lab$zscore)
df_lab <- df_lab %>% group_by(`Phenotype Description`) %>% 
  mutate(mostsigPheno = max(absZ))
df_lab <- df_lab %>% group_by(gene_name, `Phenotype Description`) %>% 
  mutate(mostsigGenePheno = max(absZ))
df_lab <- df_lab %>% group_by(gene_name, Category) %>% 
  mutate(mostsigGeneCat = max(absZ))
df_lab <- df_lab %>% group_by(Category) %>% 
  mutate(mostsigCat = max(absZ))
df_lab <- df_lab %>% group_by(Category) %>% mutate(xlabelpoint = median(Position)) %>%
  arrange(groupnum)

cutoff = 0.0008353 # FDR=0.01 
labelcutoff=2e-10
mostsigPheno <- unique(df_lab$mostsigPheno)
mostsigGenePheno <- unique(df_lab$mostsigGenePheno)
mostsigGeneCat <- unique(df_lab$mostsigGeneCat)
mostsigCat <- unique(df_lab$mostsigCat)
df_lab$LabelGenePheno <- ifelse(df_lab$absZ %in% mostsigGenePheno & df_lab$pvalue < labelcutoff, paste0(df_lab$gene_name, " - ", df_lab$`Phenotype Description`), NA)
df_lab$LabelPheno <- ifelse(df_lab$absZ %in% mostsigPheno & df_lab$pvalue < labelcutoff, paste0(df_lab$gene_name, " - ", df_lab$`Phenotype Description`), NA)
df_lab$LabelGeneCat <- ifelse(df_lab$absZ %in% mostsigGeneCat & df_lab$pvalue < labelcutoff, paste0(df_lab$gene_name, " - ", df_lab$`Phenotype Description`), NA)
df_lab$LabelCat <- ifelse(df_lab$absZ %in% mostsigCat & df_lab$pvalue < labelcutoff, paste0(df_lab$gene_name, " - ", df_lab$`Phenotype Description`), NA)

xBreaks <- unique(df_lab$xlabelpoint)
colors4graphing <- unique(df_lab$color)
colors <- df_lab %>% select(groupnum, color, Category) %>% 
  filter(color %in% colors4graphing) %>% unique() %>% arrange(groupnum)
catColors <- colors$color
names(catColors) = colors$Category
xLabels <- as.ordered(colors$Category)
df_lab$Category <- as.ordered(df_lab$Category)

# Manhattan plot
p_lab_bonferroni <- ggplot(df_lab, aes(x=Position, y = -log10(pvalue), color = Category, label=LabelPheno)) +
  geom_point(size = 1.5, position = "jitter") +
  geom_label_repel(size=3, min.segment.length = 1, max.overlaps = Inf, nudge_y = 15, force=40, force_pull=0.01) +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust = 1),
        legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ylim(0,100) +
  scale_color_manual(values=catColors) +
  scale_x_continuous(labels = xLabels, breaks = xBreaks) +
  geom_hline(yintercept=-log10(cutoff), color = "red")+
  labs(y = "-log10(p-value)")
p_lab_bonferroni
ggsave('Fig5A_Manhattan_Lab.tiff', plot=p_lab_bonferroni, units="in", width=8.5, height=5)



# For ferroptosis/cuproptosis/parthanatos genes -------
unannotateddf <- read_csv('UKBBv3_SPrediXcanResults_ExtraGenes_FINAL.csv')

# Load phenotype mapping, color, position information for graphing
phenos <- read_csv('metadata/categorycolors2phenotypes.csv')

# Annotate data with phenotype information for graphing
df <- merge(unannotateddf, phenos, by=c("Phenotype Code", "Phenotype Description", "Derivation"))
df <- select(df, -Category.x)
colnames(df)[colnames(df) == "Category.y"] <- "Category"
df_ehr <- filter(df, Derivation=="ICD10/FinnGen")
df_ehr <- df_ehr %>% group_by(gene_name, `Phenotype Description`) %>% 
  mutate(mostsigGenePheno = min(pvalue)) 
df_ehr <- df_ehr %>% group_by(gene_name, Category) %>% 
  mutate(mostsigGeneCat = min(pvalue))
df_ehr <- df_ehr %>% group_by(Category) %>% 
  mutate(mostsigCat = min(pvalue)) 
df_ehr <- df_ehr %>% group_by(Category) %>% mutate(xlabelpoint = median(Position)) %>%
  arrange(groupnum)

# Top gene-trait associations per category given FDR < 0.25 cutoff
cutoff = 7.0e-4
labelcutoff=7.0e-4
mostsigGenePheno <- unique(df_ehr$mostsigGenePheno)
mostsigGeneCat <- unique(df_ehr$mostsigGeneCat)
mostsigCat <- unique(df_ehr$mostsigCat)
df_ehr$LabelGenePheno <- ifelse(df_ehr$pvalue %in% mostsigGenePheno & df_ehr$pvalue < labelcutoff, paste0(df_ehr$gene_name, " - ", df_ehr$`Phenotype Description`), NA)
df_ehr$LabelGeneCat <- ifelse(df_ehr$pvalue %in% mostsigGeneCat & df_ehr$pvalue < labelcutoff, paste0(df_ehr$gene_name, " - ", df_ehr$`Phenotype Description`), NA)
df_ehr$LabelCat <- ifelse(df_ehr$pvalue %in% mostsigCat & df_ehr$pvalue < labelcutoff, paste0(df_ehr$gene_name, " - ", df_ehr$`Phenotype Description`), NA)

xBreaks <- unique(df_ehr$xlabelpoint)
xBreaks[15] <- 488 # requires custom nudging
xBreaks[16] <- 499
colors4graphing <- unique(df_ehr$color)
colors <- df_ehr %>% select(groupnum, color, Category) %>% 
  filter(color %in% colors4graphing) %>% unique() %>% arrange(groupnum)
catColors <- colors$color
names(catColors) = colors$Category
xLabels <- as.ordered(colors$Category)

df_ehr$Category <- as.ordered(df_ehr$Category)

# Manhattan plot
p_ehr_fdr <- ggplot(df_ehr, aes(x=Position, y = -log10(pvalue), color = Category, label=LabelGeneCat)) +
  geom_point(size = 1.5) +
  geom_label_repel(size=3, min.segment.length = 0.5, max.overlaps = Inf, nudge_y = 5, force=40, force_pull=0.2, box.padding = unit(0.35, "lines"), fill = "white") +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust = 1),
        legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(values=catColors) +
  scale_x_continuous(labels = xLabels, breaks = xBreaks) +
  geom_hline(yintercept=-log10(cutoff), color = "red")+
  labs(y = "-log10(p-value)")
p_ehr_fdr
ggsave('SuppFig3A_Manhattan_EHR.tiff', plot=p_ehr_fdr, units="in", width=8.5, height=5)


# For Lab/Continuous traits
df_lab <- filter(df, Derivation=="Lab/Continuous")

bloodtraits2omit <- c("High light scatter reticulocyte percentage", 
                        "Haematocrit percentage", "Lymphocyte percentage", 
                        "Monocyte percentage", "Neutrophil percentage",
                        "Reticulocyte percentage", "Mean sphered cell volume") # these traits are redundant

df_lab <- filter(df_lab, !`Phenotype Description` %in% bloodtraits2omit)
df_lab$absZ <- abs(df_lab$zscore)
df_lab <- df_lab %>% group_by(`Phenotype Description`) %>% 
  mutate(mostsigPheno = max(absZ))
df_lab <- df_lab %>% group_by(gene_name, `Phenotype Description`) %>% 
  mutate(mostsigGenePheno = max(absZ))
df_lab <- df_lab %>% group_by(gene_name, Category) %>% 
  mutate(mostsigGeneCat = max(absZ))
df_lab <- df_lab %>% group_by(Category) %>% 
  mutate(mostsigCat = max(absZ))
df_lab <- df_lab %>% group_by(Category) %>% mutate(xlabelpoint = median(Position)) %>%
  arrange(groupnum)

cutoff = 1.440092e-05 # Bonferroni (31 traits * 112 gene/tissue tests for the 3 added genes)
labelcutoff=1.440092e-05
mostsigPheno <- unique(df_lab$mostsigPheno)
mostsigGenePheno <- unique(df_lab$mostsigGenePheno)
mostsigGeneCat <- unique(df_lab$mostsigGeneCat)
mostsigCat <- unique(df_lab$mostsigCat)
df_lab$LabelGenePheno <- ifelse(df_lab$absZ %in% mostsigGenePheno & df_lab$pvalue < labelcutoff, paste0(df_lab$gene_name, " - ", df_lab$`Phenotype Description`), NA)
df_lab$LabelPheno <- ifelse(df_lab$absZ %in% mostsigPheno & df_lab$pvalue < labelcutoff, paste0(df_lab$gene_name, " - ", df_lab$`Phenotype Description`), NA)
df_lab$LabelGeneCat <- ifelse(df_lab$absZ %in% mostsigGeneCat & df_lab$pvalue < labelcutoff, paste0(df_lab$gene_name, " - ", df_lab$`Phenotype Description`), NA)
df_lab$LabelCat <- ifelse(df_lab$absZ %in% mostsigCat & df_lab$pvalue < labelcutoff, paste0(df_lab$gene_name, " - ", df_lab$`Phenotype Description`), NA)


xBreaks <- unique(df_lab$xlabelpoint)
colors4graphing <- unique(df_lab$color)
colors <- df_lab %>% select(groupnum, color, Category) %>% 
  filter(color %in% colors4graphing) %>% unique() %>% arrange(groupnum)
catColors <- colors$color
names(catColors) = colors$Category
xLabels <- as.ordered(colors$Category)
df_lab$Category <- as.ordered(df_lab$Category)

# Manhattan plot
p_lab_bonferroni <- ggplot(df_lab, aes(x=Position, y = -log10(pvalue), color = Category, label=LabelPheno)) +
  geom_point(size = 1.5, position = "jitter") +
  geom_label_repel(size=3, min.segment.length = 1, max.overlaps = Inf, nudge_y = 15, force=40, force_pull=0.01) +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust = 1),
        legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ylim(0,30) +
  scale_color_manual(values=catColors) +
  scale_x_continuous(labels = xLabels, breaks = xBreaks) +
  geom_hline(yintercept=-log10(cutoff), color = "red")+
  labs(y = "-log10(p-value)")
p_lab_bonferroni
ggsave('SuppFig3B_Manhattan_Lab.tiff', plot=p_lab_bonferroni, units="in", width=6, height=5)
