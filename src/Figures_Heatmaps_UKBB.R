### Heatmaps for S-PrediXcan Results
#   Updated 2024JAN19

# Load packages -------
library(tidyverse)
library(ggrepel)
library(corrplot)
library(RColorBrewer)

# Load and format data -------
df.all <- read_csv('results/UKBBv3_SPrediXcanResults_PCDGenes_FINAL.csv')
phenos2color <- read_csv('metadata/categorycolors2phenotypes.csv')

# Annotate data with trait category and color information
df.all <- merge(df.all, phenos2color, all.x=T)

# Factor data for preferred order of category, subpathways
df.all$SubPathway <- factor(df.all$SubPathway, 
                            levels = c("PROSURVIVAL", "PROAPOPTOTIC BH3-ONLY", 
                                       "PROAPOPTOTIC EFFECTOR", "INITIATOR CASPASE", 
                                       "APOPTOSOME", 
                                       "EXECUTIONER CASPASE", "CASPASE ANTAGONIST",
                                       "NECROPTOSIS", "PYROPTOSIS"))
df.all$Category <- factor(df.all$Category, 
                          levels = c("Blood/immune", "Circulatory", "Congenital", 
                                     "Dermatologic", "Digestive", "Ear", "Endocrine/metabolic",
                                     "Eye", "Genitourinary", "Infectious disease", 
                                     "Mental/behavioral", "Musculoskeletal", "Neoplasms", 
                                     "Nervous", "Respiratory", "Other", "White blood cells", "Red blood cells", 
                                     "Platelets", "Urine"))


# Subset for nominally significant traits
df.ehr <- filter(df.all, Derivation == "ICD10/FinnGen" & FDR_all < 0.25)
whittletraits4graphing <- c("High light scatter reticulocyte percentage", 
                        "Haematocrit percentage", "Lymphocyte percentage", 
                        "Monocyte percentage", "Neutrophil percentage",
                        "Reticulocyte percentage", "Mean sphered cell volume")
df.cont <- filter(df.all, Derivation == "Lab/Continuous" & !`Phenotype Description` %in% whittletraits4graphing)
df.cont.sig <- filter(df.cont, FDR_all < 0.01)


# Setup for graphing
gcolor <- c("#EFD4A5","#6993A7","#E59A98")
darkcolor <- c("#B88547","#274B77","#974643")
pathwaycolors <- data.frame(color=c("#EFD4A5","#6993A7","#E59A98"), 
                            Pathway=c("APOPTOSIS","NECROPTOSIS","PYROPTOSIS"))


# Heatmap of EHR traits with FDR<0.25; but for all all gene/trait associations -------
library(reshape2)

df.ehr.all <- filter(df.all, Derivation=="ICD10/FinnGen")
df.ehr.all$FDR25 <- ifelse(df.ehr.all$FDR_all<0.25, 1, 0)
df.ehr.med.all <-  df.ehr.all %>% 
  group_by(gene_name, `Phenotype Description`) %>% 
  mutate(`Median Beta` = median(effect_size)) %>% 
  group_by(gene_name, `Phenotype Description`) %>% slice_min(FDR_all) %>%
  select(gene_name, `Phenotype Description`, Category, Pathway, SubPathway, `Median Beta`, FDR25, FDR_all) %>% 
  unique()
sig.ehr.traits <- df.ehr.all %>% filter(FDR_all<0.25) %>% select(`Phenotype Description`) %>% unique()
sig.ehr.traits <- sig.ehr.traits$`Phenotype Description`

p.ehr.trait.all.2 <- df.ehr.med.all %>% 
  filter(`Phenotype Description` %in% sig.ehr.traits) %>%
  select(gene_name, `Phenotype Description`, `Median Beta`, SubPathway, Category, FDR25)
p.ehr.trait.all.2b = subset(p.ehr.trait.all.2, FDR25==1)

p.ehr.trait.all.heatplot <- ggplot(p.ehr.trait.all.2, aes(x = gene_name, y = `Phenotype Description`, group="FDR25")) +
  geom_tile(aes(fill=`Median Beta`), color="black", size=0.2) + # add color and size parameters to geom_tile()
  geom_text(data = p.ehr.trait.all.2b, aes(x = gene_name, y = `Phenotype Description`, label = "*", group="FDR25", size = 2, vjust = 0.5, hjust = 0.5)) +
  xlab("Gene Name") + ylab("Phenotype Description") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_bw() + # remove color, size, and fill parameters from theme_bw()
  theme(strip.text.x=element_blank(), strip.text.y=element_text(angle = 0), axis.text.x=element_text(angle=90, hjust=1, face = "italic")) + # removed strip.text.y = element_blank(),
  facet_grid(Category ~ SubPathway, scales = "free", space = "free")
p.ehr.trait.all.heatplot
ggsave('figures/Fig4B_Heatplot_EHR.tiff', 
       plot = p.ehr.trait.all.heatplot, device = "tiff", 
       dpi = 200, width = 14, height = 7, units= "in", limitsize = TRUE)


# Heatmap of continuous traits with FDR<0.01; but for all all gene/trait associations -------
df.cont.all <- filter(df.all, Derivation=="Lab/Continuous" & !`Phenotype Description` %in% whittletraits4graphing)
df.cont.all$FDR01 <- ifelse(df.cont.all$FDR_all<0.01, 1, 0)
df.cont.med.all <-  df.cont.all %>% 
  group_by(gene_name, `Phenotype Description`) %>% 
  mutate(`Median Beta` = median(effect_size)) %>% 
  group_by(gene_name, `Phenotype Description`) %>% slice_min(FDR_all) %>%
  select(gene_name, `Phenotype Description`, Category, Pathway, SubPathway, `Median Beta`, FDR01, FDR_all) %>% 
  unique()
min.median.beta <- min(df.cont.med.all$`Median Beta`)
max.median.beta <- max(df.cont.med.all$`Median Beta`)
sig.cont.traits <- df.cont.all %>% filter(FDR_all<0.01) %>% select(`Phenotype Description`) %>% unique()
sig.cont.traits <- sig.cont.traits$`Phenotype Description`

p.cont.trait.all.2 <- df.cont.med.all %>% 
  filter(`Phenotype Description` %in% sig.cont.traits) %>%
  select(gene_name, `Phenotype Description`, `Median Beta`, SubPathway, Category, FDR01)
p.cont.trait.all.2b = subset(p.cont.trait.all.2, FDR01==1)

p.cont.trait.all.heatplot <- ggplot(p.cont.trait.all.2, aes(x = gene_name, y = `Phenotype Description`, group="FDR01")) +
  geom_tile(aes(fill=`Median Beta`), color="black", size=0.2) +
  geom_text(data = p.cont.trait.all.2b, aes(x = gene_name, y = `Phenotype Description`, label = "*", group="FDR01", size = 2, vjust = 0.5, hjust = 0.5)) +
  xlab("Gene Name") + ylab("Phenotype Description") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_bw() +
  theme(strip.text.x = element_blank(), strip.text.y=element_text(angle=0), axis.text.x=element_text(angle=90, hjust=1, face = "italic")) +
  facet_grid(Category ~ SubPathway, scales = "free", space = "free")
p.cont.trait.all.heatplot
ggsave('figures/Fig5B_Heatplot_Lab.tiff', 
       plot = p.cont.trait.all.heatplot, device = "tiff", 
       dpi = 200, width = 14, height = 5, units= "in", limitsize = TRUE)
