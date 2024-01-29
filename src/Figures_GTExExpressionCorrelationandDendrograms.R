### Correlation Expression and Modules
 
# Load packages -------
library(tidyverse)
library(corrplot)
library(reshape2)
library(ggrepel)
 
# Load data
df <- read_csv('data/GTExv8_TPM_ProgrammedCellDeathGenes.csv')
 
# Define colors for graphing specific pathways
gcolor <- c("#EFD4A5","#6993A7","#E59A98")
darkcolor <- c("#B88547","#274B77","#974643")
pathwaycolors <- data.frame(color=c("#EFD4A5","#6993A7","#E59A98"), Pathway=c("APOPTOSIS","NECROPTOSIS","PYROPTOSIS"))
 
# Define pathways
pathways <- read_csv('metadata/ProgrammedCellDeathGeneList_Pathways.csv')
pathways[40,3] <- "DFNA5" # GSDME is DFNA5 in GTEx dataset

# Define directory for saving figures
figurepath="data/figures/"

# Calculate median TPM for each tissue
all <- df %>% pivot_longer(cols=c(3:46), names_to="gene_name", values_to="TPM") 
all <- all %>% group_by(gene_name, Tissue) %>% 
  mutate(medianTPM = median(TPM)) %>% 
  select(Tissue, gene_name, medianTPM) %>% 
  unique()
   
# Select only tissues in GTEx v8 analysis
tissuelist <- as.vector(unique(df$Tissue), mode="list") # create a list
tissuelist <- tissuelist[-c(7,24,25,31,35)] # remove extra tissues
tissuevec <- as.vector(tissuelist, mode="character")
all <- all %>% dplyr::filter(Tissue %in% tissuevec) %>% select(gene_name, Tissue, medianTPM) %>% unique()
 
# Annotate with pathway, subpathway information
all <- merge(all, pathways, by.x="gene_name", by.y="Gene Name", all.x=T)
all$SubPathway <- factor(all$SubPathway, 
                         levels = c("PROSURVIVAL", "PROAPOPTOTIC BH3-ONLY", 
                                    "PROAPOPTOTIC EFFECTOR", "INITIATOR CASPASE", 
                                    "APOPTOSOME",
                                    "EXECUTIONER CASPASE", "CASPASE ANTAGONIST",
                                    "NECROPTOSIS", "PYROPTOSIS"))
 
# Clean up tissue labeling
all$Tissue <- str_replace_all(all$Tissue, "_", " ")
 
 
# Graph expression values (median TPM) in a boxplot -------
p.medianTPM <- ggplot(all, aes(x=gene_name, y=medianTPM, fill=Pathway)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, face="italic", hjust = 1), strip.text.x = element_blank()) + 
  scale_fill_manual(values=gcolor) +
  xlab("Gene Name") + ylab("Median TPM") +
  ylim(0,500) +
  facet_grid(cols=vars(SubPathway), scales="free", space="free") 
p.medianTPM
ggsave(path = figurepath, filename ='GTExExpression_MedianTPM_AllTissues.tiff',
       plot=p.medianTPM, device="tiff", units = "in", width = 8.5, height = 4)
 
 
# Clustering of tissues based on median TPM values from all genes -------
alltissue <- all %>% pivot_wider(id_cols=gene_name, names_from=Tissue, values_from=medianTPM)
alltissue <- column_to_rownames(alltissue, var="gene_name")
 
# Compute Pearson's correlation
mat_tissue <- cor(alltissue)
 
# Graph using corrplot package; clustering using Ward's D; 4 modules
png(paste0(figurepath,"Corrplot_GTEx_onTissue_geneTPM.png"), width = 900, height = 700)
corrplot(mat_tissue, method="square", order="hclust", 
         hclust.method="ward.D", addrect = 4, tl.col="black", rect.lwd=3, cl.cex=1, rect.col="firebrick3")
dev.off()
 
# Pull correlation values for each tissue-tissue comparison
TissueCorr.res <- as.data.frame(mat_tissue)
TissueCorr.res <- rownames_to_column(TissueCorr.res, var = "Tissue1")
TissueCorr.res <- pivot_longer(TissueCorr.res, cols = 2:50, names_to ="Tissue2")
TissueCorr.res$Comparison <- paste0(TissueCorr.res$Tissue1, sep="-", TissueCorr.res$Tissue2) # create tissue-tissue comparison column
keycomparisons <- data.frame(TissueCorr.res, stringsAsFactors = FALSE) %>% 
  mutate(key = paste0(pmin(Tissue1, Tissue2), sep = "-", pmax(Tissue1, Tissue2))) %>% 
  distinct(key) # Pull only unique tissue-tissue comparisons without considering order
TissueCorr.res <- filter(TissueCorr.res, Comparison %in% keycomparisons$key)
write_csv(TissueCorr.res, 'data/GeneGeneCorrelations_GTExModules.csv')
 
 
# Module-specific analysis setup -------
tissuemods <- read_csv('results/GTExTissueModules.csv')
mod1 <- filter(tissuemods, Module == "1")
mod2 <- filter(tissuemods, Module == "2")
mod3 <- filter(tissuemods, Module == "3")
mod4 <- filter(tissuemods, Module == "4")

# Module-specific expression values -------
# Gene expression for mod 1 -- bonus
all.mod1 <- filter(all, Tissue %in% mod1$Tissue)
p.mod1.medianTPM <- ggplot(all.mod1, aes(x=gene_name, y=medianTPM, fill=Pathway)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, face="italic", hjust = 1), strip.text.x = element_blank()) + 
  scale_fill_manual(values=gcolor) +
  xlab("Gene Name") + ylab("Median TPM") +
  ylim(0,500) +
  facet_grid(cols=vars(SubPathway), scales="free", space="free") 
p.mod1.medianTPM
ggsave(path = figurepath, filename ='GTExExpression_MedianTPM_Mod1.tiff',
       plot=p.mod1.medianTPM, device="tiff", units = "in", width = 8.5, height = 4)
 
# Gene expression for mod 2
all.mod2 <- filter(all, Tissue %in% mod2$Tissue)
p.mod2.medianTPM <- ggplot(all.mod2, aes(x=gene_name, y=medianTPM, fill=Pathway)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, face="italic", hjust = 1), strip.text.x = element_blank()) + 
  scale_fill_manual(values=gcolor) +
  xlab("Gene Name") + ylab("Median TPM") +
  ylim(0,500) +
  facet_grid(cols=vars(SubPathway), scales="free", space="free") 
p.mod2.medianTPM
ggsave(path = figurepath, filename ='GTExExpression_MedianTPM_Mod2.tiff',
       plot=p.mod2.medianTPM, device="tiff", units = "in", width = 8.5, height = 4)
 
# Gene expression for mod 3
all.mod3 <- filter(all, Tissue %in% mod3$Tissue)
p.mod3.medianTPM <- ggplot(all.mod3, aes(x=gene_name, y=medianTPM, fill=Pathway)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, face="italic", hjust = 1), strip.text.x = element_blank()) + 
  scale_fill_manual(values=gcolor) +
  xlab("Gene Name") + ylab("Median TPM") +
  ylim(0,500) +
  facet_grid(cols=vars(SubPathway), scales="free", space="free") 
p.mod3.medianTPM
ggsave(path = figurepath, filename ='GTExExpression_MedianTPM_Mod3.tiff',
       plot=p.mod3.medianTPM, device="tiff", units = "in", width = 8.5, height = 4)
rm(p.mod1.medianTPM, p.mod2.medianTPM, p.mod3.medianTPM)
 
 
# Module-specific correlation clustering (dendrograms) -------
# Gene/gene correlations across tissues in brain module
mod1 <- filter(tissuemods, Module == "1")
all_mod1 <- all %>% filter(Tissue %in% mod1$Tissue)
all_mod1 <- dcast(all_mod1, gene_name~Tissue, value.var = "medianTPM")
all_mod1 <- column_to_rownames(all_mod1, var="gene_name")
all_mod1 <- t(all_mod1)
mat_mod1 <- cor(all_mod1)
png(paste0(figurepath,"Corrplot_GTEx_PCD_geneTPM_mod1.png"), width = 600, height = 600)
corrplot(mat_mod1, method="square", order="hclust", hclust.method="ward.D", addrect = 5, tl.col = "black")
dev.off() 
 
# Gene/gene correlations by module 2 
mod2 <- filter(tissuemods, Module=="2")
all_mod2 <- all %>% filter(Tissue %in% mod2$Tissue)
all_mod2 <- dcast(all_mod2, gene_name~Tissue, value.var = "medianTPM")
all_mod2 <- column_to_rownames(all_mod2, var="gene_name")
all_mod2 <- t(all_mod2)
mat_mod2 <- cor(all_mod2)
 
# Gene/gene correlations by module 3 
mod3 <- filter(tissuemods, Module=="3")
all_mod3 <- all %>% filter(Tissue %in% mod3$Tissue)
all_mod3 <- dcast(all_mod3, gene_name~Tissue, value.var = "medianTPM")
all_mod3 <- column_to_rownames(all_mod3, var="gene_name")
all_mod3 <- t(all_mod3)
mat_mod3 <- cor(all_mod3)
 
# Gene/gene correlations by module 4
mod4 <- filter(tissuemods, Module=="4")
all_mod4 <- all %>% filter(Tissue %in% mod4$Tissue)
all_mod4 <- dcast(all_mod4, gene_name~Tissue, value.var = "medianTPM")
all_mod4 <- column_to_rownames(all_mod4, var="gene_name")
all_mod4 <- t(all_mod4)
mat_mod4 <- cor(all_mod4)
 
 
# Dendrogram construction -------
pathways <- as.data.frame(pathways)
i=0
colLab<<-function(n){
  if(is.leaf(n)){
    a=attributes(n)
    ligne=match(attributes(n)$label,pathways[,3])
    pathway=pathways[ligne,c(1)];
    if(pathway=="APOPTOSIS"){col_pathway="#B88547"};if(pathway=="NECROPTOSIS"){col_pathway="#274B77"};if(pathway=="PYROPTOSIS"){col_pathway="#974643"}
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=20,lab.col=col_pathway, col=col_pathway, lab.font=1,lab.cex=1))
  }
  return(n)
}
sim.brain <- hclust(dist(mat_mod1))
genes_brain <- cutree(sim.brain, k=5)
sim.brain <- as.dendrogram(sim.brain)
p.sim.brain <- dendrapply(sim.brain, colLab)
plot(p.sim.brain , main="Module 1 Gene Expression Relationships")
abline(h=4.4, col="red")
legend("topleft", 
       legend = c("Apoptosis" , "Necrotposis" , "Pyroptosis"), 
       col = c("#EFD4A5", "#6993A7", "#E59A98"), 
       pch = c(20,20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0, 0)) #750 x 400 size
dev.off()
 
# rapidly dividing tissues module 2
sim.m2 <- hclust(dist(mat_mod2))
genes_mod2 <- cutree(sim.m2, k=5)
sim.m2 <- as.dendrogram(sim.m2)
p.sim.m2 <- dendrapply(sim.m2, colLab)
plot(p.sim.m2 , main="Module 2 Gene Expression Relationships")
abline(h=3, col="red")
legend("topright", 
       legend = c("Apoptosis" , "Necrotposis" , "Pyroptosis"), 
       col = c("#EFD4A5", "#6993A7", "#E59A98"), 
       pch = c(20,20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0, 0))
dev.off()
 
# module 3
sim.m3 <- hclust(dist(mat_mod3))
genes_mod3 <- cutree(sim.m3, k=5)
sim.m3 <- as.dendrogram(sim.m3)
p.sim.m3 <- dendrapply(sim.m3, colLab)
plot(p.sim.m3 , main="Module 3 Gene Expression Relationships")
abline(h=5.4, col="red")
legend("topright", 
       legend = c("Apoptosis" , "Necrotposis" , "Pyroptosis"), 
       col = c("#EFD4A5", "#6993A7", "#E59A98"), 
       pch = c(20,20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0, 0))
dev.off()
 
