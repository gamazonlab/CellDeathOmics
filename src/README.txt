This folder contains R scripts used to perform multiple hypothesis testing
analysis, enrichment analyses, and figure generation. 

The files are as follows:

"EnrichmentOfTraitAssociations.R" - Performs Fisher's Exact test for enrichment 
of PCD genes across all traits. These results are shown in Fig. 5C.

"Figures_GTExExpressionCorrelationandDendrograms.R" - This generates correlation 
plots and dendrograms for the GTEx-derived expression data deposited in 
'CellDeathOmics/data'. This script generates the graphs shown in Fig. 2B-F.

"Figures_Heatmaps_BioVUreplication.R" - This generates the tiled heatplot 
represented in Fig. 6 from results from the replication of blood traits from 
BioVU cross-referenced with significant associations identified from the 
primary UKBB blood trait analysis.

"Figures_Heatmaps_UKBB.R" - This generates the tiled heatplots represented in 
Figs. 4B and 5B with median effect sizes for clinical traits and lab-derived 
traits from the UKBB S-PrediXcan analysis. 

"Figures_ManhattanPlots.R" - This customizable script generates the Manhattan plots 
shown in Figs. 4A and 5A for clinical traits and lab-derived traits from the 
UKBB S-PrediXcan analysis. 

"HeritabilityOfCuratedUKBBTraits.R" - This script generates the boxplot shown 
in Fig. 3B by comparing LDSC-derived heritabilty estimates for UKBB traits 
curated for our analyses vs. those excluded from our analyses. 

"MultipleHypothesisTestingCorrectionandFormatting.R" - This outlines our 
approach to multiple hypothesis testing corrections (Benjamini-Hochberg) for 
UKBB S-PrediXcan analyses. 
