This folder contains files used to constrain our analyses and to present/
describe our results. 


The files are as follows:
  
"Dennis2021_LabGWAS_TraitList_v2_Annotated.csv" - The list of traits available 
for replication analysis of blood traits within the BioVU large-scale biobank, 
mapped onto equivalent traits from the UKBB lab/continuous data analysis. 

"GenesTestedForEachJTITissue_All.csv" - Lists all genes tested for each GTEx 
tissue analyzed, including both gene name and Ensembl gene ID. This serves as 
a foundation for determining the total number of tests performed for multiple 
hypothesis testing correction.
  
"ProgrammedCellDeathGeneList_Pathways.csv" - The list of genes used in the 
analysis annotated with their associated pathway, subpathway (if such a 
subdivision is relevant), and corresponding protein name. There are several 
instances where there's considerable divergence between the two (for instance,
the gene "BCL2L2" encodes the protein "BCL-W").

"ProgrammedCellDeathGeneList_Pathways_Extra.csv" - Select genes involved in 
ferroptosis, cuproptosis, and parthanatos cell death pathways incorporated into 
follow-up analysis.

"SNPsTestedForEachGeneJTI.csv" - Our JTI TWAS models select different numbers of
SNPs to model for a given gene across different tissues. This file describes the 
number of SNPs modeled for each gene in each tissue, plus provides the number of
tissues for which a given gene was modeled (and therefore tested within our 
analyses). This is essential for multiple hypothesis testing correction and 
relevant to understanding the scope of our analysis. 

"UKBBv3_AllTraits2Category_FINAL_v2_plusCaseControl.csv" - The list of traits 
which formed the basis of our analysis of the UKBB (version 3, as analyzed and 
reported by the Neale Lab  "http://www.nealelab.is/uk-biobank"). Beyond 
information relevant to the UKBB including phenotype codes and case/control 
counts for the version 3 GWAS studies, this file also assigns each trait to a 
category used in the construction of Manhattan plots.

"categorycolors2phenotypes.csv" - This file assigns a color (R colors, described 
here "https://r-charts.com/colors/") and a position to traits for constructing 
beautiful Manhattan plots.
