#With this R script we calculate genetic distances, substitution saturation, clocklikeness, and several other metrics to assess 
#suitability of genes for inclusion in the phylogenomic (ASTRAL/IQ-TREE) steps to follow.


#LOAD PACKAGES ----
library(stringr)
library(seqinr)
library(ape)
library(Rmisc)
library(phytools)
library(phangorn)
library(treeio)


#LOAD PREVIOUS RESULTS ----
genetic_distances<-readRDS(file = "results_assessment_saturation_and_clocklikeness/results_R_objects/genetic_distances.Rds")


#CALCULATE SUBSTITUTION SATURATION ----

#This analysis was done using the GUI on of PhyloMAd (https://github.com/duchene/phylomad) on an Apple Macbook Pro.
#PhyloMAd simply takes the gene alignments used in the analysis of gene distances above
#(directory "results_assessment_saturation_and_clocklikeness/input_gene_alignments_cleaned/").

#Read in the results from PhyloMAd.
saturation<-read.csv(file = "results_assessment_saturation_and_clocklikeness/results_PhyloMAd/saturation.test.results.csv", header = T)
saturation$genes<-sapply(strsplit(saturation$X, "_"), "[", 1)

#Combine results using merging of tables.
#Note that merging keeps the smallest number of genes from either table, which in this case
#is from saturation; saturation could not be calculated for all genes.
results_assessment_table<-merge(genetic_distances, saturation, by = "genes")
results_assessment_table<-results_assessment_table[, -which(names(results_assessment_table) == "X")]


#CALCULATE clocklikeness ----

#Load two R functions, kindly supplied by David Duchene.
source("results_assessment_saturation_and_clocklikeness/R_functions_David_Duchene/pathnode.R")
source("results_assessment_saturation_and_clocklikeness/R_functions_David_Duchene/rtt.cov.R")


#CACULATE NODE SUPPORT ----

#Get gene tree file names and pull gene names from those.
gene_tree_files<-list.files(path = "results_assessment_saturation_and_clocklikeness/results_iqtree_gene_trees/", pattern = ".treefile", full.names = F)
gene_tree_list<-sapply(str_split(gene_tree_files, "_"), "[", 2)

#Add a column to the dataframe to store the results.
results_assessment_table$gene_tree_mean_node_support<-NA
results_assessment_table$gene_tree_clocklikeness<-NA

#Loop through the genes already in the table and analyse.
for (j in 1:nrow(results_assessment_table)){
  #Read the tree file.
  gene_j<-results_assessment_table$genes[j]
  gene_tree_j<-read.tree(file = paste0("results_assessment_saturation_and_clocklikeness/results_iqtree_gene_trees/result_", gene_j, "_iqtree.treefile"))
  #Calculate mean bootstrap from all nodes and collect in dataframe.
  results_assessment_table$gene_tree_mean_node_support[j]<-mean(as.numeric(gene_tree_j$node.label), na.rm = T)
  #Calculate clocklikeness (def.: "the root-to-tip coefficient of variation") and collect in dataframe.
  results_assessment_table$gene_tree_clocklikeness[j]<-rtt.cov(gene_tree_j)
  #Print update to screen.
  print(paste0("Finished tree ", j, " of a total of ", nrow(results_assessment_table), " gene trees."))
}


#ADD SNP PROP. RESULTS ----

#We first copy in a table from HybPhaser that we can use to add results on mean 
#Collect results from the HybPhaser SNP table saved as R object and take the mean for each gene.
HybPhaser_Table_SNPs<-readRDS("results_HybPhaser_allGenes/output_HybPhaser_allGenes/00_R_objects/inclusive/Table_SNPs.Rds")
genes_mean_SNP_prop<-data.frame(genes=names(rowMeans(as.data.frame(HybPhaser_Table_SNPs), na.rm = T)),
                               mean_SNP_prop=rowMeans(as.data.frame(HybPhaser_Table_SNPs), na.rm = T))

#Combine with previous results.
results_assessment_table<-merge(results_assessment_table, genes_mean_SNP_prop, by = "genes")

#SET RULES TO SUBSET GENES ----

#Suggestions from David Duchene (by email):

# Removing saturated genes will definitely help with resolving deep divergences. 
# It is an intuitive and effective first filter. There are two other filters that are a bit 
# more difficult to use, but which we have studied in detail and found to be very effective across 
# a large number of data sets (Vankan et al. 2022, Syst Biol, 71(2) 490-500).
# 
# The first is simply the mean branch support across branches for each gene
# tree (can be bootstrap or aLRT). Gene trees with high support have lots of signal about 
# deeper branches, and also there is no reason to believe that the gene trees with low supports
# have a particular signal in a direction that could bias the results when removed. You can choose 
# a threshold of 50 or perhaps 75 but only if the first threshold removes very few or no genes.
# 
# The second is the root-to-tip coefficient of variation (aka clocklikeness). 
# A value of 0 indicates an even signal about ancestors for each of your samples. 
# What we find is that genes with very high values have a signal that is restricted to
# only a few samples, so that their overall information about deep branches is poor. 
# The threshold for this metric is not simple to define, but you might want to exclude 
# a 5 or 10% of genes with the highest values. I attach the R function rtt.cov that receives 
# a tree (depends on pathnode, also attached).

#Set required thresholds.
remove_saturated_genes<-"yes" #Set whether to remove genes flagged by PhyloMAd as 'high risk' of saturation.
remove_low_BS_genes<-75 #Set the lower threshold (bootstrap in %) of genes to be removed.
remove_non_clock_genes<-0.5 #Set the proportion of genes with the least clocklikeness to be removed.
remove_high_SNP_prop_genes<-0.02 #Set an upper threshold of genes to be removed (cf. HybPhaser removal of likely paralogs).
remove_short_alignments<-250 #Set minimum number of bp of alignment length.

#Find genes to keep after applying the above thresholds.
genes_non_saturated<-if(remove_saturated_genes=="yes"){results_assessment_table[!(results_assessment_table$Risk_Entropy=='high.risk' | results_assessment_table$Risk_Entropyvar=='high.risk'), ]$genes}else{results_assessment_table$genes}
genes_high_BS<-results_assessment_table[results_assessment_table$gene_tree_mean_node_support>remove_low_BS_genes, ]$genes
genes_clocklike<-results_assessment_table[rank(results_assessment_table$gene_tree_clocklikeness)/length(results_assessment_table$gene_tree_clocklikeness)<(1-remove_non_clock_genes), ]$genes
genes_SNP_prop_low<-results_assessment_table[results_assessment_table$mean_SNP_prop<remove_high_SNP_prop_genes,]$genes
genes_alignments_long<-results_assessment_table[results_assessment_table$gene_length>remove_short_alignments,]$genes

#We only keep the genes that comply with all of the above rules.
#Find the intersection of the five vectors created above.
genes_to_keep<-Reduce(intersect, list(genes_non_saturated, genes_high_BS, genes_clocklike, genes_SNP_prop_low, genes_alignments_long))
#The number of genes to keep is as follows.
length(genes_to_keep)

#Add the results from the above five criteria and the combined vector to the results dataframe.
results_assessment_table$genes_non_saturated<-as.numeric(results_assessment_table$genes %in% genes_non_saturated)
results_assessment_table$genes_high_BS<-as.numeric(results_assessment_table$genes %in% genes_high_BS)
results_assessment_table$genes_clocklike<-as.numeric(results_assessment_table$genes %in% genes_clocklike)
results_assessment_table$genes_SNP_prop_low<-as.numeric(results_assessment_table$genes %in% genes_SNP_prop_low)
results_assessment_table$genes_alignments_long<-as.numeric(results_assessment_table$genes %in% genes_alignments_long)
results_assessment_table$genes_to_keep<-as.numeric(results_assessment_table$genes %in% genes_to_keep)


#SAVE TABLE AS R OBJECT ----
saveRDS(results_assessment_table, file = "results_assessment_saturation_and_clocklikeness/results_R_objects/results_assessment_table.Rds")

#OR LOAD PREVIOUS RESULTS ----
results_assessment_table<-readRDS(file = "results_assessment_saturation_and_clocklikeness/results_R_objects/results_assessment_table.Rds")


#Check if there is much overlap in genes that are (non) saturated and those that have high/low SNP props.  
addmargins(prop.table(table(Non_saturated=results_assessment_table$genes_non_saturated, 
                            SNP_prop_low=results_assessment_table$genes_SNP_prop_low)))
#We find this is not the case, which suggests that these metrics really do score different things.
#(Note that in the case of correlation, we would expect high values on the diagonal and low values outside.)

#EXPORT LISTS WITH GENES AND TREES TO KEEP ----
write.table(paste0("result_", genes_to_keep, "_iqtree.treefile"), file = "4d.trees_to_keep.txt", quote = F, row.names = F, col.names = F)
write.table(paste0(genes_to_keep, "_cleaned.fasta"), file = "4d.genes_to_keep.txt", quote = F, row.names = F, col.names = F)


#Export the assessment table.
# write.csv(results_assessment_table, file = "4e.assessment_saturation_clocklikeness_table.csv")




