#With this script we set prepare the files that serve as input for the tool rate4site, with which we will
#calculate mutation rates from gene alignments, constrained on the species tree.

#Load packages.
library(ape)
library(phytools)

#Load the species tree that we want to use to make the assessment of phylogenetic informativeness for each gene.
#We chose to work with the ML tree based on a supermatrix of strict genes.
species_tree<-read.tree(file = "results_assessment_phylogenetic_informativeness/input_species_tree/9b.IQ-TREE_supermatrix_approach.tree")
#Root the tree to the outgroup.
species_tree<-root(species_tree, outgroup = "PAFTOL_019361")
#Make the tree ultrametric.
species_tree <- chronos(species_tree)
#It might be needed to reset the class of the tree object to "phylo".
class(species_tree)<-"phylo"
#Remove bootstrap values from the tree.
species_tree$node.label<-NULL

genelist<-read.csv("8c.genelist.txt", header = F)

#For each gene alignment, create a species tree pruned to the entries in the gene alignment.
for (i in 1:nrow(genelist)){
  gene<-genelist[i,]
  #Load the alignment fasta file.
  fasta_file<-read.FASTA(file = paste0("results_assessment_phylogenetic_informativeness/input_gene_alignments/", gene, "_taper_final_inclusive.fasta"))
  #Read the sample names from this alignment.
  fasta_file_names<-names(fasta_file)
  #Create a copy of the species tree.
  species_tree_pruned<-drop.tip(species_tree, species_tree$tip.label[!species_tree$tip.label %in% fasta_file_names])
  #Save this copy of the pruned species tree for use in rate4site.
  write.tree(species_tree_pruned, file = paste0("results_assessment_phylogenetic_informativeness/input_species_trees_pruned/",gene,"_pruned_species_tree.newick"))
}

