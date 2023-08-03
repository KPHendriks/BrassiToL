#This script will load the supermatrix species tree and prune any tip labels not in the
#most clocklike gene trees used for calibration in treePL.

library(ape)

#Load gene alignments.
clockline_gene_alignments_concatenated<-read.FASTA("results_treePL_calibration/output_concatenated_sequences/allseqs.fas")

#Find samples from alignemnt.
samples<-names(clockline_gene_alignments_concatenated)

#Load species tree.
species_tree_to_fix_topology<-read.tree("9b.IQ-TREE_supermatrix_approach.tree")

#Prune this species tree.
tips_to_drop<-species_tree_to_fix_topology$tip.label[! species_tree_to_fix_topology$tip.label %in% samples]
species_tree_to_fix_topology_pruned<-drop.tip(species_tree_to_fix_topology, tips_to_drop)

#Set all branch lengths to a random value of 1.
species_tree_to_fix_topology_pruned$edge.length<-rep(1, length(species_tree_to_fix_topology_pruned$edge.length))

#Write pruned species tree.
write.tree(species_tree_to_fix_topology_pruned,
           file = "results_treePL_calibration/input_supermatrix_species_tree/species_tree_to_fix_topology_pruned.tree")
