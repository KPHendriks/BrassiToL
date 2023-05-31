#With this script we prune gene trees and gene alignments in order for the tip labels/samples to correspond to the species tree
#for which we want to calculate concordance factors in IQ-TREE.

#LOAD PACKAGES ----
library(ape)
library(treeio)
library(seqinr)


#GET SPECIES TREE TIP LABELS ----
tree_IQ_TREE_supermatrix_calibrated<-read.beast(file = "results_treePL_calibration/output_treePL_dating/treePL_calibrated_phylogeny_summary.tre")
species_tree_tip_labels<-tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label


#PRUNE GENE TREES ----

#Load gene trees.
gene_trees_input<-read.tree(file = "results_treePL_calibration/input_tree_files/gene_trees_combined.tree")

#Create an empty tree object for storage of pruned gene trees.
gene_trees_input_pruned<-list()
class(gene_trees_input_pruned)<-"multiPhylo"

#Loop through gene trees and prune all tip labels not in the species tree.
for (i in 1:length(gene_trees_input)){
  #Get gene tree.
  gene_tree<-gene_trees_input[[i]]
  #Drop tips.
  gene_tree_pruned<-drop.tip(gene_tree, which(!gene_tree$tip.label %in% species_tree_tip_labels))
  #Store pruned gene tree in overall gene trees object.
  gene_trees_input_pruned[[i]]<- gene_tree_pruned
}


#SAVE GENE TREES ----
write.tree(gene_trees_input_pruned,
           file = "results_treePL_calibration/input_tree_files/gene_trees_combined_pruned.tree")


#PRUNE GENE ALIGNMENTS ----

#Get gene list.
gene_list<-read.csv("results_IQ-TREE_supermatrix/output_concatenated_sequences/list.genes.to.keep", header = F)
gene_list<-sapply(strsplit(gene_list[,1], "/"), "[", 3)
gene_list<-sapply(strsplit(gene_list, "_"), "[", 1)

#Calculate genetic distances from alignments. ----
for (i in 1:length(gene_list)){
  # for (i in 1:10){
  tryCatch({
    #Load gene alignment.
    gene<-gene_list[i]
    fasta<-read.alignment(paste0("results_treePL_calibration/input_multiple_sequence_alignments/", gene_list[i], "_taper_final_inclusive.fasta"), format="fasta")
    #Get sequence index for samples to be removed.
    keep_seqs<-which(fasta$nam %in% species_tree_tip_labels)
    #Remove sequences for these samples (or actually: keep samples that need not be removed).
    fasta_cleaned<-fasta$seq[keep_seqs]
    #Save the updated fasta file to work with; note the subsetting of sequence names needs to follow the subsetting of genes,
    #or otherwise read data will be connected to the wrong samples!
    write.fasta(fasta_cleaned, file.out = paste0("results_treePL_calibration/input_multiple_sequence_alignments_pruned/" ,gene, ".fasta"), names = fasta$nam[keep_seqs])
    #Print a progress update to the screen.
    print(paste0("Finished calculations for gene ",i ," of in total ", length(gene_list)," genes (",round(100*i/length(gene_list),1),"%)."))
  })
}


