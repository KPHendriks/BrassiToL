#With this R script we update gene multiple sequence alignments by pruning rogue and outgroup taxa (except Cleome).

#LOAD PACKAGES ----
library(stringr)
library(seqinr)
library(ape)
library(Rmisc)
library(phytools)
library(phangorn)
library(treeio)

#READ IN METADATA ----
metadata<-read.csv("2e.metadata.csv", header=T, stringsAsFactors=FALSE, fileEncoding="latin1")

#DEFINE THE POSSIBLE CONFOUNDING TAXA ----
prune_samples_by_tribe<-metadata[metadata$Tribe %in% c("Anastaticeae", "Biscutelleae", "Megacarpaeeae"),]$Library_ID
prune_samples_by_outgroup<-metadata$Library_ID[grepl("Outgroup_", metadata$Tribe)]
samples_to_prune<-unique(c(prune_samples_by_tribe, prune_samples_by_outgroup))  
#Make sure to include a single outgroup: PAFTOL_019361 (Cleome albescens).
samples_to_prune<-samples_to_prune[!samples_to_prune == "PAFTOL_019361"]

#LOAD LIST OF GENES TO KEEP ----
genes_to_keep<-read.csv(file = "4d.genes_to_keep.txt", header = F)
  
#LOAD GENES ALIGNMENTS AND PRUNE POSSIBLE CONFOUNDING TAXA. ----
for (i in 1:nrow(genes_to_keep)){
  # for (i in 1:10){
  tryCatch({
    #Load gene alignment.
    gene<-strsplit(genes_to_keep[i,], "_")[[1]][1]
    #Read the fasta file.
    fasta<-read.fasta(paste0("results_ASTRALIII_superstrict_and_clocklike_and_pruned/input_gene_alignments/", genes_to_keep[i,]))
    #Prune the pasta file.
    fasta_pruned<-fasta[which(!names(fasta) %in% samples_to_prune)]
    #Save the pruned version of the fasta file.
    write.fasta(sequences = fasta_pruned, 
                names = names(fasta_pruned),
                file = paste0("results_ASTRALIII_superstrict_and_clocklike_and_pruned/input_gene_alignments_pruned/", gene, "_pruned.fasta"))
    #Print progress update to the screen.
    print(paste0("Finished pruning alignment for gene ",gene ," (", i," of in total ", nrow(genes_to_keep)," genes)."))
  })
}
