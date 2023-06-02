#With this script we subset fasta files based on a list of samples that need to be removed.


#Load packages.
library(ape)
library(seqinr)
library(stringr)

#Read sample name from shell script.
args = commandArgs(trailingOnly=TRUE)
gene<-args[1]

#Read list of samples to be removed.
samples_to_remove<-read.table(file = "2n.list_of_hybrids_and_rogue_taxa.txt")

#Read in the fasta file for this gene.
fasta<-read.fasta(paste0("results_HybPhaser_allGenes/output_HybPhaser_allGenes/03_sequence_lists_superstrict_excluding_hybrids/loci_consensus/", gene, "_consensus.fasta"), seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)

#Subset the fasta file by removing samples from the list.
fasta_subset<-fasta[-which(names(fasta) %in% samples_to_remove[,1])]

#Overwrite original fasta file with the subsetted fasta file.
write.fasta(as.list(fasta_subset),
            names=names(fasta_subset),
            file=paste0("results_HybPhaser_allGenes/output_HybPhaser_allGenes/03_sequence_lists_superstrict_excluding_hybrids/loci_consensus/", gene, "_consensus.fasta"))
