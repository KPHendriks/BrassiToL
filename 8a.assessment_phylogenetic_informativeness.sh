#!/bin/bash
#SBATCH -J R4S

#[[THE FIRST PART OF THIS SCRIPT WAS WRITTEN TO BE RAN ON A LOCAL MACBOOK MACHINE]]

#This script can be used to assess the phylogenetic informativeness of our gene data, and compare these among bait sets.

mkdir -p results_assessment_phylogenetic_informativeness/
mkdir -p results_assessment_phylogenetic_informativeness/input_species_tree/
mkdir -p results_assessment_phylogenetic_informativeness/input_gene_alignments/
mkdir -p results_assessment_phylogenetic_informativeness/input_species_trees_pruned/
mkdir -p results_assessment_phylogenetic_informativeness/results_rate4site/

#Copy the species tree that serves as the reference.
cp 9b.IQ-TREE_supermatrix_approach.tree results_assessment_phylogenetic_informativeness/input_species_tree/

#Copy the gene alignments that we want to use.
###NOTE UPDATE TO TAKE STRICT DATA!
cp results_ASTRALIII_inclusive/results_taper_final_inclusive/*.fasta results_assessment_phylogenetic_informativeness/input_gene_alignments/

#Make an overview of all the genes included.
find results_assessment_phylogenetic_informativeness/input_gene_alignments/ -type f -path 'results_assessment_phylogenetic_informativeness/input_gene_alignments/*.fasta' -exec basename {} ';' | sort -u -o 8c.genelist.txt
sed -i.bak 's/_taper_final_inclusive.fasta//g' 8c.genelist.txt

#Use the following R script to create versions of the species tree that are subsetted to each alignment, i.e. all samples not in the phylogeny need to be pruned.
###NOTE still need to update the R script to read the 8c.genelist from this terminal variable!!
Rscript 8b.assessment_phylogenetic_informativeness_prt1.R 8c.genelist.txt

#[[THE NEXT PART OF THIS SCRIPT WAS WRITTEN TO BE RAN ON A HPC, SUCH AS THE NATURALIS HIGH MEM CLUSTER]]

#Analyse sequence alignments in parallel.
parallel -j 47 bash 8d.assessment_phylogenetic_informativeness_by_gene.sh {} :::: 8c.genelist.txt

#[[THE NEXT PART OF THIS SCRIPT WAS WRITTEN TO BE RAN ON A LOCAL MACBOOK MACHINE]]

#Finally run the following R script to collect and further analyse the results and create an informative graph.
Rscript 8e.assessment_phylogenetic_informativeness_prt2.R
