#!/bin/bash

#This script is used to calculate several important gene aligment/tree metrics that need to be assessed to decide
#what genes should be removed from further analysis.

#[[The following part of th script was written to be ran locally on Apple Macbook pro.]]

#Create several directories to work in.
mkdir -p results_assessment_saturation_and_clocklikeliness/input_gene_alignments_cleaned/
mkdir -p results_assessment_saturation_and_clocklikeliness/results_R_objects/
mkdir -p results_assessment_saturation_and_clocklikeliness/results_iqtree_gene_trees/


#Calculate genetic distances from gene alignments using the following R script.
Rscript 4b.assessment_genetic_distances.R

#[[The following part of th script was written to be ran on a HPC such as the Naturalis HighMem cluster.]]

#In the previous R script, samples with short reads were removed from the gene alignments.
#Now infer new bootstrapped ML trees with IQ-TREE.
#Chosen settings are taken from Baker et al. (2021) A Comprehensive Phylogenomic Platform for Exploring the Angiosperm Tree of Life.

#Create a file with all gene alignments to run in parallel below.
ls results_assessment_saturation_and_clocklikeliness/input_gene_alignments_cleaned/ > 4c.gene_alignments_cleaned_list.txt
sed -i.bak 's/_cleaned.fasta//g' 4c.gene_alignments_cleaned_list.txt

#Run iq-tree in parallel on a computing cluster.
parallel -j 47 '~/iqtree-2.1.3-Linux/bin/iqtree2 -s results_assessment_saturation_and_clocklikeliness/input_gene_alignments_cleaned/{}_cleaned.fasta -B 1000 -nm 1000 -m GTR+F+R -nt AUTO --prefix results_assessment_saturation_and_clocklikeliness/results_iqtree_gene_trees/result_{}_iqtree' :::: 4c.gene_alignments_cleaned_list.txt

#[[The following part of th script was written to be ran locally on Apple Macbook Pro.]]

#Calculate substitution saturation from alignments.
#This analysis was done using the GUI on of PhyloMAd (https://github.com/duchene/phylomad) on an Apple Macbook Pro.
#PhyloMAd simply takes the gene alignments used in the analysis of gene distances above
#(directory "results_assessment_saturation_and_clocklikeliness/input_gene_alignments_cleaned/").
#The following R script will read in the results from the PhyloMAd output table.


#Calculate gene tree mean node support (from bootstrap values) and clocklikeness using the following R script.
#Make sure to set several important selection criteria in this R script!
Rscript 4c.assessment_saturation_and_clocklikeness.R

#Using only the genes that passed the selection criteria, combine gene trees and run ASTRAL-III.

#Create a directory to collect all iqtree2 gene trees from pipeline.
mkdir results_ASTRALIII_superstrict_and_clocklike/

#Copy iqtree2 gene trees; use the text file written by the previous R script.
cat 4d.trees_to_keep.txt | while read FILENAME; \
do find results_assessment_saturation_and_clocklikeliness/results_iqtree_gene_trees/ -name "$FILENAME" -type f -print -exec cp '{}' results_ASTRALIII_superstrict_and_clocklike \; ; done;

#Move to the newly created directory.
cd results_ASTRALIII_superstrict_and_clocklike/

#First combine all gene tree files into a single file.
cat *.treefile > gene_trees_combined.tre

#Now remove iqtree2 .treefiles again.
rm *.treefile

#Now run ASTRALIII.
java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -o species_tree.tre 2> species_tree.log
java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 2 -o species_tree_t2.tre 2> species_tree_t2.log
#java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 16 -o species_tree_t16.tre 2> species_tree_t16.log
#java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 32 -o species_tree_t32.tre 2> species_tree_t32.log
#java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 1 -o species_tree_t1.tre 2> species_tree_t1.log

#Move up to main directory.
cd ../


#In addition, we take this subset of genes but now remove all all possibly confounding taxa
#(all outgroups except sample PAFTOL_019361 (Cleome albescens) and all ‘rogue taxa’, viz. all members of tribes Anastaticeae, Biscutelleae, and Megacarpaeeae.
#After removing samples from the alignments, we run MAFFT once more to possibly improve the alignments.

#Create some directories to store results in.
mkdir results_ASTRALIII_superstrict_and_clocklike_and_pruned/
mkdir results_ASTRALIII_superstrict_and_clocklike_and_pruned/input_gene_alignments/
mkdir results_ASTRALIII_superstrict_and_clocklike_and_pruned/input_gene_alignments_pruned/
mkdir results_ASTRALIII_superstrict_and_clocklike_and_pruned/results_mafft/
mkdir results_ASTRALIII_superstrict_and_clocklike_and_pruned/results_iqtree_gene_trees/

#Copy gene alignments from the genes to keep.
cat 4d.genes_to_keep.txt | while read FILENAME; \
do find results_assessment_saturation_and_clocklikeliness/input_gene_alignments_cleaned/ -name "$FILENAME" -type f -print -exec cp '{}' results_ASTRALIII_superstrict_and_clocklike_and_pruned/input_gene_alignments/ \; ; done;

#Use the following R script to remove all possible confounding taxa.
Rscript 4f.update_gene_alignments.R

#Use MAFFT to update the alignments.
cp 4d.genes_to_keep.txt 4d.genes_to_keep_short.txt
sed -i.bak 's\_cleaned.fasta\\g' 4d.genes_to_keep_short.txt
parallel -j 45 'mafft --quiet results_ASTRALIII_superstrict_and_clocklike_and_pruned/input_gene_alignments_pruned/{}_pruned.fasta > results_ASTRALIII_superstrict_and_clocklike_and_pruned/results_mafft/{}_mafft.fasta' :::: 4d.genes_to_keep_short.txt

#Run IQ-TREE to create gene trees. We use the same settings as for the ASTRAL-III gene tree for species tree analysis before.
parallel -j 45 '~/iqtree-2.1.3-Linux/bin/iqtree2 -s results_ASTRALIII_superstrict_and_clocklike_and_pruned/results_mafft/{}_mafft.fasta -B 1000 -nm 1000 -m GTR+F+R -nt AUTO --prefix results_ASTRALIII_superstrict_and_clocklike_and_pruned/results_iqtree_gene_trees/result_{}_iqtree' :::: 4d.genes_to_keep_short.txt

#Combine all gene tree files into a single file for ASTRAL-III to work with.
cat results_ASTRALIII_superstrict_and_clocklike_and_pruned/results_iqtree_gene_trees/*.treefile > results_ASTRALIII_superstrict_and_clocklike_and_pruned/gene_trees_combined.tre

#Move to the ASTRAL-III folder.
cd results_ASTRALIII_superstrict_and_clocklike_and_pruned/

#Now run ASTRALIII.
java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -o species_tree.tre 2> species_tree.log
java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 2 -o species_tree_t2.tre 2> species_tree_t2.log
#java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 16 -o species_tree_t16.tre 2> species_tree_t16.log
#java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 32 -o species_tree_t32.tre 2> species_tree_t32.log
#java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 1 -o species_tree_t1.tre 2> species_tree_t1.log

#Move up to main directory.
cd ../
