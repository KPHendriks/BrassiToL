#!/bin/bash

#This file contains the true analyses script.

gene=$1

#Create a directory for each gene to work in and move there.
mkdir -p results_assessment_phylogenetic_informativeness/${gene}/; cd results_assessment_phylogenetic_informativeness/${gene}/
#Remove any files from previous runs in case such a previous run was not finished and cleared.
rm ${gene}*

#Run rate4site on this gene alignment.
#[[IF RUNNING LOCALLY ON MY MAC:]]
#~/Software/rate4site/rate4site.3.2.source/sourceMar09/rate4site -i ml -mn -s ../input_gene_alignments/${gene}_taper_final_superstrict.fasta
#[[IF RUNNING ON NATURALIS HIGHMEM COMPUTING CLUSTER:]]
~/rate4site.3.2.source/sourceMar09/rate4site -i ml -mn -s ../input_gene_alignments/${gene}_taper_final_inclusive.fasta -t ../input_species_trees_pruned/${gene}_pruned_species_tree.newick

#Rename and copy the result to the results directory.
mv r4sOrig.res ${gene}_result_rate4site.res
cp ${gene}_result_rate4site.res ../results_rate4site/

#Go one directory up and remove this directory.
cd ..
rm -r ${gene}/
