#!/bin/bash

#Create a directory to collect all iqtree2 gene trees from pipeline.
mkdir results_ASTRALIII_inclusive/
cd results_ASTRALIII_inclusive/

#Copy iqtree2 gene trees.
cp ../results_iqtree_final_inclusive/*.treefile .

#First combine all gene tree files into a single file.
cat *.treefile > gene_trees_combined.tre

#Now remove iqtree2 .treefiles again.
rm *.treefile

#Remove empty lines with just a ";", which means keeping only odd lines and moving these to new file.
##NOTE: this was needed before when I concatenated tree files from pasta which had unexpected line breaks.
#sed 'n; d' gene_trees_combined.tre > gene_trees_combined2.tre
#sed "s/'//g" gene_trees_combined2.tre > gene_trees_combined3.tre

#Now run ASTRALIII.
java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -o species_tree.tre 2> species_tree.log
java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 2 -o species_tree_t2.tre 2> species_tree_t2.log
#java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 16 -o species_tree_t16.tre 2> species_tree_t16.log
#java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 32 -o species_tree_t32.tre 2> species_tree_t32.log
#java -jar ~/astral/ASTRAL/Astral/astral.5.7.8.jar -i gene_trees_combined.tre -t 1 -o species_tree_t1.tre 2> species_tree_t1.log

#Move up to main directory.
cd ../
