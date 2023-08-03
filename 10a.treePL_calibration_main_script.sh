#!/bin/bash

#Create some directories to work in.
mkdir -p results_treePL_calibration/
mkdir -p results_treePL_calibration/data_gene_trees_clocklike/
mkdir -p results_treePL_calibration/output_concatenated_sequences/
mkdir -p results_treePL_calibration/output_iqtree/
mkdir -p results_treePL_calibration/output_treePL_prime/
mkdir -p results_treePL_calibration/output_treePL_cross_validation/
mkdir -p results_treePL_calibration/output_treePL_dating/
mkdir -p results_treePL_calibration/input_tree_files/
mkdir -p results_treePL_calibration/input_multiple_sequence_alignments_pruned/
mkdir -p results_treePL_calibration/output_iqtree_condordance_factors/

###DELETE FROM HERE#####
#Start by copying the supermatrix approach IQ-TREE phylogeny (plus bootstrapped phylogenies) to work with in treePL.
#We will use this tree as our topology.
#Note that this tree is already rooted by setting an outgroup in IQ-TREE before.
cp 9b.IQ-TREE_supermatrix_approach.tree results_treePL_calibration/input_supermatrix_species_tree/
cp results_IQ-TREE_supermatrix/output_iqtree/iqtree_bootstrap_bootstrapped_trees_from_fixed_topology.tree results_treePL_calibration/input_supermatrix_species_tree/
###DELETE TO HERE#####

#We now run IQ-TREE once again, now only using the 20 most 'clocklike' genes.
#The phylogeny from the full supermatrix approach will be used to fix the topology in IQ-TREE.

#The 20 most 'clocklike' genes were found in script 4 ("4a.assessment_saturation_and_clocklikeliness.sh"), and are as follows.
#Copy the alignments for these genes.
while read gene
do
  cp results_ASTRALIII_inclusive/results_taper_final_inclusive/${gene}_taper_final_inclusive.fasta results_treePL_calibration/data_gene_trees_clocklike/
done < 10b.genes_clocklike.txt

#Create a list of fasta files.
ls results_treePL_calibration/data_gene_trees_clocklike/*.fasta > 10b.genes_clocklike_filenames.txt

#Use the cat sequences tool (https://github.com/ChrisCreevey/catsequences) to concatenate the multiple sequence alignments into a single supermatrix.
~/catsequences/catsequences 10b.genes_clocklike_filenames.txt

#Move the results to a directory with a sensible name.
mv allseqs.fas results_treePL_calibration/output_concatenated_sequences/
mv allseqs.partitions.txt results_treePL_calibration/output_concatenated_sequences/

#Use an R script to remove any tip labels from the species tree used to contrain in case these are not in the clocklike dataset.
Rscript 10c.treePL_prune_tip_labels_species_treeR.R

#Run IQ-TREE on this clocklike dataset, but use topology from full dataset, as in script 9 ("9a.IQ-TREE_supermatrix_approach.sh").
~/iqtree-2.1.3-Linux/bin/iqtree2 -s results_treePL_calibration/output_concatenated_sequences/allseqs.fas -o "S1321" -g results_treePL_calibration/input_supermatrix_species_tree/species_tree_to_fix_topology_pruned.tree -m GTR+F+R --prefix results_treePL_calibration/output_iqtree/ -bb 1000 -wbtl -nt 46 #-wbtl makes sure to write bootstrap trees including branch lengths.
####!!!! TEMP THE ABOVE LINE RESULTS IN 1,000 EXACTLY THE SAME BOOTSTRAP TREES, WITH EXACTLY THE SAME BRANCH LENGTHS EACH...!
~/iqtree-2.1.3-Linux/bin/iqtree2 -s results_treePL_calibration/output_concatenated_sequences/allseqs.fas -o "S1321" -m GTR+F+R --prefix results_treePL_calibration/output_iqtree/ -bb 1000 -wbtl -nt 46 #-wbtl makes sure to write bootstrap trees including branch lengths.
##NEXT ATTEMPT: the older slower parametric bootstrapped
~/iqtree-2.1.3-Linux/bin/iqtree2 -s results_treePL_calibration/output_concatenated_sequences/allseqs.fas -o "S1321" -g results_treePL_calibration/input_supermatrix_species_tree/species_tree_to_fix_topology_pruned.tree -m GTR+F+R --prefix results_treePL_calibration/output_iqtree_parametric_bootstrap/ -b 100 -wbtl -nt 46


#IQ-TREE saves some of the files we need as hidden files, so make these visible first.
find results_treePL_calibration/output_iqtree/ -type f -name '.*' -execdir sh -c 'mv -i "$0" "./${0#./.}"' {} \;

#And rename the files we need with a sensible name.
mv results_treePL_calibration/output_iqtree/iqtree results_treePL_calibration/output_iqtree/iqtree_report.txt
mv results_treePL_calibration/output_iqtree/treefile results_treePL_calibration/output_iqtree/iqtree_ML_tree.tree
mv results_treePL_calibration/output_iqtree/splits.nex results_treePL_calibration/output_iqtree/iqtree_bootstrap_support_values.nex
mv results_treePL_calibration/output_iqtree/contree results_treePL_calibration/output_iqtree/iqtree_bootstrap_concensus_tree.tree
mv results_treePL_calibration/output_iqtree/ufboot results_treePL_calibration/output_iqtree/iqtree_bootstrap_trees.tree
mv results_treePL_calibration/output_iqtree/log results_treePL_calibration/output_iqtree/iqtree_bootstrap_logfile.log


#Run treePL first to prime the analysis and save output from screen to a log file for later review.
#[Note that this script assumes the use of treePL on the Naturalis HighMem computing cluster, which requires sudo to run.]
sudo treePL 10e.treePL_config_prime > results_treePL_calibration/output_treePL_prime/treePL_prime.log

#Run treePL for cross validation, using the "best optimization parameters" suggested by treePL prime.
sudo treePL 10e.treePL_config_cross_validation > results_treePL_calibration/output_treePL_cross_validation/treePL_cross_validation1.log

#Run treePL to date the tree using the smoothing value from cross validation.
#As input we now use the 1,000 bootstrapped trees from IQ-TREE, such that we can calculate confidence intervals (CIs) on the calibration later.
sudo treePL 10e.treePL_config_dating > results_treePL_calibration/output_treePL_dating/treePL_dating.log

#Now turn to TreeAnnotator and load the 1,000 treePL calibrated bootstrapped trees and use the input IQ-TREE as target tree type > user target tree.
#These are the files to use:
#results_treePL_calibration/output_treePL_dating/treePL_calibrated_phylogeny.tre
#results_treePL_calibration/output_iqtree/iqtree_ML_tree.tree
#In TreeAnnotator, save the result as results_treePL_calibration/output_treePL_dating/treePL_calibrated_phylogeny_summary.tre

#Make a copy of this result for publication.
cp results_treePL_calibration/output_treePL_dating/treePL_calibrated_phylogeny_summary.tre 10f.treePL_supermatrix_approach_calibrated.tree

#And copy the calibrated tree to the main directory.
cp results_treePL_calibration/output_treePL_dating/treePL_calibrated_phylogeny.tre 10f.treePL_calibrated_phylogeny.tre

#Now we run IQ-TREE once more to calculate concordance factors for this specific tree.
#Note that this calibrated species tree is in topology the same as results_IQ-TREE_supermatrix/output_iqtree/iqtree_ML_tree.tree, but lacks
#three tip labels from samples not present in the data from the 20 most clocklike genes.
#Therefore, we need to update gene alignments and gene trees.
#We simply take the same genes as used in the previous calculation of concordance factors (script 9a.IQ-TREE_supermatrix_approach.sh), which already
#for genes with SNP proportions <0.02 only.

#Copy the gene alignments and gene trees from the previous script.
cp -r results_IQ-TREE_supermatrix/input_multiple_sequence_alignments results_treePL_calibration/
cp -r results_IQ-TREE_supermatrix/input_tree_files/gene_trees_combined.tree results_treePL_calibration/input_tree_files/

#Run the following R script to prune samples from gene trees and gene alignments.
Rscript 10g.prune_gene_trees_and_alignments.R

#Start IQ-TREE to run calculation of concordance factors from species tree + gene trees + gene alignments.
#We followed suggestions from http://www.iqtree.org/doc/Concordance-Factor.
~/iqtree-2.1.3-Linux/bin/iqtree2 -t results_treePL_calibration/output_iqtree/iqtree_ML_tree.tree --gcf results_treePL_calibration/input_tree_files/gene_trees_combined_pruned.tree -p results_treePL_calibration/input_multiple_sequence_alignments_pruned/ --scf 1000 --prefix results_treePL_calibration/output_iqtree_condordance_factors/ -nt 47

#Like before, IQ-TREE writes out hidden files. Let's rename these with sensible names and make them visible.
mv results_treePL_calibration/output_iqtree_condordance_factors/.cf.tree results_treePL_calibration/output_iqtree_condordance_factors/iqtree_ML_tree_with_concordance_factors.tree
mv results_treePL_calibration/output_iqtree_condordance_factors/.cf.tree.nex results_treePL_calibration/output_iqtree_condordance_factors/iqtree_ML_tree_with_concordance_factors.nex
mv results_treePL_calibration/output_iqtree_condordance_factors/.cf.branch results_treePL_calibration/output_iqtree_condordance_factors/iqtree_ML_tree_with_branch_IDs.tree
mv results_treePL_calibration/output_iqtree_condordance_factors/.cf.stat results_treePL_calibration/output_iqtree_condordance_factors/iqtree_ML_tree_concordance_factors_per_branch.stat

#IQ-TREE writes out hidden files, so let's make all remaining hidden files visible.
find results_treePL_calibration/output_iqtree_condordance_factors/ -type f -name '.*' -execdir sh -c 'mv -i "$0" "./${0#./.}"' {} \;
