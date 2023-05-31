#!/bin/bash

#Create several directories to store results.
#In case a directory already exists we'll simply remove it first, to make sure new and old results will not get mixed up.
rm -r results_mafft_preliminary_strict; mkdir results_mafft_preliminary_strict
rm -r results_trimal_preliminary_strict; mkdir results_trimal_preliminary_strict
rm -r results_taper_preliminary_strict; mkdir results_taper_preliminary_strict
rm -r results_iqtree_preliminary_strict; mkdir results_iqtree_preliminary_strict
rm -r results_treeshrink_strict; mkdir results_treeshrink_strict
rm -r results_mafft_final_strict; mkdir results_mafft_final_strict
rm -r results_trimal_final_strict; mkdir results_trimal_final_strict
rm -r results_taper_final_strict; mkdir results_taper_final_strict
rm -r results_iqtree_final_strict; mkdir results_iqtree_final_strict

#Create a genelist_strict from HybPhaser R script 1c output.
ls results_HybPhaser_allGenes/output_HybPhaser_allGenes/03_sequence_lists_strict/loci_consensus/ > genelist_consensus_strict.txt
sed 's/_consensus.fasta//g' genelist_consensus_strict.txt > 3c.genelist_strict.txt
rm genelist_consensus_strict.txt
#Note that at this point, the above file can be replaced to run a user-provided subset of genes from HybPhaser.


#For the 'strict' routine, start by adding PAFTOL sequence data.
#This is data for which we had no access to raw sequence reads and which we therefore could not map ourselves.
#Instead we'll have to make do with the (max 353) target sequence reads that could be downloaded from the PAFTOL Tree of Life Explorer (https://treeoflife.kew.org/).

# Move the new consensus sequences from HybPhaser and the PAFTOL sequences to new directories, which makes it easier to concatenate new with old results.
mkdir -p sequences_one/ sequences_two/ sequences_three/ sequences_four/ sequences_five/ sequences_six/ sequences_seven/ sequences_combined/
cp results_HybPhaser_allGenes/output_HybPhaser_allGenes/03_sequence_lists_strict/loci_consensus/*.fasta sequences_one/
cp PAFTOL_sequence_data/CRNC/*.FNA sequences_two/
cp PAFTOL_sequence_data/HYZL/*.FNA sequences_three/
cp PAFTOL_sequence_data/MYZV/*.FNA sequences_four/
cp PAFTOL_sequence_data/SRR6441723/*.FNA sequences_five/
cp PAFTOL_sequence_data/ERR4210213/*.FNA sequences_six/
cp PAFTOL_sequence_data/ERR2789774/*.FNA sequences_seven/
#Rename the PAFTOL_sequence_data files to match the new sequence files.
cd sequences_two; for f in *.FNA; do mv -- "$f" "${f%.FNA}_consensus.fasta"; done; cd ..
cd sequences_three; for f in *.FNA; do mv -- "$f" "${f%.FNA}_consensus.fasta"; done; cd ..
cd sequences_four; for f in *.FNA; do mv -- "$f" "${f%.FNA}_consensus.fasta"; done; cd ..
cd sequences_five; for f in *.FNA; do mv -- "$f" "${f%.FNA}_consensus.fasta"; done; cd ..
cd sequences_six; for f in *.FNA; do mv -- "$f" "${f%.FNA}_consensus.fasta"; done; cd ..
cd sequences_seven; for f in *.FNA; do mv -- "$f" "${f%.FNA}_consensus.fasta"; done; cd ..
#Find the intersection of files from all these folders and save in a list.
find sequences_one/ sequences_two/ sequences_three/ sequences_four/ sequences_five/ sequences_six/ sequences_seven/ -type f -path '*/*.fasta' -exec basename {} ';' | sort -u -o file.list
# Then concatenate fasta files with same names, file by file. We pull only genes from the PAFTOL samples that we also have from HybPhaser.
#(I.e. any genes that HybPhaser removed, are now also ignored from the PAFTOL data, even if PAFTOL retrieved these sequences.)
while read -r name; do find sequences_two sequences_three sequences_four sequences_five sequences_six sequences_seven sequences_one -type f -path "*/$name" -exec cat {} + > "sequences_combined/$name"; done < file.list
#Copy the results back to the HybPhaser output folder.
cp sequences_combined/*.fasta results_HybPhaser_allGenes/output_HybPhaser_allGenes/03_sequence_lists_strict/loci_consensus/
#Remove temporary directories and files.
rm -r sequences_two/ sequences_three/ sequences_four/ sequences_five/ sequences_six/ sequences_seven/ sequences_combined/ sequences_one/ file.list

#Can use the following line to rename all sequences inside a fasta file to the name of the sample, which needs to be done once for data from PAFTOL sequences.
# for i in *.FNA; do sed "1s/.*/>SRR6441723/" $i > temp; mv temp $i; done

#Start running a script in parallel that for each gene creates an alignment using MAFFT and improved the alignment using.
#TrimAl, then creates a gene tree using IQ-TREE.
parallel -j 47 bash 3b.gene_tree_cleaning_prt1.sh {} :::: 3c.genelist_strict.txt


#List the genes for which IQ-TREE did actually find a tree (meaning that a fasta file was also successfully exported from TrimAl and TAPER);
#Create an updated genelist_strict file to work with next.
ls results_iqtree_preliminary_strict/*.treefile > 3c.genelist_strict_succesfull_temp1.txt
sed 's\results_iqtree_preliminary_strict/\\g' 3c.genelist_strict_succesfull_temp1.txt > 3c.genelist_strict_succesfull_temp2.txt
sed 's\_iqtree_preliminary_strict.treefile\\g' 3c.genelist_strict_succesfull_temp2.txt > 3c.genelist_strict_succesfull.txt
rm 3c.genelist_strict_succesfull_temp*.txt

#Use TreeShrink (https://github.com/uym2/TreeShrink) to detect and remove outliers based on long branches in gene trees.
#TreeShrink takes as input all the gene trees together and weighs these to detect outliers.
#The output consists of pruned gene trees and alignments.

#Create a temporary directory to work in.
rm -r temp_treeshrink_strict/
mkdir -p temp_treeshrink_strict/
cd temp_treeshrink_strict/
#Create an input directory structure as required by TreeShrink.
#Loop through the genelist_strict and create a directory for each gene with related gene tree and gene alignment.
while read gene
do
  mkdir -p ${gene}
  cp ../results_iqtree_preliminary_strict/${gene}_iqtree_preliminary_strict.treefile ${gene}/
  mv ${gene}/${gene}_iqtree_preliminary_strict.treefile ${gene}/input.tree
  cp ../results_taper_preliminary_strict/${gene}_taper_preliminary_strict.fasta ${gene}/
  mv ${gene}/${gene}_taper_preliminary_strict.fasta ${gene}/input.fasta
done < ../3c.genelist_strict_succesfull.txt

#Now we can run TreeShrink on the created set of gene trees and alignments.
python ~/TreeShrink/run_treeshrink.py -i . -t input.tree -a input.fasta -f > treeshrinklog.txt

#Rename the output from TreeShrink such that file names relate to the relevant gene.
while read gene
do
  mv ${gene}/input.tree ${gene}/${gene}_treeshrink_input.tree
  mv ${gene}/input.fasta ${gene}/${gene}_treeshrink_input.fasta
  mv ${gene}/output.tree ${gene}/${gene}_treeshrink_output.tree
  mv ${gene}/output.fasta ${gene}/${gene}_treeshrink_output.fasta
  mv ${gene}/output.txt ${gene}/${gene}_treeshrink_output.txt
done < ../3c.genelist_strict_succesfull.txt

#Save the output from TreeShrink for later use.
cp */*_treeshrink_*put.* ../results_treeshrink_strict/
cp treeshrinklog.txt ../results_treeshrink_strict/
mv output_summary.txt treeshrink_output_summary.txt
cp treeshrink_output_summary.txt ../results_treeshrink_strict/

#Move up one directory and remove the temporary directory for TreeShrink.
cd ..
rm -r temp_treeshrink_strict/

#Run another script in parallel that for each gene creates an updates alignment using MAFFT
#and then creates an updated gene tree using IQ-TREE.
parallel -j 47 bash 3b.gene_tree_cleaning_prt2.sh {} :::: 3c.genelist_strict_succesfull.txt
