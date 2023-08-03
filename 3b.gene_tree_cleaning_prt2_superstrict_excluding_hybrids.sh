#!/bin/bash

#This file contains the true analyses script.

gene=$1

#Create a directory for each gene to work in and move there.
mkdir -p ${gene}; cd ${gene}
#Remove any files from previous runs in case such a previous run was not finished and cleared.
rm ${gene}*

#Copy the aligned fasta file for this gene from TreeShrink to the current directory.
cp ../results_treeshrink_superstrict_by_tribe/${gene}_treeshrink_output.fasta .

#Use MAFFT to align.
mafft ${gene}_treeshrink_output.fasta > ${gene}_mafft_final_superstrict_by_tribe.fasta
cp ${gene}_mafft_final_superstrict_by_tribe.fasta ../results_mafft_final_superstrict_by_tribe/

#Repeat TrimAl (https://github.com/inab/trimal) to trim the MAFFT alignments before creating gene trees.
#trimal -in ${gene}_mafft_final_superstrict_by_tribe.fasta -out ${gene}_trimal_final_superstrict_by_tribe.fasta -gappyout
trimal -in ${gene}_mafft_final_superstrict_by_tribe.fasta -out ${gene}_trimal_final_superstrict_by_tribe.fasta -resoverlap 0.75 -seqoverlap 0.90 -gt 0.90

#Save results from TrimAl.
cp ${gene}_trimal_final_superstrict_by_tribe.fasta ../results_trimal_final_superstrict_by_tribe/

#Mask possible remaining errors from multiple sequence alignments using TAPER (https://github.com/chaoszhang/TAPER).
~/julia-1.7.0/bin/julia ~/TAPER/correction_multi.jl -m N -a N ${gene}_trimal_final_superstrict_by_tribe.fasta > ${gene}_taper_final_superstrict_by_tribe.fasta

#Save results from TAPER.
cp ${gene}_taper_final_superstrict_by_tribe.fasta ../results_taper_final_superstrict_by_tribe/

#Run ultrafast bootstrapping with IQ-TREE.
#Chosen settings are taken from Baker et al. (2021) A Comprehensive Phylogenomic Platform for Exploring the Angiosperm Tree of Life.
~/iqtree-2.1.3-Linux/bin/iqtree2 -s ${gene}_taper_final_superstrict_by_tribe.fasta -B 1000 -nm 1000 -m GTR+F+R -nt AUTO --prefix ${gene}_iqtree_final_superstrict_by_tribe

#Save all IQ-TREE output for later use.
cp ${gene}_iqtree_final_superstrict_by_tribe.* ../results_iqtree_final_superstrict_by_tribe/

#Go one directory up and remove this gene's directory.
cd ..
rm -r ${gene}/
