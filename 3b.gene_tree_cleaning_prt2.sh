#!/bin/bash

#This file contains the true analyses script.

#Read in variables.
routine=$1
gene=$2

#Create a directory for each gene to work in and move there.
mkdir -p ${gene}; cd ${gene}
#Remove any files from previous runs in case such a previous run was not finished and cleared.
rm ${gene}*

#Copy the aligned fasta file for this gene from TreeShrink to the current directory.
cp ../results_${routine}/results_treeshrink_${routine}/${gene}_treeshrink_output.fasta .

#Use MAFFT to align.
mafft ${gene}_treeshrink_output.fasta > ${gene}_mafft_final_${routine}.fasta
cp ${gene}_mafft_final_${routine}.fasta ../results_${routine}/results_mafft_final_${routine}/

#Repeat TrimAl (https://github.com/inab/trimal) to trim the MAFFT alignments before creating gene trees.
#trimal -in ${gene}_mafft_final_${routine}.fasta -out ${gene}_trimal_final_${routine}.fasta -gappyout
~/trimal/source/trimal -in ${gene}_mafft_final_${routine}.fasta -out ${gene}_trimal_final_${routine}.fasta -resoverlap 0.75 -seqoverlap 0.90 -gt 0.90

#Save results from TrimAl.
cp ${gene}_trimal_final_${routine}.fasta ../results_${routine}/results_trimal_final_${routine}/

#Mask possible remaining errors from multiple sequence alignments using TAPER (https://github.com/chaoszhang/TAPER).
julia ~/TAPER/correction_multi.jl -m N -a N ${gene}_trimal_final_${routine}.fasta > ${gene}_taper_final_${routine}.fasta

#Save results from TAPER.
cp ${gene}_taper_final_${routine}.fasta ../results_${routine}/results_taper_final_${routine}/

#Run ultrafast bootstrapping with IQ-TREE.
#Chosen settings are taken from Baker et al. (2021) A Comprehensive Phylogenomic Platform for Exploring the Angiosperm Tree of Life.
~/iqtree-2.2.2-Linux/bin/iqtree2 -s ${gene}_taper_final_${routine}.fasta -B 1000 -nm 1000 -m GTR+F+R -nt AUTO --prefix ${gene}_iqtree_final_${routine}

#Save all IQ-TREE output for later use.
cp ${gene}_iqtree_final_${routine}.* ../results_${routine}/results_iqtree_final_${routine}/

#Go one directory up and remove this gene's directory.
cd ..
rm -r ${gene}/
