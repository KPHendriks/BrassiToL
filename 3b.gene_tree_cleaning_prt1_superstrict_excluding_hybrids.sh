#!/bin/bash

#This file contains the true analyses script.

#Read in variables.
routine=$1
gene=$2

#Create a directory for each gene to work in and move there.
mkdir -p ${gene}; cd ${gene}
#Remove any files from previous runs in case such a previous run was not finished and cleared.
rm ${gene}*

#Copy the fasta file for this gene to the current directory.
cp ../results_HybPhaser_allGenes/output_HybPhaser_allGenes/03_sequence_lists_${routine}/loci_consensus/${gene}*.fasta .

#Subset the fasta file to only those samples we want to include in our further analysis.
#Note that the subsetting was done prior to running HybPhaser, but here any accidentally included samples can still be pruned from the dataset.
#~/seqtk/seqtk subseq ${gene}_consensus.fasta ../3d.samplelist.txt > ${gene}_subset.fasta
mv ${gene}_consensus.fasta ${gene}_subset.fasta

#Use MAFFT to align.
mafft ${gene}_subset.fasta > ${gene}_mafft_preliminary_${routine}.fasta
cp ${gene}_mafft_preliminary_${routine}.fasta ../results_${routine}/results_mafft_preliminary_${routine}/

#Run TrimAl (https://github.com/inab/trimal) to trim the MAFFT alignments before creating gene trees.
#trimal -in ${gene}_mafft_preliminary_strict.fasta -out ${gene}_trimal_preliminary_strict.fasta -gappyout
 ~/trimal/source/trimal -in ${gene}_mafft_preliminary_${routine}.fasta -out ${gene}_trimal_preliminary_${routine}.fasta -resoverlap 0.75 -seqoverlap 0.90 -gt 0.90

#Save results from TrimAl.
cp ${gene}_trimal_preliminary_${routine}.fasta ../results_${routine}/results_trimal_preliminary_${routine}/

#Mask possible remaining errors from multiple sequence alignments using TAPER (https://github.com/chaoszhang/TAPER).
julia ~/TAPER/correction_multi.jl -m N -a N ${gene}_trimal_preliminary_${routine}.fasta > ${gene}_taper_preliminary_${routine}.fasta

#Save results from TAPER.
cp ${gene}_taper_preliminary_${routine}.fasta ../results_${routine}/results_taper_preliminary_${routine}/

#Create preliminary ML tree in IQ-TREE.
#Chosen settings are taken from Baker et al. (2021) A Comprehensive Phylogenomic Platform for Exploring the Angiosperm Tree of Life.
~/iqtree-2.2.2-Linux/bin/iqtree2 -s ${gene}_taper_preliminary_${routine}.fasta -nm 1000 -m GTR+F+R -nt AUTO --prefix ${gene}_iqtree_preliminary_${routine}

#Save all IQ-TREE output for later use.
cp ${gene}_iqtree_preliminary_${routine}.* ../results_${routine}/results_iqtree_preliminary_${routine}/

#Go one directory up and remove this gene's directory.
cd ..
rm -r ${gene}/
