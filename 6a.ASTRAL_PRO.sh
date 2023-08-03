#!/bin/bash

# When running on Stampede2 cluster, detail the following:
#   *** Serial Job on Normal Queue ***
#SBATCH	-A TG-TRA170019	                # Account to charge for this run
#SBATCH -J A-PRO         								# Job name
#SBATCH -o A-PRO.o%j     								# Name of stdout output file
#SBATCH -e A-PRO.e%j   									# Name of stderr error file
#SBATCH -p normal          							# Queue (partition) name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=68
#SBATCH -t 48:00:00        							# Run time (hh:mm:ss)
#SBATCH --mail-user=kasper.hendriks@naturalis.nl
#SBATCH --mail-type=none    						# Send no emails

#This script retrieves paralog sequences from HybPiper output and uses ASTRAL-PRO to infer a species tree.

#[[THE FIRST PART OF THE SCRIPT SHOULD BE RAN ON A CLUSTER THAT HOLDS ALL MAPPED HYBPIPER DATA IN TARBAL FILES, E.G. IN OUR CASE THE TACC STAMPEDE2 CLUSTER]]

#Load needed versions of modules.
module load intel/17.0.4
module load gnuparallel
#module load biocontainers/0.1.0 	#With this module, various other tools/modules can be loaded
#module load bwa/0.7.16a
#module load blast/2.6.0
#module load samtools/1.10
#module load python2
module list
pwd
date

#Create a directory to store results later.
#This will overwrite any previous directory with the same name to avoid the risk of merging unrelated data!
rm -r results_ASTRAL_PRO/; mkdir -p results_ASTRAL_PRO/
rm -r results_ASTRAL_PRO/results_HybPiper_paralog_sequences/; mkdir -p results_ASTRAL_PRO/results_HybPiper_paralog_sequences/

#Create a temporary working directory.
rm -r paralog_retrieval/; mkdir -p paralog_retrieval/

#Untar HybPiper mapped data for the selected samples.
parallel -j 68 bash 6b.ASTRAL_PRO_copy_samples.sh {} :::: 6c.ASTRAL_PRO_namelist.txt

#Move to the tmp directory.
cd paralog_retrieval/

#Retrieve paralogs for the samples and genes of our interest.
parallel -j 68 "python ../../HybPiper/paralog_retriever.py ../6c.ASTRAL_PRO_namelist.txt {} > {}.paralogs.fasta" :::: ../6d.ASTRAL_PRO_genelist.txt 2> ../6e.ASTRAL_PRO_paralog_counts.txt

#Save the paralog fasta files to our results directory.
cp *.fasta ../results_ASTRAL_PRO/results_HybPiper_paralog_sequences/

#Remove the temporary directory.
cd ..
rm -r paralog_retrieval/


#[[THE BELOW PART OF THE SCRIPT CAN BE RAN ON NATURALIS HIGHMEM CLUSTER]]
#Create some directories to store intermediate results.
rm -r results_ASTRAL_PRO/results_mafft_preliminary/; mkdir -p results_ASTRAL_PRO/results_mafft_preliminary/
rm -r results_ASTRAL_PRO/results_trimal_preliminary/; mkdir -p results_ASTRAL_PRO/results_trimal_preliminary/
#rm -r results_ASTRAL_PRO/results_trimal_preliminary_renamed/; mkdir -p results_ASTRAL_PRO/results_trimal_preliminary_renamed/
rm -r results_ASTRAL_PRO/results_iqtree_gene_trees_preliminary/; mkdir -p results_ASTRAL_PRO/results_iqtree_gene_trees_preliminary/
rm -r results_ASTRAL_PRO/results_mafft_final/; mkdir -p results_ASTRAL_PRO/results_mafft_final/
rm -r results_ASTRAL_PRO/results_trimal_final/; mkdir -p results_ASTRAL_PRO/results_trimal_final/
rm -r results_ASTRAL_PRO/results_iqtree_gene_trees_final/; mkdir -p results_ASTRAL_PRO/results_iqtree_gene_trees_final/
rm -r results_ASTRAL_PRO/results_treeshrink/; mkdir -p results_ASTRAL_PRO/results_treeshrink/

#Add PAFTOL A353 sequence data from a set of samples for which we have no raw data available.
# Move the new consensus sequences from HybPhaser and the PAFTOL sequences to new directories, which makes it easier to concatenate new with old results.
mkdir -p sequences_one/ sequences_two/ sequences_three/ sequences_four/ sequences_five/ sequences_six/ sequences_seven/ sequences_combined/
cp results_ASTRAL_PRO/results_HybPiper_paralog_sequences/*paralogs.fasta sequences_one/
cp PAFTOL_sequence_data/CRNC/*.FNA sequences_two/
cp PAFTOL_sequence_data/HYZL/*.FNA sequences_three/
cp PAFTOL_sequence_data/MYZV/*.FNA sequences_four/
cp PAFTOL_sequence_data/SRR6441723/*.FNA sequences_five/
cp PAFTOL_sequence_data/ERR4210213/*.FNA sequences_six/
cp PAFTOL_sequence_data/ERR2789774/*.FNA sequences_seven/
#Rename the PAFTOL_sequence_data files to match the new sequence files.
cd sequences_two; for f in *.FNA; do mv -- "$f" "${f%.FNA}.paralogs.fasta"; done; cd ..
cd sequences_three; for f in *.FNA; do mv -- "$f" "${f%.FNA}.paralogs.fasta"; done; cd ..
cd sequences_four; for f in *.FNA; do mv -- "$f" "${f%.FNA}.paralogs.fasta"; done; cd ..
cd sequences_five; for f in *.FNA; do mv -- "$f" "${f%.FNA}.paralogs.fasta"; done; cd ..
cd sequences_six; for f in *.FNA; do mv -- "$f" "${f%.FNA}.paralogs.fasta"; done; cd ..
cd sequences_seven; for f in *.FNA; do mv -- "$f" "${f%.FNA}.paralogs.fasta"; done; cd ..
#Find the intersection of files from all these folders and save in a list.
find sequences_one/ sequences_two/ sequences_three/ sequences_four/ sequences_five/ sequences_six/ sequences_seven/ -type f -path '*/*.fasta' -exec basename {} ';' | sort -u -o file.list
# Then concatenate fasta files with same names, file by file. We pull only genes from the PAFTOL samples that we also have from HybPhaser.
#(I.e. any genes that HybPhaser removed, are now also ignored from the PAFTOL data, even if PAFTOL retrieved these sequences.)
while read -r name; do find sequences_two sequences_three sequences_four sequences_five sequences_six sequences_seven sequences_one -type f -path "*/$name" -exec cat {} + > "sequences_combined/$name"; done < file.list
#Copy the results back to the HybPhaser output folder.
cp sequences_combined/*.fasta results_ASTRAL_PRO/results_HybPiper_paralog_sequences/
#Remove temporary directories and files.
rm -r sequences_two/ sequences_three/ sequences_four/ sequences_five/ sequences_six/ sequences_seven/ sequences_combined/ sequences_one/ file.list

#Remove the B764 copy of each of the 36 duplicate genes (i.e. genes present in both B764 and A353 bait sets) from the genelist file that we'll work with next.
sed -i.bak '/1G03190/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/1G07010/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/1G08490/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/1G17760/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/1G31780/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/1G43860/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/1G49540/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/1G55880/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/1G75330/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/2G21280/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/2G31955/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/2G39090/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/3G02660/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/3G06950/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/3G08650/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/3G12080/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/3G20780/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/3G48500/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/3G51050/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/3G53760/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/4G02790/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/4G29380/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/4G29830/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/4G33030/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G06830/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G10920/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G13680/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G14760/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G17290/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G18570/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G46180/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G48470/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G48600/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G61770/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G64050/d' 6d.ASTRAL_PRO_genelist.txt
sed -i.bak '/5G67530/d' 6d.ASTRAL_PRO_genelist.txt


#Use MAFFT to align all fasta files.
#First make sure that the gene list file does not contain special characters after export from Excel, so just in case, remove any.
sed $'s/[^[:print:]\t]//g' 6d.ASTRAL_PRO_genelist.txt > 6d.ASTRAL_PRO_genelist_any_special_characters_removed.txt
parallel -j 45 'mafft --quiet results_ASTRAL_PRO/results_HybPiper_paralog_sequences/{}.paralogs.fasta > results_ASTRAL_PRO/results_mafft_preliminary/{}.paralogs.mafft.fasta' :::: 6d.ASTRAL_PRO_genelist_any_special_characters_removed.txt

#Use trimal to trim mafft sequence alignments.
parallel -j 45 'trimal -in results_ASTRAL_PRO/results_mafft_preliminary/{}.paralogs.mafft.fasta -out results_ASTRAL_PRO/results_trimal_preliminary/{}.paralogs.trimal.fasta -resoverlap 0.75 -seqoverlap 0.90 -gt 0.90' :::: 6d.ASTRAL_PRO_genelist_any_special_characters_removed.txt
#Note that trimal updates sequence names in fasta files just as we need them, with anything after the space removed.

#Update names in the fasta files to shorter names without spaces.
#This is not needed if trimal was just used.
#parallel -j 45 "cat results_ASTRAL_PRO/results_trimal_preliminary/{}.paralogs.trimal.fasta | awk -F '. ' '{print $1}' > results_ASTRAL_PRO/results_trimal_preliminary_renamed/{}.paralogs.trimal.renamed.fasta" :::: 6d.ASTRAL_PRO_genelist_any_special_characters_removed.txt

#Run IQ-TREE to create gene trees. We use the same settings as for the ASTRAL-III gene tree for species tree analysis before.
parallel -j 45 '~/iqtree-2.1.3-Linux/bin/iqtree2 -s results_ASTRAL_PRO/results_trimal_preliminary/{}.paralogs.trimal.fasta -B 1000 -nm 1000 -m GTR+F+R -nt AUTO --prefix results_ASTRAL_PRO/results_iqtree_gene_trees_preliminary/result_{}_iqtree' :::: 6d.ASTRAL_PRO_genelist_any_special_characters_removed.txt

#A353 gene 6128 has shown extreme numbers of paralogs (1817 for a single sample!); this we cannot trust, so we remove that gene tree first.
rm results_ASTRAL_PRO/results_iqtree_gene_trees_preliminary/result_6128_iqtree.treefile




#######START ADDED###########


#List the genes for which IQ-TREE did actually find a tree (meaning that a fasta file was also successfully exported from TrimAl and TAPER);
#Create an updated genelist file to work with next.
ls results_ASTRAL_PRO/results_iqtree_gene_trees_preliminary/*.treefile > 6d.ASTRAL_PRO_genelist_succesfull_temp1.txt
sed 's\results_ASTRAL_PRO/results_iqtree_gene_trees_preliminary/result_\\g' 6d.ASTRAL_PRO_genelist_succesfull_temp1.txt > 6d.ASTRAL_PRO_genelist_succesfull_temp2.txt
sed 's\_iqtree.treefile\\g' 6d.ASTRAL_PRO_genelist_succesfull_temp2.txt > 6d.ASTRAL_PRO_genelist_succesfull.txt
rm 6d.ASTRAL_PRO_genelist_succesfull_temp*.txt

#Again, make sure to remove gene 6128, which we want to exclude from further analyses because of extremely high paralog numbers.
sed -i.bak 's\6128\\g' 6d.ASTRAL_PRO_genelist_succesfull.txt
#Remove empty line.
sed -i.bak '/^$/d' 6d.ASTRAL_PRO_genelist_succesfull.txt

#Use TreeShrink (https://github.com/uym2/TreeShrink) to detect and remove outliers based on long branches in gene trees.
#TreeShrink takes as input all the gene trees together and weighs these to detect outliers.
#The output consists of pruned gene trees and alignments.

#Create a temporary directory to work in.
rm -r results_ASTRAL_PRO/temp_treeshrink_ASTRAL_Pro/
mkdir -p results_ASTRAL_PRO/temp_treeshrink_ASTRAL_Pro/
cd results_ASTRAL_PRO/temp_treeshrink_ASTRAL_Pro/
#Create an input directory structure as required by TreeShrink.
#Loop through the genelist_strict and create a directory for each gene with related gene tree and gene alignment.
while read gene
do
  mkdir -p ${gene}
  cp ../results_iqtree_gene_trees_preliminary/result_${gene}_iqtree.treefile ${gene}/
  mv ${gene}/result_${gene}_iqtree.treefile ${gene}/input.tree
  cp ../results_trimal_preliminary/${gene}.paralogs.trimal.fasta ${gene}/
  mv ${gene}/${gene}.paralogs.trimal.fasta ${gene}/input.fasta
done < ../../6d.ASTRAL_PRO_genelist_succesfull.txt

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
done < ../../6d.ASTRAL_PRO_genelist_succesfull.txt

#Save the output from TreeShrink for later use.
cp */*_treeshrink_*put.* ../results_treeshrink/
cp treeshrinklog.txt ../results_treeshrink/
mv output_summary.txt treeshrink_output_summary.txt
cp treeshrink_output_summary.txt ../results_treeshrink/

#Move up one directory and remove the temporary directory for TreeShrink.
cd ..
rm -r temp_treeshrink_ASTRAL_Pro/
cd ..

#Again use MAFFT to align all updated fasta files.
parallel -j 45 'mafft --quiet results_ASTRAL_PRO/results_treeshrink/{}_treeshrink_output.fasta > results_ASTRAL_PRO/results_mafft_final/{}.paralogs.mafft.fasta' :::: 6d.ASTRAL_PRO_genelist_succesfull.txt

#Again use trimal to trim mafft sequence alignments.
parallel -j 45 'trimal -in results_ASTRAL_PRO/results_mafft_final/{}.paralogs.mafft.fasta -out results_ASTRAL_PRO/results_trimal_final/{}.paralogs.trimal.fasta -resoverlap 0.75 -seqoverlap 0.90 -gt 0.90' :::: 6d.ASTRAL_PRO_genelist_succesfull.txt
#Note that trimal updates sequence names in fasta files just as we need them, with anything after the space removed.

#Update names in the fasta files to shorter names without spaces.
#This is not needed if trimal was just used.
#parallel -j 45 "cat results_ASTRAL_PRO/results_trimal_final/{}.paralogs.trimal.fasta | awk -F '. ' '{print $1}' > results_ASTRAL_PRO/results_trimal_final_renamed/{}.paralogs.trimal.renamed.fasta" :::: 6d.ASTRAL_PRO_genelist_succesfull.txt

#Run IQ-TREE to create gene trees. We use the same settings as for the ASTRAL-III gene tree for species tree analysis before.
parallel -j 45 '~/iqtree-2.1.3-Linux/bin/iqtree2 -s results_ASTRAL_PRO/results_trimal_final/{}.paralogs.trimal.fasta -B 1000 -nm 1000 -m GTR+F+R -nt AUTO --prefix results_ASTRAL_PRO/results_iqtree_gene_trees_final/result_{}_iqtree' :::: 6d.ASTRAL_PRO_genelist_succesfull.txt

#Combine all gene tree files into a single file for ASTRAL-PRO to work with.
cat results_ASTRAL_PRO/results_iqtree_gene_trees_final/*.treefile > 6f.ASTRAL_PRO_input_gene_trees.tre

#Collect all paralog names from all samples and gene trees, translate to sample name, and create an association file for use in ASTRAL-PRO.
Rscript 6g.ASTRAL_PRO_paralogs_to_samples.R

#We now have all the right input to run ASTRAL-PRO.
#Make sure to refer to the correct directories for the Djava library and the actual jar file with exact and full file paths!
java -D"java.library.path=/home/kasper.hendriks/A-pro/ASTRAL-MP/lib" -jar ~/A-pro/ASTRAL-MP/astral.1.1.6.jar -i 6f.ASTRAL_PRO_input_gene_trees.tre -a 6h.ASTRAL_PRO_paralogs_to_samples.txt -o 6i.ASTRAL_PRO_species_tree.tre -t 2 2>6i.ASTRAL_PRO_species_tree.log
