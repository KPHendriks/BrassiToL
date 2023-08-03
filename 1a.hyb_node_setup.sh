#!/bin/bash

# When running on Stampede2 cluster, detail the following:
#   *** Serial Job on Normal Queue ***
#SBATCH	-A TG-TRA170019	                # Account to charge for this run
#SBATCH -J map353too         								# Job name
#SBATCH -o map353too.o%j     								# Name of stdout output file
#SBATCH -e map353too.e%j   									# Name of stderr error file
#SBATCH -p normal          							# Queue (partition) name
#SBATCH --nodes=1
#SBATCH --ntasks=4                      # Value as suggested by Alexandre Zuntini of PAFTOL/Kew
#SBATCH --cpus-per-task=17              # Value as suggested by Alexandre Zuntini of PAFTOL/Kew
#SBATCH -t 48:00:00        							# Run time (hh:mm:ss)
#SBATCH --mail-user=kasper.hendriks@naturalis.nl
#SBATCH --mail-type=none    						# Send no emails

#When working on the Stampede2 cluster, first load relevant modules that contain the required tools.
#module load ooops #See https://portal.tacc.utexas.edu/tutorials/managingio for an explanation
#module load python_cacher #See https://portal.tacc.utexas.edu/tutorials/managingio for an explanation
module load intel/17.0.4
module load gnuparallel
module load biocontainers/0.1.0 	#With this module, various other tools/modules can be loaded
module load bwa/0.7.16a
module load blast/2.6.0
module load samtools/1.10
module load sratoolkit/2.8.2 #Check version! Before I used version 2.10.8
module load trimmomatic/ctr-0.38--1
module load fastq-tools/ctr-0.8--1
#Load Python2 and install Biopython to my ".local" directory in this directory (See email from Erik Ferlanti from TACC Consulting)
module load python2
#pip install --user biopython==1.76 #Run this line at least once before the rest of the script to install Biopython.
module list
pwd
date

#When working on the Stampede2 $SCRATCH partition, first copy some necessary files and tools from the $WORK2 partition to the /tmp directory.
#SPAdes.
cp -r /work2/08079/tg873782/stampede2/SPAdes-3.14.1-Linux.tar.gz /tmp/SPAdes-3.14.1-Linux.tar.gz
tar -xzvf /tmp/SPAdes-3.14.1-Linux.tar.gz -C /tmp/
PATH=/tmp/SPAdes-3.14.1-Linux/bin/:$PATH #HybPiper appears unable to find SPAdes from the module, so I simply put the scripts in a private folder.
#Exonerate.
cp -r /work2/08079/tg873782/stampede2/exonerate-2.2.0-x86_64.tar.gz /tmp/exonerate-2.2.0-x86_64.tar.gz
tar -xzvf /tmp/exonerate-2.2.0-x86_64.tar.gz -C /tmp/
PATH=/tmp/exonerate-2.2.0-x86_64/bin/:$PATH #HybPiper appears unable to find Exonerate from the module, so I simply put the scripts in a private folder.
#HybPiper.
cp -r /work2/08079/tg873782/stampede2/HybPiper.tar.gz /tmp/HybPiper.tar.gz
tar -xzvf /tmp/HybPiper.tar.gz -C /tmp/
#And the folder containing all reference genomes (also copy to SCRATCH for use in Trimmomatic)
cp -r /work2/08079/tg873782/stampede2/Reference_genomes/ /tmp/
cp -r /work2/08079/tg873782/stampede2/Reference_genomes/ /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/

#Save current path as a variable
batch_path="$PWD"

#Start mapping data using HybPiper. Use GNU parallels to do so in parallel, four samples at a time.
#Mapping four samples at a time, 17 CPUs per sample, was suggested by Alexandre Zuntini of PAFTOL/Kew.
parallel -j 4 bash 1b.hyb_sample.sh {} :::: namelist.txt

#And make a note that this batch of samples to map was completed.
#echo $PWD > Finished_batch.txt
#cat /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Finished_batches.txt Finished_batch.txt >> /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Finished_batches.txt
