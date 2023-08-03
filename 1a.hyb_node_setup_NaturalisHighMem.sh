#!/bin/bash

#SBATCH -J mapSRR         								# Job name
#SBATCH -o mapSRR.o%j     								# Name of stdout output file
#SBATCH -e mapSRR.e%j   									# Name of stderr error file
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12

#When working on the Naturalis HighMem cluster, first run these lines.
#First make sure that the tools and their correct (latest) versions can be found.
#Some were not on the Naturalis HighMem cluster yet, so these I installed to my own home folder as follows.
PATH=~/SPAdes-3.14.1-Linux/bin/:$PATH
PATH=~/samtools-1.10/bin/:$PATH
PATH=~/sratoolkit.2.10.8-ubuntu64/bin/:$PATH
PATH=~/exonerate-2.2.0-x86_64/bin/:$PATH

#It is possible to define a namelist.txt or a namelist_sl.txt.
#The first has a single column only with libraries that need to be mapped.
#The second has in the first column the "sl" library name, and in subsequent columns all libraries that need to be combined after trimming and before mapping.
#In the second case, create a simple namelist.txt from the namelist_sl.txt file.
cat namelist_sl.txt | sed 's/|/ /' | awk '{print $1}' > namelist.txt

#Start mapping data using HybPiper. Use GNU parallel to do so in parallel, four samples at a time.
#Mapping x samples at a time, ~17 CPUs per sample, was suggested by Alexandre Zuntini of PAFTOL/Kew.
parallel -j 1 bash 1b.hyb_sample_NaturalisHighMem.sh {} :::: namelist.txt

#And make a note that this batch of samples to map was completed.
echo $PWD > Finished_batch.txt
