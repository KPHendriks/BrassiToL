#!/bin/bash

i=$1

#Create a directory for this sample to work in and move there.
mkdir -p ${i}; cd ${i}

#Save the location's path.
batch_path2="$PWD"

#Create a single-sample text file for the current sample.
echo "$i" > namelist_i.txt

####################################################################################
#0. Download from SRA or copy from cluster storage the raw data to map, and sort data.

echo "Start downloading/copying and sort data for sample $i"

if  test -f /scratch/08079/tg873782/all_data_trimmed_moveToScrath/"$i"_R1_paired.fastq.gz;
#First check if data were downloaded/copied and trimmed before. If so, take trimmed data.
then
  cp /scratch/08079/tg873782/all_data_trimmed_moveToScrath/"$i"*.fastq.gz .
  #In case we use only non-"-2"-type libraries, remove again the files which contain "-" in their names.
  if [[ $i != *"-"* ]]
  then
    rm "$i"*"-2"*.fastq.gz
  fi
  echo $"Trimmed_before" > Raw_reads_count_sample_i.txt
else
  if [[ $i == "SRR"* ]]
  #Check if the raw data is on GenBank's SRA database, and if so, pull data from there directly to the cluster.
  then
    #Use SRA-toolkit to download data and immediately split so that fastq files become available directly.
    module load sratoolkit/2.8.2 #Check version! Before I used version 2.10.8
    fastq-dump "$i" --split-files
    #Note that the faster tool fasterq-dump is available only from version 2.9.1 only so need to ask TACC to add that to the sra module.
    #fasterq-dump --split-files -f -O ./ "$i"

    #Zip these files for proper use in Trimmomatic.
    gzip -c "$i"_1.fastq > "$i"_R1.fastq.gz
    gzip -c "$i"_2.fastq > "$i"_R2.fastq.gz
  elif [[ $i == PAFTOL_* ]]
  then
    #Copy all raw data files that belong to sample i to our current directory.
    cp /work2/08079/tg873782/stampede2/data/all_data/"$i"*.fastq.gz .
    #Rename files for use in Trimmomatic
    mv "$i"*_1.fastq.gz "$i"_R1.fastq.gz
    mv "$i"*_2.fastq.gz "$i"_R2.fastq.gz
  elif [[ $i == ERR* ]]
  then
    #Copy all raw data files that belong to sample i to our current directory.
    cp /work2/08079/tg873782/stampede2/data/all_data/"$i"*.fastq.gz .
    #Rename files for use in Trimmomatic
    mv "$i"*_1.fastq.gz "$i"_R1.fastq.gz
    mv "$i"*_2.fastq.gz "$i"_R2.fastq.gz
  else
    #Copy all raw data files that belong to sample i to our current directory.
    cp /work2/08079/tg873782/stampede2/data/all_data/"$i"*.fastq.gz .
    #In case we use only non-"-2"-type libraries, remove again the files which contain "-" in their names.
    if [[ $i != *"-"* ]]
    then
      rm "$i"*"-2"*.fastq.gz
    fi
    #Cat all forward and reverse files and create a new file.
    cat "$i"*_R1.fastq.gz > "$i"_1.fastq.gz
    cat "$i"*_R2.fastq.gz > "$i"_2.fastq.gz
    #Remove the original input files.
    rm "$i"*_R1.fastq.gz "$i"*_R2.fastq.gz
    #And rename these for proper use in Trimmomatic.
    mv "$i"_1.fastq.gz "$i"_R1.fastq.gz
    mv "$i"_2.fastq.gz "$i"_R2.fastq.gz

    #BELOW CAN BE DELETED?
    #Sort fastq files. This can be needed when raw files from different sequencing runs were combined.
    #This has happened e.g. when sequencing company BaseClear ran a sequencing pool twice to add data, but then added data in different orders.
    #Unpack the files before sorting.
    #gunzip "$i"_1.fastq.gz "$i"_2.fastq.gz
    #Load module fastq-tools again. Strangely, it seems to be unloaded when starting a parallel run...
    #module load fastq-tools/ctr-0.8--1
    #fastq-sort "$i"_1.fastq > "$i"_R1_sorted.fastq
    #fastq-sort "$i"_2.fastq > "$i"_R2_sorted.fastq
    #Zip these new files again for proper use in Trimmomatic.
    #gzip -c "$i"_R1_sorted.fastq > "$i"_R1.fastq.gz
    #gzip -c "$i"_R2_sorted.fastq > "$i"_R2.fastq.gz

  fi
  echo "Finished downloading/copying and sort data for sample $i"

  #Count number of read pairs and log.
  echo $(zcat "$i"_R1.fastq.gz|wc -l)/4|bc > Raw_reads_count_sample_i.txt

  echo "Start trimming data for sample $i using Trimmomatic"

  # Trim reads using Trimmomatic.
  start=$(date +"%T") #Log time it takes to run Trimmomatic.
  # Settings for trimming as used by Hendriks, ..., Bailey (2021) APPS.
  module load trimmomatic/ctr-0.38--1
  trimmomatic PE -phred33 "$i"_R1.fastq.gz "$i"_R2.fastq.gz "$i"_R1_paired.fastq.gz "$i"_R1_unpaired.fastq.gz "$i"_R2_paired.fastq.gz "$i"_R2_unpaired.fastq.gz ILLUMINACLIP:../Reference_genomes/illumina_adapters_for_trimmomatic_normal_and_palindrome_mode.fasta:2:30:10:2:true LEADING:10  TRAILING:10  SLIDINGWINDOW:4:20  MINLEN:40 2> log_trim_${i}.txt
  #OLD trimmomatic PE -phred33 "$i"_R1.fastq.gz "$i"_R2.fastq.gz "$i"_R1_paired.fastq.gz "$i"_R1_unpaired.fastq.gz "$i"_R2_paired.fastq.gz "$i"_R2_unpaired.fastq.gz ILLUMINACLIP:../Reference_genomes/TruSeq3-PE.fa:2:30:10:2:true LEADING:10  TRAILING:10  SLIDINGWINDOW:4:20  MINLEN:40
  end=$(date +"%T") #Again log time it takes to run Trimmomatic.
  echo $i " ,Time taken to run Trimmomatic," $start "," $end >> /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Logged_times.txt
  #Save log file for sample i.
  mv log_trim_*.txt /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Log_trimmomatic/
  #Save the files of trimmed and (un)paired reads for later use in further mapping activities, so that Trimmomatic step can then be skipped.
  cp *paired.fastq.gz /scratch/08079/tg873782/all_data_trimmed_moveToScrath/

fi

#Combine unpaired read files for use in HybPiper later (as suggested by Alexandre Zuntini by email on 28 July 2021).
cat "$i"_R*_unpaired.fastq.gz > "$i"_unpaired.fastq.gz

#Unzip the (un)paired read files for proper use in HybPiper.
gunzip -c "$i"_R1_paired.fastq.gz > "$i"_R1_paired.fastq
gunzip -c "$i"_R2_paired.fastq.gz > "$i"_R2_paired.fastq
gunzip -c "$i"_unpaired.fastq.gz > "$i"_unpaired.fastq

#Count number of (un)paired reads from Trimmomatic and log.
echo $(cat "$i"_R1_paired.fastq|wc -l)/4|bc > Raw_reads_paired_count_sample_i.txt
echo $(cat "$i"_unpaired.fastq|wc -l)/4|bc > Raw_reads_unpaired_count_sample_i.txt
paste namelist_i.txt Raw_reads_count_sample_i.txt Raw_reads_paired_count_sample_i.txt Raw_reads_unpaired_count_sample_i.txt > Read_counts_i.txt
cat Read_counts_i.txt >> /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Read_counts.txt

#Cleanup
rm *.fastq.gz
rm Raw_reads_paired_count_sample_i.txt Raw_reads_count_sample_i.txt Raw_reads_unpaired_count_sample_i.txt Read_counts_i.txt

echo "Finished trimming data for sample $i using Trimmomatic"

#Create another sample-specific directory within the node's /tmp/ directory.
mkdir -p /tmp/${i}/
#Copy files and move to directory /tmp/ to run HybPiper to prevent excessive IO-issues.
cp "$i"_R*_paired.fastq /tmp/${i}/
cp "$i"_unpaired.fastq /tmp/${i}/
cd /tmp/${i}/


####################################################################################
#1. Map data against Nikolov et al. (2019) reference genome (764 genes from 1827 exons, reference file as created by NikH) translated genes and save results.

#Create a Nikolov-specific directory to work in.
mkdir -p Nikolov2019/
cd Nikolov2019/
cp ../*.fastq .

echo "Start mapping data from bait set Nikolov et al. (2019) for sample $i"

#Log starting time for mapping.
start=$(date +"%T")

module load bwa/0.7.16a
module load samtools/1.10

#Run HybPiper to map new read data against the Nikolov et al. (2019) reference genome.
python /tmp/HybPiper/reads_first.py -r "$i"_R*_paired.fastq --unpaired "$i"_unpaired.fastq -b /tmp/Reference_genomes/NikHay_Nikolov_genes_reference.fasta --prefix $i --bwa

#Log finish time for mapping and write to file.
end=$(date +"%T")
echo $i " ,Time taken to map Nikolov genes using bwa," $start "," $end >> /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Logged_times.txt

#Study and log possible paralogs for this sample.
python /tmp/HybPiper/paralog_investigator.py $i

#Collect results into single directory and create tar.gz file to store these.
mkdir -p ${i}_paralog_results/
cp ${i}/genes_with_paralog_warnings.txt ${i}_paralog_results
cat ${i}/*/${i}/paralog_warning.txt > ${i}_paralog_results/${i}_paralog_warnings.txt
cp ${i}/*/${i}/paralogs/*_paralogs.fasta ${i}_paralog_results
#Create tar.gz and copy to SCRATCH to store these intronerate results for possible later use.
tar -czvf ${i}_paralog_results.tar.gz ${i}_paralog_results
cp ${i}_paralog_results.tar.gz /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/results_paralogs_Nikolov2019/

#Study and log introns for sample i.
python /tmp/HybPiper/intronerate.py --prefix $i

#Collect results into single directory and create tar.gz file to store these.
mkdir -p ${i}_intronerate_results
cp ${i}/${i}_genes.gff ${i}_intronerate_results
cp ${i}/intron_stats.txt ${i}_intronerate_results
cp ${i}/exonerate_genelist.txt ${i}_intronerate_results
cp ${i}/genes_with_seqs.txt ${i}_intronerate_results
cp ${i}/genes_with_paralog_warnings.txt ${i}_intronerate_results
cp ${i}/*/${i}/sequences/intron/*.fasta ${i}_intronerate_results
#Create tar.gz and copy to SCRATCH to store these intronerate results for possible later use.
tar -czvf ${i}_intronerate_results.tar.gz ${i}_intronerate_results
cp ${i}_intronerate_results.tar.gz /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/results_intronerate_Nikolov2019/

#Cleanup mapping results from HybPiper reads_first.py output.
python /tmp/HybPiper/cleanup.py ${i}

#Copy results back to SCRATCH partition.
tar -czvf ${i}_mapped_Nikolov2019.tar.gz ${i}
cp ${i}_mapped_Nikolov2019.tar.gz /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/mapped_Nikolov2019/
#And unpack the tar file in the right directory on SCRATCH.
batch_path3="$PWD"
cd /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/mapped_Nikolov2019/
tar -xf ${i}_mapped_Nikolov2019.tar.gz
cd $batch_path3/

#Move up a directory before next mapping round.
cd ..

echo "Finished mapping data for sample $i against reference EXONS genome from Nikolov et al. (2019)"


####################################################################################
#2. Map data against Angiosperms353 New Targets reference genome and save results.

#Create a Angiosperms353-specific directory to work in.
mkdir -p Angiosperms353/
cd Angiosperms353/
cp ../*.fastq .

echo "Start mapping data from bait set Angiosperms353 for sample $i"

start=$(date +"%T")

module load bwa/0.7.16a
module load samtools/1.10

# Run HybPiper to map new read data against the Angiosperms353 reference genome.
python /tmp/HybPiper/reads_first.py -r "$i"_R*_paired.fastq --unpaired "$i"_unpaired.fastq -b /tmp/Reference_genomes/A353_NewTargets_refgenome.fasta --prefix $i --bwa

end=$(date +"%T")
echo $i ",Time taken to map Angiosperms353 New Targets using bwa," $start "," $end >> /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Logged_times.txt

#Study and log possible paralogs for this sample.
python /tmp/HybPiper/paralog_investigator.py $i

#Collect results into single directory and create tar.gz file to store these.
mkdir -p ${i}_paralog_results
cp ${i}/genes_with_paralog_warnings.txt ${i}_paralog_results
cat ${i}/*/${i}/paralog_warning.txt > ${i}_paralog_results/${i}_paralog_warnings.txt
cp ${i}/*/${i}/paralogs/*_paralogs.fasta ${i}_paralog_results
#Create tar.gz and copy to SCRATCH to store these intronerate results for possible later use.
tar -czvf ${i}_paralog_results.tar.gz ${i}_paralog_results
cp ${i}_paralog_results.tar.gz /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/results_paralogs_Angiosperms353NewTargets/

#Study and log introns for sample i.
python /tmp/HybPiper/intronerate.py --prefix $i

#Collect results into single directory and create tar.gz file to store these.
mkdir -p ${i}_intronerate_results
cp ${i}/${i}_genes.gff ${i}_intronerate_results
cp ${i}/intron_stats.txt ${i}_intronerate_results
cp ${i}/exonerate_genelist.txt ${i}_intronerate_results
cp ${i}/genes_with_seqs.txt ${i}_intronerate_results
cp ${i}/genes_with_paralog_warnings.txt ${i}_intronerate_results
cp ${i}/*/${i}/sequences/intron/*.fasta ${i}_intronerate_results
#Create tar.gz and copy to SCRATCH to store these intronerate results for possible later use.
tar -czvf ${i}_intronerate_results.tar.gz ${i}_intronerate_results
cp ${i}_intronerate_results.tar.gz /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/results_intronerate_Angiosperms353NewTargets/

#Cleanup mapping results from HybPiper reads_first.py output.
python /tmp/HybPiper/cleanup.py ${i}

#Copy results back to SCRATCH partition.
tar -czvf ${i}_mapped_Angiosperms353Newtargets.tar.gz ${i}
cp ${i}_mapped_Angiosperms353Newtargets.tar.gz /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/mapped_Angiosperms353NewTargets/
#And unpack the tar file in the right directory on SCRATCH.
batch_path4="$PWD"
cd /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/mapped_Angiosperms353NewTargets/
tar -xf ${i}_mapped_Angiosperms353Newtargets.tar.gz
cd $batch_path4/

#Move up a directory before next mapping round.
cd ..

echo "Finished mapping data from bait set Angiosperms353 for sample $i"

:'
####################################################################################
#3. Map data against chloroplast reference created by Nikolai Hai

#Create a chloroplast-specific directory to work in.
mkdir -p chloroplast/
cd chloroplast/
cp ../*.fastq .

echo "Start mapping data against chloroplast reference for sample $i"

#Log starting time for mapping.
start=$(date +"%T")

module load bwa/0.7.16a
module load samtools/1.10

#Run HybPiper to map new read data against the Nikolov et al. (2019) reference genome.
python /tmp/HybPiper/reads_first.py -r "$i"_R*_paired.fastq --unpaired "$i"_unpaired.fastq -b /tmp/Reference_genomes/NikHay_chloro_bait_brassicaceae_new.fasta --prefix $i --bwa

end=$(date +"%T")
echo $i ",Time taken to map chloroplast genes using bwa," $start "," $end >> /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Logged_times.txt

#Study and log possible paralogs for this sample.
python /tmp/HybPiper/paralog_investigator.py $i

#Collect results into single directory and create tar.gz file to store these.
mkdir -p ${i}_paralog_results
cp ${i}/genes_with_paralog_warnings.txt ${i}_paralog_results
cat ${i}/*/${i}/paralog_warning.txt > ${i}_paralog_results/${i}_paralog_warnings.txt
cp ${i}/*/${i}/paralogs/*_paralogs.fasta ${i}_paralog_results
#Create tar.gz and copy to SCRATCH to store these intronerate results for possible later use.
tar -czvf ${i}_paralog_results.tar.gz ${i}_paralog_results
cp ${i}_paralog_results.tar.gz /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/results_paralogs_chloroplast_HybPiper/

#Study and log introns for sample i.
python /tmp/HybPiper/intronerate.py --prefix $i

#Collect results into single directory and create tar.gz file to store these.
mkdir -p ${i}_intronerate_results
cp ${i}/${i}_genes.gff ${i}_intronerate_results
cp ${i}/intron_stats.txt ${i}_intronerate_results
cp ${i}/exonerate_genelist.txt ${i}_intronerate_results
cp ${i}/genes_with_seqs.txt ${i}_intronerate_results
cp ${i}/genes_with_paralog_warnings.txt ${i}_intronerate_results
cp ${i}/*/${i}/sequences/intron/*.fasta ${i}_intronerate_results
#Create tar.gz and copy to SCRATCH to store these intronerate results for possible later use.
tar -czvf ${i}_intronerate_results.tar.gz ${i}_intronerate_results
cp ${i}_intronerate_results.tar.gz /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/results_intronerate_chloroplast_HybPiper/

#Cleanup mapping results from HybPiper reads_first.py output.
python /tmp/HybPiper/cleanup.py ${i}

#Copy results back to SCRATCH partition.
tar -czvf ${i}_mapped_chloroplast_HybPiper.tar.gz ${i}
cp ${i}_mapped_chloroplast_HybPiper.tar.gz /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/mapped_chloroplast_HybPiper/
#Unpack the tar file in the right directory on SCRATCH.
batch_path5="$PWD"
cd /scratch/08079/tg873782/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/mapped_chloroplast_HybPiper/
tar -xf ${i}_mapped_chloroplast_HybPiper.tar.gz
cd $batch_path5/

#Move up a directory before next mapping round.
cd ..

echo "Finished mapping data against chloroplast reference for sample $i"
'

####################################################################################

#Move up one directory and remove sample i's directory from /tmp/
cd ..
rm -r ${i}/

#Move back to batch directory on SCRATCH.
cd $batch_path2

#Now also remove sample i's directory from here.
cd ..
rm -r ${i}/

#In case data were downloaded from sra, remove files created by sra-toolkit on $HOME.
rm /home1/08079/tg873782/ncbi/public/sra/${i}.sra
rm /home1/08079/tg873782/ncbi/public/sra/${i}.sra.cache

echo "Finished all mapping analyses for sample $i"
