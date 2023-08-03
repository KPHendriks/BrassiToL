#!/bin/bash

i=$1

#Create a directory for this sample to work in and move there.
mkdir -p ${i}; cd ${i}

####################################################################################
#0. Only use data trimmed already on Stampede2 cluster before.
#For downloading from SRA or running Trimmomatic first, see the script prepared for the Stampede2 cluster.

echo "Start downloading/copying and sort data for sample $i"

if [[ $i == *sl ]]
then
  echo "Combine multiple libraries for updated mapping in HybPiper."
  #First copy all libraries from the list.
  grep -F $i ../namelist_sl.txt > libraries_to_combine.txt
  #Now copy trimmed data files to current directory.
  lib1=$(cat libraries_to_combine.txt | sed 's/|/ /' | awk '{print $2}')
  cp ~/all_data_trimmed/${lib1}_R*paired.fastq.gz .
  lib2=$(cat libraries_to_combine.txt | sed 's/|/ /' | awk '{print $3}')
  cp ~/all_data_trimmed/${lib2}_R*paired.fastq.gz .
  lib3=$(cat libraries_to_combine.txt | sed 's/|/ /' | awk '{print $4}')
  cp ~/all_data_trimmed/${lib3}_R*paired.fastq.gz .
  lib4=$(cat libraries_to_combine.txt | sed 's/|/ /' | awk '{print $5}')
  cp ~/all_data_trimmed/${lib4}_R*paired.fastq.gz .
  lib5=$(cat libraries_to_combine.txt | sed 's/|/ /' | awk '{print $6}')
  cp ~/all_data_trimmed/${lib5}_R*paired.fastq.gz .
  #Now combine all forward and all reverse files.
  cat *_R1_paired.fastq.gz > ${i}_R1_paired.fastq.gz
  cat *_R2_paired.fastq.gz > ${i}_R2_paired.fastq.gz
  cat *_R1_unpaired.fastq.gz > ${i}_R1_unpaired.fastq.gz
  cat *_R2_unpaired.fastq.gz > ${i}_R2_unpaired.fastq.gz
else
  if  test -f ~/all_data_trimmed/"$i"_R1_paired.fastq.gz;
  #First check if data were downloaded/copied and trimmed before. If so, take trimmed data.
  then
    cp ~/all_data_trimmed/"$i"*.fastq.gz .
    #In case we use only non-"-2"-type libraries, remove again the files which contain "-" in their names.
    if [[ $i != *"-"* ]]
    then
      rm "$i"*"-2"*.fastq.gz
    fi
  fi
fi

#Combine unpaired read files for use in HybPiper later (as suggested by Alexandre Zuntini by email on 28 July 2021).
cat ${i}_R*_unpaired.fastq.gz > ${i}_unpaired.fastq.gz

#Unzip the (un)paired read files for proper use in HybPiper.
gunzip -c "$i"_R1_paired.fastq.gz > "$i"_R1_paired.fastq
gunzip -c "$i"_R2_paired.fastq.gz > "$i"_R2_paired.fastq
gunzip -c "$i"_unpaired.fastq.gz > "$i"_unpaired.fastq

#Cleanup
rm *.fastq.gz

#Count number of (un)paired reads from Trimmomatic and log.
echo "$i" > namelist_i.txt
echo $"Trimmed_before" > Raw_reads_count_sample_i.txt
echo $(cat "$i"_R1_paired.fastq|wc -l)/4|bc > Raw_reads_paired_count_sample_i.txt
echo $(cat "$i"_unpaired.fastq|wc -l)/4|bc > Raw_reads_unpaired_count_sample_i.txt
paste namelist_i.txt Raw_reads_count_sample_i.txt Raw_reads_paired_count_sample_i.txt Raw_reads_unpaired_count_sample_i.txt > Read_counts_i.txt
cat Read_counts_i.txt >> /home/kasper.hendriks/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Read_counts.txt

####################################################################################
#1. Map data against Nikolov et al. (2019) reference genome (764 genes from 1827 exons, reference file as created by NikH) translated genes and save results.

#Create a Nikolov-specific directory to work in.
mkdir -p Nikolov2019/
cd Nikolov2019/
mv ../*.fastq .

echo "Start mapping data from bait set Nikolov et al. (2019) for sample $i"

#Log starting time for mapping.
start=$(date +"%T")

#Run HybPiper to map new read data against the Nikolov et al. (2019) reference genome.
~/HybPiper/reads_first.py -r "$i"_R*_paired.fastq --unpaired "$i"_unpaired.fastq -b ../../../Reference_genomes/NikHay_Nikolov_genes_reference.fasta --prefix $i --bwa

#Log finish time for mapping and write to file.
end=$(date +"%T")
echo $i " ,Time taken to map Nikolov genes using bwa," $start "," $end >> /home/kasper.hendriks/2021-08-09_Map_for_genus_level_phylogeny_on_Stampede2/Logged_times.txt

#Study and log possible paralogs for this sample.
python ~/HybPiper/paralog_investigator.py $i

#Collect results into single directory and create tar.gz file to store these.
mkdir -p ${i}_paralog_results/
cp ${i}/genes_with_paralog_warnings.txt ${i}_paralog_results
cat ${i}/*/${i}/paralog_warning.txt > ${i}_paralog_results/${i}_paralog_warnings.txt
cp ${i}/*/${i}/paralogs/*_paralogs.fasta ${i}_paralog_results
#Create tar.gz and copy to SCRATCH to store these intronerate results for possible later use.
tar -czvf ${i}_paralog_results.tar.gz ${i}_paralog_results
cp ${i}_paralog_results.tar.gz ../../../results_paralogs_Nikolov2019/

#Study and log introns for sample i.
python ~/HybPiper/intronerate.py --prefix $i

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
cp ${i}_intronerate_results.tar.gz ../../../results_intronerate_Nikolov2019/

#Cleanup mapping results from HybPiper reads_first.py output.
python ~/HybPiper/cleanup.py ${i}

#Copy results back to SCRATCH partition.
tar -czvf ${i}_mapped_Nikolov2019.tar.gz ${i}
#cp ${i}_mapped_Nikolov2019.tar.gz ../../../mapped_Nikolov2019/

mv *.fastq ../
#Move up a directory before next mapping round.
cd ..

echo "Finished mapping data for sample $i against reference EXONS genome from Nikolov et al. (2019)"


####################################################################################
#2. Map data against Angiosperms353 New Targets reference genome and save results.

#Create a Angiosperms353-specific directory to work in.
mkdir -p Angiosperms353/
cd Angiosperms353/
mv ../*.fastq .

echo "Start mapping data from bait set Angiosperms353 for sample $i"

start=$(date +"%T")

# Run HybPiper to map new read data against the Angiosperms353 reference genome.
~/HybPiper/reads_first.py -r "$i"_R*_paired.fastq --unpaired "$i"_unpaired.fastq -b ../../../Reference_genomes/A353_NewTargets_refgenome.fasta --prefix $i --bwa

end=$(date +"%T")
echo $i ",Time taken to map Angiosperms353 New Targets using bwa," $start "," $end >> ../../../Logged_times.txt

#Study and log possible paralogs for this sample.
python ~/HybPiper/paralog_investigator.py $i

#Collect results into single directory and create tar.gz file to store these.
mkdir -p ${i}_paralog_results
cp ${i}/genes_with_paralog_warnings.txt ${i}_paralog_results
cat ${i}/*/${i}/paralog_warning.txt > ${i}_paralog_results/${i}_paralog_warnings.txt
cp ${i}/*/${i}/paralogs/*_paralogs.fasta ${i}_paralog_results
#Create tar.gz and copy to SCRATCH to store these intronerate results for possible later use.
tar -czvf ${i}_paralog_results.tar.gz ${i}_paralog_results
cp ${i}_paralog_results.tar.gz ../../../results_paralogs_Angiosperms353NewTargets/

#Study and log introns for sample i.
python ~/HybPiper/intronerate.py --prefix $i

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
cp ${i}_intronerate_results.tar.gz ../../../results_intronerate_Angiosperms353NewTargets/

#Cleanup mapping results from HybPiper reads_first.py output.
python ~/HybPiper/cleanup.py ${i}

#Copy results back to SCRATCH partition.
tar -czvf ${i}_mapped_Angiosperms353Newtargets.tar.gz ${i}
#cp ${i}_mapped_Angiosperms353Newtargets.tar.gz ../../../mapped_Angiosperms353NewTargets/

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

#Run HybPiper to map new read data against the Nikolov et al. (2019) reference genome.
~/HybPiper/reads_first.py -r "$i"_R*_paired.fastq --unpaired "$i"_unpaired.fastq -b ../../../Reference_genomes/NikHay_chloro_bait_brassicaceae_new.fasta --prefix $i --bwa

end=$(date +"%T")
echo $i ",Time taken to map chloroplast genes using bwa," $start "," $end >> ../../../Logged_times.txt

#Study and log possible paralogs for this sample.
python ~/HybPiper/paralog_investigator.py $i

#Collect results into single directory and create tar.gz file to store these.
mkdir -p ${i}_paralog_results
cp ${i}/genes_with_paralog_warnings.txt ${i}_paralog_results
cat ${i}/*/${i}/paralog_warning.txt > ${i}_paralog_results/${i}_paralog_warnings.txt
cp ${i}/*/${i}/paralogs/*_paralogs.fasta ${i}_paralog_results
#Create tar.gz and copy to SCRATCH to store these intronerate results for possible later use.
tar -czvf ${i}_paralog_results.tar.gz ${i}_paralog_results
cp ${i}_paralog_results.tar.gz ../../../results_paralogs_chloroplast_HybPiper/

#Study and log introns for sample i.
python ~/HybPiper/intronerate.py --prefix $i

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
cp ${i}_intronerate_results.tar.gz ../../../results_intronerate_chloroplast_HybPiper/

#Cleanup mapping results from HybPiper reads_first.py output.
python ~/HybPiper/cleanup.py ${i}

#Copy results back to SCRATCH partition.
tar -czvf ${i}_mapped_chloroplast_HybPiper.tar.gz ${i}
cp ${i}_mapped_chloroplast_HybPiper.tar.gz ../../../mapped_chloroplast_HybPiper/

#Move up a directory before next mapping round.
cd ..

echo "Finished mapping data against chloroplast reference for sample $i"


####################################################################################
'
#Move up one directory and remove sample i's directory from /tmp/
cd ..
rm -r ${i}/
:'
#To clear space on Naturalis HighMem cluster, now also remove trimmed data used as input.
rm ~/all_data_trimmed/${lib1}_R*paired.fastq.gz
rm ~/all_data_trimmed/${lib2}_R*paired.fastq.gz
rm ~/all_data_trimmed/${lib3}_R*paired.fastq.gz
rm ~/all_data_trimmed/${lib4}_R*paired.fastq.gz
rm ~/all_data_trimmed/${lib5}_R*paired.fastq.gz
'
echo "Finished all mapping analyses for sample $i"
