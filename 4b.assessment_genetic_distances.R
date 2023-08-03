#With this R script we calculate genetic distances, substitution saturation, clocklikeness, and several other metrics to assess 
#suitability of genes for inclusion in the phylogenomic (ASTRAL/IQ-TREE) steps to follow.


#Load packages. ----
library(stringr)
library(seqinr)
library(ape)
library(Rmisc)

#Get data. ----
#In the calculations below we will use the 1,018 gene alignments that have been trimmed and masked following the
#script "3a.ASTRAL_III_inclusive.sh". However, to assess saturation properly, samples with short reads needed to be removed from the alignments.

#Get gene alignment file names and pull gene names from those.
gene_alignment_files<-list.files(path = "results_ASTRALIII_inclusive/results_taper_final_inclusive", full.names = F)
gene_alignment_list<-sapply(str_split(gene_alignment_files, "_"), "[", 1)

#Create a dataframe to store results.
genetic_distances<-data.frame(genes=gene_alignment_list,
                              no_samples_before=NA,
                              no_samples_cleaned=NA,
                              no_samples_removed=NA,
                              fraction_samples_removed=NA,
                              gene_length=NA,
                              dist_raw_mean=NA,
                              dist_raw_min=NA,
                              dist_raw_max=NA,
                              dist_raw_stdev=NA,
                              dist_raw_CI_95_lower=NA,
                              dist_raw_CI_95_upper=NA)

#Calculate genetic distances from alignments. ----
for (i in 1:length(gene_alignment_list)){
  # for (i in 1:10){
    tryCatch({
    #Load gene alignment.
    gene<-gene_alignment_list[i]
    fasta<-read.alignment(paste0("results_ASTRALIII_inclusive/results_taper_final_inclusive/", gene_alignment_files[i]), format="fasta")
    #Remove sequences with >20% missing data from fasta file and save as new fast file.
    keep_seqs<-which(!unlist(lapply(fasta$seq, function(x) str_count(x[1], "-|n|N")))/nchar(fasta$seq[[1]])>0.2)
    fasta_cleaned<-as.alignment(fasta$seq[keep_seqs])
    #Save the updated fasta file to work with later; note the subsetting of sequence names needs to follow the subsetting of genes,
    #or otherwise read data will be connected to the wrong samples!
    write.fasta(as.list(fasta_cleaned$seq), file.out = paste0("results_assessment_saturation_and_clocklikeness/input_gene_alignments_cleaned/" ,gene, "_cleaned.fasta"), names = fasta$nam[keep_seqs])
    #Calculate pairwise genetic distances.
    dist<-dist.dna(as.DNAbin(fasta_cleaned), pairwise.deletion = T, model="raw")
    #List number of samples before and after cleaning.
    genetic_distances$no_samples_before[i]<-length(fasta$seq)
    genetic_distances$no_samples_cleaned[i]<-length(fasta_cleaned$seq)
    genetic_distances$no_samples_removed[i]<-genetic_distances$no_samples_before[i]-genetic_distances$no_samples_cleaned[i]
    genetic_distances$fraction_samples_removed[i]<-1-(genetic_distances$no_samples_cleaned[i]/genetic_distances$no_samples_before[i])
    genetic_distances$gene_length[i]<-nchar(fasta_cleaned$seq[1])
    #Calculate summary metrics and store in overall dataframe.
    genetic_distances$dist_raw_mean[i]<-mean(na.omit(dist[!is.infinite(dist)]))
    genetic_distances$dist_raw_min[i]<-min(na.omit(dist[!is.infinite(dist)]))
    genetic_distances$dist_raw_max[i]<-max(na.omit(dist[!is.infinite(dist)]))
    genetic_distances$dist_raw_stdev[i]<-sd(na.omit(dist[!is.infinite(dist)]))
    genetic_distances$dist_raw_CI_95_lower[i]<-CI(na.omit(dist[!is.infinite(dist)]), ci=0.95)[3]
    genetic_distances$dist_raw_CI_95_upper[i]<-CI(na.omit(dist[!is.infinite(dist)]), ci=0.95)[1]
    #Print a progress update to the screen.
    print(paste0("Finished calculations for gene ",i ," of in total ", length(gene_alignment_list)," genes (",round(100*i/length(gene_alignment_list),1),"%)."))
  })
}

#Save the table as an R object. ----
saveRDS(genetic_distances, file = "results_assessment_saturation_and_clocklikeness/results_R_objects/genetic_distances.Rds")


