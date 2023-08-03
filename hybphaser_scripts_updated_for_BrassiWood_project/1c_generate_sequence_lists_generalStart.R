####################################
### Generation of sequence lists ###
####################################

####Added by Kasper Hendriks, January 2022.
#When running on a by sample mode on a HPC using GNU Parallel, load locus name from genelist_allGenes_*.txt file as follows.
args = commandArgs(trailingOnly=TRUE)
config_file<-args[1]

# load config
#if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)

# load packages
library(ape)
library(seqinr)
library(stringr)


if(name_for_dataset_optimization_subset != ""){
  folder_subset_add <- paste("_",name_for_dataset_optimization_subset, sep="")
} else {
  folder_subset_add <- ""
} 

output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)
output_sequences <- file.path(path_to_output_folder,paste("03_sequence_lists", folder_subset_add, sep=""))

####KasperH added: load objects failed_loci and failed_samples saved from script 1b.
failed_loci <- readRDS(file=file.path(output_Robjects,"failed_loci.Rds"))
failed_samples <- readRDS(file=file.path(output_Robjects,"failed_samples.Rds"))



if(intronerated_contig=="yes"){
  intronerated_name <- "intronerated" 
  intronerated_underscore <- "_"
} else {
  intronerated_name <- ""
  intronerated_underscore <- ""
}

targets <- read.fasta(fasta_file_with_targets, as.string=TRUE, set.attributes = FALSE)
targets_name <- unique(gsub(".*-","",labels(targets)))
samples <- readLines(path_to_namelist)

###KasperH: Export targets_name for use in following 1c scripts.
write.csv(unlist(targets_name), paste0("genelist_allGenes_",name_for_dataset_optimization_subset,".txt"))

#Load data from previous step 1 scripts a and b.
outsamples_missing <- readRDS(file=file.path(output_Robjects,"outsamples_missing.Rds"))
outloci_missing <- readRDS(file=file.path(output_Robjects,"outloci_missing.Rds"))
####KasperH: added this object for a by-tribe paralog analysis.
if (remove_paralogs == "by_full_dataset") {
  outloci_para_all <- readRDS(file=file.path(output_Robjects,"outloci_para_all.Rds"))
} 
if (remove_paralogs == "by_tribe") {
  outloci_para_by_tribe <- readRDS(file=file.path(output_Robjects,"outloci_para_by_tribe.Rds"))
  outloci_para_each_from_tribe <- readRDS(file=file.path(output_Robjects,"outloci_para_each_from_tribe.Rds"))
}
outloci_para_each <- readRDS(file=file.path(output_Robjects,"outloci_para_each.Rds"))
tab_snps_cl2b <- readRDS(file=file.path(output_Robjects,"Table_SNPs_cleaned.Rds"))

####KasperH added: load objects failed_loci and failed_samples saved from script 1b.
failed_loci <- readRDS(file=file.path(output_Robjects,"failed_loci.Rds"))
failed_samples <- readRDS(file=file.path(output_Robjects,"failed_samples.Rds"))



#############################################

folder4seq_consensus_loci <- file.path(output_sequences,"loci_consensus")
folder4seq_consensus_samples <- file.path(output_sequences,"samples_consensus")
folder4seq_contig_loci <- file.path(output_sequences,"loci_contigs")
folder4seq_contig_samples <- file.path(output_sequences,"samples_contigs")


unlink(c(folder4seq_consensus_loci, folder4seq_consensus_samples,  folder4seq_contig_loci, folder4seq_contig_samples),recursive = T) # delete directory, if it existed in order to prevent errors

dir.create(output_sequences, showWarnings = F)
dir.create(folder4seq_consensus_loci, showWarnings = F)
dir.create(folder4seq_consensus_samples, showWarnings = F)
dir.create(folder4seq_contig_loci, showWarnings = F)
dir.create(folder4seq_contig_samples, showWarnings = F)


