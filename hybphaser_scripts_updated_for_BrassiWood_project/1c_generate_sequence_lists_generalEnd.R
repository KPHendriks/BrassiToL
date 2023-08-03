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


#########################################################################
### concatenate consensus files across all samples to lists per locus ###
#########################################################################
# this extracts sequences from all subfolders in the HybPiper folder and collates them into one file per locus


####KasperH: this loop was moved to another R script that can be run in parallel for all loci.


# changing interleaved HybPiper files to non-interleaved fasta files
for(file in list.files(folder4seq_contig_loci, full.names = T)){
  file <- list.files(folder4seq_contig_loci, full.names = T)[1] ####KasperH: changed [3] to [1]
  lines <- readLines(file)
  conx <- file(file)
  writeLines(str_split(paste(gsub("(>.*)",":\\1:",lines),collapse =""), pattern = ":")[[1]][-1], conx)
  close(conx)
}


# check whether accessions are in the hybpiper folder but not in the samples list.
# if the sample list is smaller some sequences have to be removed from the loci lists

hybpiper_result_dirs <- list.dirs(file.path(path_to_output_folder,"01_data"), full.names = FALSE, recursive = FALSE)
dirs_not_in_sample_list <- hybpiper_result_dirs[which(!(hybpiper_result_dirs %in% samples))]

if(length(dirs_not_in_sample_list) !=0 ){
  # consensus
  for(raw_consensus_file in list.files(folder4seq_consensus_loci, full.names = TRUE)){
    locus_consensus <- readLines(raw_consensus_file)
    lines_with_samplename <- which(gsub(">","",locus_consensus) %in% dirs_not_in_sample_list)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_red <- locus_consensus[-lines_to_remove]
      conn <- file(raw_consensus_file)
      writeLines(locus_file_red, conn)
      close(conn)
    }
  }
  # contig
  for(raw_contig_file in list.files(folder4seq_contig_loci, full.names = TRUE)){
    locus_contig <- readLines(raw_contig_file)
    lines_with_samplename <- which(gsub(">","",locus_contig) %in% dirs_not_in_sample_list)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_hp_red <- locus_contig[-lines_to_remove]
      conn <- file(raw_contig_file)
      writeLines(locus_file_hp_red, conn)
      close(conn)
    }
  }
}


### remove loci from dataset optimization (missing data and paralogs)
######################################################################

## remove loci (failed/missing data/paralogs for all) for all samples

loci_files_consensus <- list.files(path = folder4seq_consensus_loci, full.names = T )
loci_files_contig <- list.files(path = folder4seq_contig_loci, full.names = T )

####KasperH: updated the loci_to_remove vector below.
if (remove_paralogs == "by_full_dataset") {
  loci_to_remove <- c(names(failed_loci), outloci_missing, outloci_para_all)
} else {
  loci_to_remove <- c(names(failed_loci), outloci_missing)
}
  

if(length(loci_to_remove)!=0){
  loci_files_to_remove_consensus <- loci_files_consensus[which(gsub(".*/(.*)_consensus.fasta","\\1",loci_files_consensus) %in% loci_to_remove)]
  loci_files_to_remove_contig <- loci_files_contig[which(gsub(".*/(.*)_contig.fasta","\\1",loci_files_contig) %in% loci_to_remove)]
  file.remove(loci_files_to_remove_consensus)
  file.remove(loci_files_to_remove_contig)
}


# get vector of all samples that should be removed from every locus

samples_to_remove_4all <- vector()

if(length(failed_samples) > 0 ){
  samples_to_remove_4all <- names(failed_samples)
} 

if(length(outsamples_missing) > 0 ){
  samples_to_remove_4all <-  unique(c(samples_to_remove_4all, outsamples_missing))
} 


## remove samples to be removed from all and sequences in each locus file from paralogs for each sample


for(locus in rownames(tab_snps_cl2b)){
  print(locus)
  if(length(grep(paste("\\b",locus,"\\b",sep=""),outloci_para_each)) >0 ){
    samples_to_remove <- c(samples_to_remove_4all, names(outloci_para_each[grep(paste("\\b",locus,"\\b",sep=""),outloci_para_each)]))
    ####KasperH: add another check if this locus is in the outloci_para_each_from_tribe object.
    if(length(grep(paste("\\b",locus,"\\b",sep=""),outloci_para_each_from_tribe)) >0 ){
      samples_to_remove <- c(samples_to_remove_4all, names(outloci_para_each[grep(paste("\\b",locus,"\\b",sep=""),outloci_para_each_from_tribe)]))
    }
  } else {
    samples_to_remove <- samples_to_remove_4all
  }
  
  
  if(length(samples_to_remove) !=0 ){
    
    # consensus
    locus_consensus <- readLines(file.path(folder4seq_consensus_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"consensus.fasta",sep="")))
    lines_with_samplename <- which(gsub(">","",locus_consensus) %in% samples_to_remove)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_red <- locus_consensus[-lines_to_remove]
      conn <- file(file.path(folder4seq_consensus_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"consensus.fasta",sep="")))
      writeLines(locus_file_red, conn)
      close(conn)
    }
    
    # contig
    locus_contig <- readLines(file.path(folder4seq_contig_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"contig.fasta",sep="")))
    lines_with_samplename <- which(gsub(">","",locus_contig) %in% samples_to_remove)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_hp_red <- locus_contig[-lines_to_remove]
      conn <- file(file.path(folder4seq_contig_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"contig.fasta",sep="")))
      writeLines(locus_file_hp_red, conn)
      close(conn)
    }
    
  }
} 

