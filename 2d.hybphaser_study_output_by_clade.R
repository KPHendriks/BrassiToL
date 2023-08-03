### Kasper Hendriks et al. (2022)
### Brassicaceae genus level family phylogeny project.


### INTRODUCTION

#This script is used to study the output from HybPhaser in a more detailed way than normally done by the HybPhaser R-scripts.
#This is needed because preliminary checks with the standard HybPhaser scripts showed that several tribes within the
#Brassicaceae family are characterised by very high proportions of SNPs, resulting in very high levels of allele divergence (AD) and
#locus heterozygosity (LH). These are signs of past hybridization events and the presence of many polyploids in our dataset.

#The aim is to filter filter out genes from our dataset that are prone to high levels of paralogy, because paralogs will obscure the
#true phylogenetic relationships within the tribes and family.


### STEP 1: LOAD RESULTS FROM HYBPHASER STEP 1

#The loaded results were created with HybPhaser scripts 1a_count_snps.R and 1b_assess_dataset.R.

#Clear environment
rm(list=ls()) 

#Load config
config_file <- "results_HybPhaser_allGenes/config_inclusive.txt"
source(config_file)

#Load packages
library(ape)
library(ggplot2)
library(cowplot)

#Set paths to results from script 1b_assess_dataset.R
if(name_for_dataset_optimization_subset != ""){
  folder_subset_add <- paste("_",name_for_dataset_optimization_subset, sep="")
}else {
  folder_subset_add <- ""
} 
output_Robjects <- file.path("results_HybPhaser_allGenes",path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)
output_assess <- file.path("results_HybPhaser_allGenes",path_to_output_folder,paste("02_assessment", folder_subset_add, sep=""))

#Load metadata for use in analyses and plotting
metadata <- read.csv(file = "2e.metadata.csv", header = T, stringsAsFactors=FALSE, fileEncoding="latin1")


### STEP 2: CREATE A HEATMAP TO STUDY SNP PROPORTIONS

#Load Rds objects
tab_snps <- readRDS(file=file.path(output_Robjects,"Table_SNPs.Rds"))

#Turn tab_snps table into long format for use with ggplot2
library(reshape2)
library(forcats) #Load to bring order easily to the plot below

#tab_snps <- tab_snps[1:100,1:100] #test on a subset
tab_snps$Gene <- rownames(tab_snps)
tab_snps_long <- melt(tab_snps, id.vars = "Gene", variable.name = "Sample", value.name = "SNP_proportion")

#Add taxonomic details for use in plotting.
tab_snps_long$Tribe <- metadata$Tribe[match(tab_snps_long$Sample, metadata$Library_ID)]
tab_snps_long$Species <- metadata$Name[match(tab_snps_long$Sample, metadata$Library_ID)]
tab_snps_long$Tribe_Species_Sample <- paste(tab_snps_long$Tribe, tab_snps_long$Species, tab_snps_long$Sample)

#Reorder by tribe for plotting.
tab_snps_long <- tab_snps_long[order(tab_snps_long$Tribe_Species_Sample),]

#Reorder by mean SNP proportion for each gene.
library(dplyr)
snps_mean_by_gene<-tab_snps_long %>%
  group_by(Gene) %>%
  summarise(Mean=mean(SNP_proportion, na.rm = T))
tab_snps_long$genes_mean_snp_proportion<-snps_mean_by_gene$Mean[match(tab_snps_long$Gene, snps_mean_by_gene$Gene)]
tab_snps_long <- tab_snps_long[order(tab_snps_long$genes_mean_snp_proportion),]

# tab_snps_long$Sample <- as.character(tab_snps_long$Sample)
# tab_snps_long$Tribe <- as.character(tab_snps_long$Tribe)

#Create a heatmap.
snps_heatmap <- ggplot(data=tab_snps_long, aes(x=fct_inorder(Gene), y=fct_inorder(Tribe_Species_Sample), fill=SNP_proportion))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red", name = "SNP proportion")+
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+
  #coord_fixed(ratio = 1)+
  theme(axis.text.y = element_text(size=3), axis.text.x = element_text(angle = 45, hjust=1, size=0.2),
        axis.title=element_text(size=14,face="bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  xlab("Gene (ordered by mean SNP proportion over all samples)")+
  ylab("Sample (clustered by tribe)")

#And save the heatmap as an image file.
ggsave(file = "2f.hybphaser_result_snsps_heatmap.pdf", snps_heatmap, width = 12, height = 15)
ggsave(file = "2f.hybphaser_result_snsps_heatmap.png", snps_heatmap, width = 12, height = 15)


### STEP 3: RE-CREATE DATA RECOVERED OVERVIEW FROM HYBPHASER, NOW SPLIT BY BAIT SET AND IN/OUTGROUP

#Load Rds objects
seq_per_sample_prop <- readRDS(file=file.path(output_Robjects,"seq_per_sample_prop.Rds"))
prop_target_length_per_sample <- readRDS(file=file.path(output_Robjects,"prop_target_length_per_sample.Rds"))
seq_per_locus_prop <- readRDS(file=file.path(output_Robjects,"seq_per_locus_prop.Rds"))
prop_target_length_per_locus <- readRDS(file=file.path(output_Robjects,"prop_target_length_per_locus.Rds"))

#Recreate the default figures in "1_Data_recovered_overview.pdf" now with annotation of bait set and in- vs. outgroup

#Put results in a dataframes.
data_recovered_per_sample <- data.frame(Sample = names(seq_per_sample_prop), seq_per_sample_prop = array(seq_per_sample_prop), prop_target_length_per_sample = array(prop_target_length_per_sample))
data_recovered_per_locus <- data.frame(Locus = names(seq_per_locus_prop), seq_per_locus_prop = array(seq_per_locus_prop), prop_target_length_per_locus = array(prop_target_length_per_locus))

#Annotate dataframes.
data_recovered_per_sample$Tribe <- metadata$Tribe[match(data_recovered_per_sample$Sample, metadata$Library_ID)]
data_recovered_per_sample$Family <- data_recovered_per_sample$Tribe
data_recovered_per_sample$Family[grepl("Outgroup", data_recovered_per_sample$Tribe, fixed = T)] <- "non-Brassicaceae"
data_recovered_per_sample$Family[!grepl("Outgroup", data_recovered_per_sample$Tribe, fixed = T)] <- "Brassicaceae"

data_recovered_per_locus$Bait_set <- NA
data_recovered_per_locus$Bait_set[nchar(data_recovered_per_locus$Locus)<5] <- "A353"
data_recovered_per_locus$Bait_set[nchar(data_recovered_per_locus$Locus)>5] <- "B764"

#Calculate mean values by locus and sample for the two bait sets, for mention in the publication's main text.
library(dplyr)

data_recovered_per_locus %>%
  group_by(Bait_set) %>%
  summarise(across(c("seq_per_locus_prop", "prop_target_length_per_locus"), ~ mean(.x, na.rm = TRUE)),
            n = n())

data_recovered_per_sample %>%
  group_by(Family) %>%
  summarise(across(c("seq_per_sample_prop", "prop_target_length_per_sample"), ~ mean(.x, na.rm = TRUE)),
            n = n())

#Now we can recreate the boxplots and point plots.
sample_loci_recovered <- ggplot(data=data_recovered_per_sample, aes(x=Family, y=seq_per_sample_prop, colour=Family))+
  geom_boxplot()+
  xlab("Family")+
  ylab("Proportion of genes per sample")
sample_length_recovered <- ggplot(data=data_recovered_per_sample, aes(x=Family, y=prop_target_length_per_sample, colour=Family))+
  geom_boxplot()+
  xlab("Family")+
  ylab("Proportion of total target length per sample")
sample_loci_vs_length_recovered <- ggplot(data=data_recovered_per_sample, aes(y=seq_per_sample_prop, x=prop_target_length_per_sample, colour=Family))+
  geom_point()+
  xlab("Proportion of total target length per sample")+
  ylab("Proportion of genes per sample")

locus_samples_recovered <- ggplot(data=data_recovered_per_locus, aes(x=Bait_set, y=seq_per_locus_prop, colour=Bait_set))+
  geom_boxplot()+
  scale_colour_manual(values=c("#CC79A7", "#009E73"), name = "Bait set")+
  xlab("Bait set")+
  ylab("Proportion of genes per sample")
locus_length_recovered <- ggplot(data=data_recovered_per_locus, aes(x=Bait_set, y=prop_target_length_per_locus, colour=Bait_set))+
  geom_boxplot() +
  ylim(c(0,1))+
  #geom_text(x=1.4, y=0.1, label="11 genes from B764 had prop. target length >1", colour= "black", size=3)+
  scale_colour_manual(values=c("#CC79A7", "#009E73"), name = "Bait set")+
  xlab("Bait set")+
  ylab("Proportion of total target length per sample")
locus_loci_vs_length_recovered <- ggplot(data=data_recovered_per_locus, aes(y=seq_per_locus_prop, x=prop_target_length_per_locus, colour=Bait_set))+
  geom_point() +
  xlim(c(0,1)) +
  #geom_text(x=0.4, y=0.15, label="11 genes from B764 had prop. target length >1", colour= "black", size=3)+
  scale_colour_manual(values=c("#CC79A7", "#009E73"), name = "Bait set")+
  xlab("Proportion of total target length per sample")+
  ylab("Proportion of genes per sample")

#Combine these six graphs into one figure and save figure.
data_recovered_overview <- plot_grid(locus_samples_recovered + theme(legend.position="none"), 
                                     locus_length_recovered + theme(legend.position="none"), 
                                     locus_loci_vs_length_recovered + theme(legend.position="none"),
                                     get_legend(locus_loci_vs_length_recovered),
                                     sample_loci_recovered + theme(legend.position="none"), 
                                     sample_length_recovered + theme(legend.position="none"), 
                                     sample_loci_vs_length_recovered + theme(legend.position="none"),
                                     get_legend(sample_loci_vs_length_recovered),
                                     ncol = 4, nrow = 2, rel_widths = c(5, 5, 7, 3),
                                     labels = c('(a)', '(b)', '(c)', '' ,'(d)', '(e)', '(f)', ''))

ggsave(data_recovered_overview, filename = "2g.hybphaser_result_data_recovery.pdf", width = 12, height = 8)
ggsave(data_recovered_overview, filename = "2g.hybphaser_result_data_recovery.png", width = 12, height = 8)


### STEP 4: RE-CREATE LH-AD GRAPHS AND HIGHLIGHT DIFFERENCES BETWEEN BAIT SETS AND TAXONOMIC CLADES

#For these graphs we use the 'strict' routine, i.e. we do not consider the samples for this data recovery was too poor to calculate trustworthy LH and AD values.
#Load config
config_file <- "results_HybPhaser_allGenes/config_strict.txt"
source(config_file)

#Set paths to results from script 1b_assess_dataset.R
if(name_for_dataset_optimization_subset != ""){
  folder_subset_add <- paste("_",name_for_dataset_optimization_subset, sep="")
}else {
  folder_subset_add <- ""
} 
output_Robjects <- file.path("results_HybPhaser_allGenes",path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)
output_assess <- file.path("results_HybPhaser_allGenes",path_to_output_folder,paste("02_assessment", folder_subset_add, sep=""))

#Load table with locus heterozygosity and allele divergence data
tab_het_ad <- readRDS(file=file.path(output_Robjects,"Summary_table.Rds"))

#Add taxonomic details for use in plotting
tab_het_ad$Tribe <- metadata$Tribe[match(tab_het_ad$sample, metadata$Library_ID)]
tab_het_ad$Genus <- metadata$Genus[match(tab_het_ad$sample, metadata$Library_ID)]
tab_het_ad$Species <- metadata$Name[match(tab_het_ad$sample, metadata$Library_ID)]

#Create a version for LH and AD of the dataframe in long format for plotting by bait set in ggplot.
tab_het_ad_long <- rbind(data.frame(Sample = tab_het_ad$sample, Tribe = tab_het_ad$Tribe, Species = tab_het_ad$Species, Bait_set = "B764", allele_divergence_bait_set = tab_het_ad$allele_divergence_Nikolov2019, locus_heterozygosity_bait_set = tab_het_ad$locus_heterozygosity_Nikolov2019),
                         data.frame(Sample = tab_het_ad$sample, Tribe = tab_het_ad$Tribe, Species = tab_het_ad$Species, Bait_set = "A353", allele_divergence_bait_set = tab_het_ad$allele_divergence_Angiosperms353, locus_heterozygosity_bait_set = tab_het_ad$locus_heterozygosity_Angiosperms353))
  
#Calculate sample sizes for tribes to be printed over boxplots
library(dplyr)
tab_het_ad_tribe_totals <- tab_het_ad %>%
  count(Tribe)
tab_het_ad_tribe_totals_bait_set <- tab_het_ad_long %>%
  count(Tribe, Bait_set)

#Create boxplots for LH and AD and sort by tribe to understand what tribes have high values and need special attention
LH_boxplot <- ggplot(data = tab_het_ad, aes(x = reorder(Tribe, -locus_heterozygosity, FUN = median), y = locus_heterozygosity))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(data = tab_het_ad_tribe_totals, aes(x = Tribe, y = 120, label = n))+
  xlab("Tribe")+
  ylab("Locus heterozygosity")

AD_boxplot <- ggplot(data = tab_het_ad, aes(x = reorder(Tribe, -allele_divergence, FUN = median), y = allele_divergence))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(data = tab_het_ad_tribe_totals, aes(x = Tribe, y = 11, label = n))+
  xlab("Tribe")+
  ylab("Allele divergence")

library(cowplot)
LH_and_AD_boxplots <- plot_grid(LH_boxplot, AD_boxplot, nrow = 2)
ggsave(file = "2h.hybphaser_result_LH_and_AD_boxplots.pdf", LH_and_AD_boxplots, width = 12, height = 10)
ggsave(file = "2h.hybphaser_result_LH_and_AD_boxplots.png", LH_and_AD_boxplots, width = 12, height = 10)


#We learn that some apparent outlier tribes are the Brassicaceae, Thelypodieae, Lepidieae and Stevenieae.
#Create a new variable for plotting purposes
tab_het_ad$Tribe_outlier <- tab_het_ad$Tribe
tab_het_ad$Tribe_outlier[!tab_het_ad$Tribe_outlier %in% c("Brassiceae", "Alysseae", "Anastatiaceae", "Cardamineae", "Thelypodieae", "Lepidieae", "Microlepidieae", "Arabideae")] <- "Other_tribe"
#Also add a column to later highlight the rogue taxa in a plot.
tab_het_ad$Rogue_taxa <- ""
tab_het_ad$Rogue_taxa[tab_het_ad$Tribe %in% c("Anastaticeae", "Biscutelleae", "Cochlearieae", "Megacarpaeeae")] <- tab_het_ad$Tribe[tab_het_ad$Tribe %in% c("Anastaticeae", "Biscutelleae", "Cochlearieae", "Megacarpaeeae")]
tab_het_ad$Rogue_taxa[tab_het_ad$Genus %in% c("Iberis", "Teesdalia")] <- tab_het_ad$Genus[tab_het_ad$Genus %in% c("Iberis", "Teesdalia")]

#Add another variable for highlighting the larger tribes
large_tribes <- data.frame(tab_het_ad %>%
  count(Tribe))
large_tribes <- large_tribes[order(-large_tribes$n),]
tab_het_ad$Tribe_large <- tab_het_ad$Tribe
tab_het_ad$Tribe_large[!tab_het_ad$Tribe_large %in% large_tribes$Tribe[1:28]] <- "Other_tribe"

#Add a variable for annotation with species name.
tab_het_ad$annotate <- tab_het_ad$Species
tab_het_ad$annotate[!tab_het_ad$Tribe %in% large_tribes$Tribe[1:28]]<-paste0(tab_het_ad$Species[!tab_het_ad$Tribe %in% large_tribes$Tribe[1:28]]," (",
                                                                             tab_het_ad$Tribe[!tab_het_ad$Tribe %in% large_tribes$Tribe[1:28]],")")
#Plot using facets by tribe.
(LH_AD_plot <- ggplot()+
    annotate("rect", xmin=1, xmax=4, ymin=80, ymax=Inf, fill="#00AFBB", alpha=0.5)+ #Hybrids?
    annotate("rect", xmin=4, xmax=Inf, ymin=80, ymax=Inf, fill="#E7B800", alpha=0.5)+ #Highly polyploids?
    annotate("rect", xmin=1, xmax=4, ymin=60, ymax=80, fill="mediumseagreen", alpha=0.5)+ #Old polyploids?
    annotate("rect", xmin=4, xmax=Inf, ymin=60, ymax=80, fill="#FC4E07", alpha=0.5)+ #Old & highly polyploids?
    geom_point(data = tab_het_ad, aes(x = allele_divergence, y = locus_heterozygosity), size = 0.1, color="black")+
    geom_text(data = tab_het_ad, aes(label = annotate, x = allele_divergence, y = locus_heterozygosity), size = 0.7, hjust = 0, vjust = 0, position=position_jitter())+
    theme_bw()+
    facet_wrap(. ~ Tribe_large, ncol=5)+
    xlab("Allele divergence")+
    ylab("Locus heterozygosity"))

ggsave(file = "2i.hybphaser_result_LH_vs_AD_plot.pdf", LH_AD_plot, width = 7, height = 7)
ggsave(file = "2i.hybphaser_result_LH_vs_AD_plot.png", LH_AD_plot, width = 7, height = 7)

#Plot showing possible differences by bait set (but there are none).
(LH_AD_plot_bait_set <- ggplot(data = tab_het_ad_long)+
    annotate("rect", xmin=1, xmax=4, ymin=80, ymax=Inf, fill="#00AFBB", alpha=0.5)+ #Hybrids?
    annotate("rect", xmin=4, xmax=Inf, ymin=80, ymax=Inf, fill="#E7B800", alpha=0.5)+ #Highly polyploids?
    annotate("rect", xmin=1, xmax=4, ymin=60, ymax=80, fill="mediumseagreen", alpha=0.5)+ #Old polyploids?
    annotate("rect", xmin=4, xmax=Inf, ymin=60, ymax=80, fill="#FC4E07", alpha=0.5)+ #Old & highly polyploids?
    geom_point(aes(x = allele_divergence_bait_set, y = locus_heterozygosity_bait_set, color = Bait_set), size = 0.3)+
    theme_bw()+
    scale_color_manual(values = c("black","grey"))+
    xlab("Allele divergence")+
    ylab("Locus heterozygosity"))

ggsave(file = "2j.hybphaser_result_LH_vs_AD_plot_by_bait_set.pdf", LH_AD_plot_bait_set, width = 12, height = 8)

#Plot for all data together, now using different shapes to highlight most extreme tribes.
(LH_AD_plot_main <- ggplot()+
    annotate("rect", xmin=1, xmax=4, ymin=80, ymax=Inf, fill="#00AFBB", alpha=0.5)+ #Hybrids?
    annotate("rect", xmin=4, xmax=Inf, ymin=80, ymax=Inf, fill="#E7B800", alpha=0.5)+ #Highly polyploids?
    annotate("rect", xmin=1, xmax=4, ymin=60, ymax=80, fill="mediumseagreen", alpha=0.5)+ #Old polyploids?
    annotate("rect", xmin=4, xmax=Inf, ymin=60, ymax=80, fill="#FC4E07", alpha=0.5)+ #Old & highly polyploids?
    # geom_point(data = tab_het_ad[tab_het_ad$Tribe_outlier!="Other_tribe",], aes(x = allele_divergence, y = locus_heterozygosity, shape=Tribe_outlier), size=1, colour="red")+
    # scale_shape_manual(values=c(0, 8, 1, 2, 5,6, 9, 13))+
    stat_ellipse(data = tab_het_ad[tab_het_ad$Tribe_outlier!="Other_tribe",], type = 't', geom = "polygon", aes(x = allele_divergence, y = locus_heterozygosity, linetype = Tribe_outlier), fill="grey", alpha=0, size = 0.5, colour = "black")+
    scale_linetype_manual(name="Selected tribes", values=c("dotdash", "blank", "longdash", "dotted", "dashed", "twodash", "solid"))+
    geom_point(data = tab_het_ad, aes(x = allele_divergence, y = locus_heterozygosity), size=1, shape=20)+
    # geom_text(data = tab_het_ad, aes(label = annotate, x = allele_divergence, y = locus_heterozygosity), size = 0.7, hjust = 0, vjust = 0, position=position_jitter())+
    geom_point(data = tab_het_ad[tab_het_ad$Rogue_taxa != "", ], aes(x = allele_divergence, y = locus_heterozygosity, shape=Tribe_outlier), size=1, shape=20)+
    theme_bw()+
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 120)) +
    xlab("Allele divergence [%]")+
    ylab("Locus heterozygosity [%]"))

ggsave(file = "2k.hybphaser_result_LH_vs_AD_plot_main.pdf", LH_AD_plot_main, width = 6.5, height = 4.5)

#Summarise LH and AD values by tribe for publication in a table.
library(dplyr)
AD_and_LH_tribe_summary<-tab_het_ad %>%
  group_by(Tribe) %>%
  summarise(n=n(),
            LH_mean_all_genes=mean(locus_heterozygosity, na.rm = T),
            LH_min_all_genes=min(locus_heterozygosity, na.rm = T),
            LH_max_all_genes=max(locus_heterozygosity, na.rm = T),
            LH_mean_B764=mean(locus_heterozygosity_Nikolov2019, na.rm = T),
            LH_min_B764=min(locus_heterozygosity_Nikolov2019, na.rm = T),
            LH_max_B764=max(locus_heterozygosity_Nikolov2019, na.rm = T),
            LH_mean_A353=mean(locus_heterozygosity_Angiosperms353, na.rm = T),
            LH_min_A353=min(locus_heterozygosity_Angiosperms353, na.rm = T),
            LH_max_A353=max(locus_heterozygosity_Angiosperms353, na.rm = T),
            AD_mean_all_genes=mean(allele_divergence, na.rm = T),
            AD_min_all_genes=min(allele_divergence, na.rm = T),
            AD_max_all_genes=max(allele_divergence, na.rm = T),
            AD_mean_B764=mean(allele_divergence_Nikolov2019, na.rm = T),
            AD_min_B764=min(allele_divergence_Nikolov2019, na.rm = T),
            AD_max_B764=max(allele_divergence_Nikolov2019, na.rm = T),
            AD_mean_A353=mean(allele_divergence_Angiosperms353, na.rm = T),
            AD_min_A353=min(allele_divergence_Angiosperms353, na.rm = T),
            AD_max_A353=max(allele_divergence_Angiosperms353, na.rm = T)) %>%
  arrange(desc(LH_mean_all_genes))

#Update table for print.
AD_and_LH_tribe_summary_upd<-data.frame(Tribe = AD_and_LH_tribe_summary$Tribe,
                                        n = AD_and_LH_tribe_summary$n,
                                        LH_all_genes = paste0(format(round(AD_and_LH_tribe_summary$LH_mean_all_genes, 1), nsmall = 1), " (",
                                                              paste0(format(round(AD_and_LH_tribe_summary$LH_min_all_genes, 1), nsmall = 1)), " - ",
                                                              paste0(format(round(AD_and_LH_tribe_summary$LH_max_all_genes, 1), nsmall = 1)), ")"),
                                        LH_B764_genes = paste0(format(round(AD_and_LH_tribe_summary$LH_mean_B764, 1), nsmall = 1), " (",
                                                               paste0(format(round(AD_and_LH_tribe_summary$LH_min_B764, 1), nsmall = 1)), " - ",
                                                               paste0(format(round(AD_and_LH_tribe_summary$LH_max_B764, 1), nsmall = 1)), ")"),
                                        LH_A353_genes = paste0(format(round(AD_and_LH_tribe_summary$LH_mean_A353, 1), nsmall = 1), " (",
                                                               paste0(format(round(AD_and_LH_tribe_summary$LH_min_A353, 1), nsmall = 1)), " - ",
                                                               paste0(format(round(AD_and_LH_tribe_summary$LH_max_A353, 1), nsmall = 1)), ")"),
                                        AD_all_genes = paste0(format(round(AD_and_LH_tribe_summary$AD_mean_all_genes, 1), nsmall = 1), " (",
                                                              paste0(format(round(AD_and_LH_tribe_summary$AD_min_all_genes, 1), nsmall = 1)), " - ",
                                                              paste0(format(round(AD_and_LH_tribe_summary$AD_max_all_genes, 1), nsmall = 1)), ")"),
                                        AD_B764_genes = paste0(format(round(AD_and_LH_tribe_summary$AD_mean_B764, 1), nsmall = 1), " (",
                                                               paste0(format(round(AD_and_LH_tribe_summary$AD_min_B764, 1), nsmall = 1)), " - ",
                                                               paste0(format(round(AD_and_LH_tribe_summary$AD_max_B764, 1), nsmall = 1)), ")"),
                                        AD_A353_genes = paste0(format(round(AD_and_LH_tribe_summary$AD_mean_A353, 1), nsmall = 1), " (",
                                                               paste0(format(round(AD_and_LH_tribe_summary$AD_min_A353, 1), nsmall = 1)), " - ",
                                                               paste0(format(round(AD_and_LH_tribe_summary$AD_max_A353, 1), nsmall = 1)), ")"))
AD_and_LH_tribe_summary_upd<- AD_and_LH_tribe_summary_upd[!grepl("Outgroup", AD_and_LH_tribe_summary_upd$Tribe, fixed = T),]

#Write results to Excel for export and print in publication.
library(openxlsx)
write.xlsx(AD_and_LH_tribe_summary_upd, '2m.hybphaser_results_LH_and_AD_summary_by_tribe.xlsx')


### STEP 5: CREATE AN OVERVIEW OF LIKELY HYBRIDS AND ROGUE TAXA TO BE EXCLUDED FROM ADDITIONAL BRASSITOL ANALYSIS IN NEXT STEP

#Add details on rogue taxa for subsetting later.
rogue_tribes<-c("Anastaticeae", "Biscutelleae", "Cochlearieae", "Iberideae I", "Iberideae II", "Megacarpaeeae")
rogue_samples<-tab_het_ad$sample[tab_het_ad$Tribe %in% rogue_tribes]

#Add details on likely hybrids for subsetting later.
#Although rather arbitrary, we will select all samples in the upper 50% of LH and AD ranges.
sample_hybrid_LH<-tab_het_ad$sample[tab_het_ad$locus_heterozygosity > median(tab_het_ad$locus_heterozygosity)]
sample_hybrid_AD<-tab_het_ad$sample[tab_het_ad$allele_divergence > median(tab_het_ad$allele_divergence)]

#Create a vector which is the sum of the above vectors.
samples_to_remove<-sort(unique(c(rogue_samples, sample_hybrid_LH, sample_hybrid_AD)))
#Note that this means we'll remove over half of all samples in our later analyses excluding possible hybrids and rogue taxa.

#Export this vector for later use.
write.table(samples_to_remove, file = "2n.list_of_hybrids_and_rogue_taxa.txt",
            row.names = F,
            col.names= F,
            sep="\t",
            quote = FALSE)
