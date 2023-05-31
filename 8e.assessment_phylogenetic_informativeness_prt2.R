#This script collects rate4site results and creates phylogenetic informativeness profiles for all genes.

#Run the following lines once the first time to install the correct version of PhyInformR.
#library(devtools)
#install_github("carolinafishes/PhyInformR")

#Load packages.
library(ape)
library(phytools)
library(tidyr)
library(PhyInformR)
library(ggplot2)
library(dplyr)
library(forcats)

###MAKE SURE THE BASH SCRIPT PIPES THE GENELIST TO THIS R SCRIPT LATER.
genelist<-read.csv("8c.genelist.txt", header = F)

#Create a list to store all results.
results_site_rates<-list()

#Load the results from the tool rate4site.
for (i in 1:nrow(genelist)){
#for (i in 1:4){
  gene<-genelist[i,]
  #Check if gene file exits, otherwise move on to the next gene in the list.
  if(file.exists(paste0("results_assessment_phylogenetic_informativeness/results_rate4site/",gene,"_result_rate4site.res"))){
    #Load the results file.
    result<-read.table(file = paste0("results_assessment_phylogenetic_informativeness/results_rate4site/",gene,"_result_rate4site.res"))
    result_list<-list(result[,3])
    names(result_list)<-as.character(gene)
    #Collect the site rates and store in the list.
    results_site_rates<-append(results_site_rates, result_list)
  }
}

#Load the species tree that we want to use to make the assessment of phylogenetic informativeness for each gene.
#We chose to work with the ML tree based on a supermatrix of strict genes.
species_tree<-read.tree(file = "results_assessment_phylogenetic_informativeness/input_species_tree/9b.IQ-TREE_supermatrix_approach.tree")
#Root the tree to the outgroup.
species_tree<-root(species_tree, outgroup = "PAFTOL_019361")
#Make the tree ultrametric.
species_tree <- chronos(species_tree)
#Remove bootstrap values from the tree.
species_tree$node.label<-NULL

#Create a dataframe in which we can store the phylogenetic informativeness results for plotting.
phylogenetic_informativeness_profiles<-data.frame(time=NA,
                                                  PI=NA,
                                                  gene=NA,
                                                  bait_set=NA)

#We calculate phylogenetic informativeness using the functions that come with package PhyInformR, which we extracted for our goals.
#We loop through the genes in our site rates result list.
for (j in 1:length(results_site_rates)){
  gene<-names(results_site_rates)[j]
  #Set correct input data for this gene.
  tree<-species_tree
  rate.vector<-as.matrix(results_site_rates[[j]])
  #Run PhyInformR function informativeness.profile(.
  btimes <- branching.times(tree)
  btimes2 <- c(0, btimes)
  #Run hidden PhyInformR function inform.profile.generator.
  sorted.btimes <- sort(btimes2)
  branching.points <- length(btimes2)
  calculation.length <- length(rate.vector)
  inform.at.time <- matrix(ncol = branching.points)
  for (k in 1:branching.points) {
    btime <- sorted.btimes[k]
    #Run hidden PhyInformR function site.summer.
    at.site <- matrix(ncol = calculation.length)
    for (m in 1:calculation.length) {
      current <- rate.vector[m]
      at.site[m] <- 16 * current * current * btime * exp(-4 * current * btime)
    }
    inform.at.time[k] <- sum(at.site)
  }
  close <- rbind(sorted.btimes, inform.at.time)
  #Store results in the dataframe for plotting later.
  PI_results<-NULL
  PI_results<-data.frame(time=as.numeric(sorted.btimes),
                         PI=as.numeric(inform.at.time),
                         gene=as.character(gene),
                         bait_set=if (nchar(gene)<5) {"A353"} else {"B764"})
  phylogenetic_informativeness_profiles<-rbind(phylogenetic_informativeness_profiles, PI_results)
}

#Remove NA values from dataframe.
phylogenetic_informativeness_profiles<-phylogenetic_informativeness_profiles[!(is.na(phylogenetic_informativeness_profiles$time)),]

#We want to make a distinction between genes with a mean SNP proportion <0.02 and >0.02; the latter were highlighted by HybPhaser step 1b in the 'superstrict' routine.
#Collect results from the HybPhaser SNP table saved as R object and take the mean for each gene.
HybPhaser_Table_SNPs<-readRDS("results_HybPhaser_allGenes/output_HybPhaser_allGenes/00_R_objects/inclusive/Table_SNPs.Rds")
genes_high_SNP_prop<-rownames(HybPhaser_Table_SNPs)[(rowMeans(as.data.frame(HybPhaser_Table_SNPs), na.rm = T)>0.02)]

#Update the dataframe.
phylogenetic_informativeness_profiles$mean_SNP_prop<-NA
phylogenetic_informativeness_profiles$mean_SNP_prop[phylogenetic_informativeness_profiles$gene %in% genes_high_SNP_prop]<-">0.02"
phylogenetic_informativeness_profiles$mean_SNP_prop[!phylogenetic_informativeness_profiles$gene %in% genes_high_SNP_prop]<-"<=0.02"

#We also want to know the length of the target gene alignment used; short genes will be unable to convey the information properly and may need to be removed.
phylogenetic_informativeness_profiles$target_length<-NA
for (k in 1:length(unique(phylogenetic_informativeness_profiles$gene))){
  #Read gene sequence alignment and length.
  gene_k<-unique(phylogenetic_informativeness_profiles$gene)[k]
  gene_alignment_k<-read.FASTA(paste0("results_assessment_phylogenetic_informativeness/input_gene_alignments/", gene_k, "_taper_final_inclusive.fasta"))
  gene_alignment_length_k<-length(gene_alignment_k[[1]])
  #Add result to dataframe.
  phylogenetic_informativeness_profiles[phylogenetic_informativeness_profiles$gene==gene_k,]$target_length<-gene_alignment_length_k
}

#Plot the results (exclude alignments <100 bp).
phylogenetic_informativeness_profiles_plot<-ggplot(data=phylogenetic_informativeness_profiles[phylogenetic_informativeness_profiles$target_length>99,],aes(x=time, y=PI, group=gene, colour=fct_relevel(mean_SNP_prop,"<=0.02",">0.02"), size=fct_relevel(mean_SNP_prop,"<=0.02",">0.02")))+
  geom_line()+
  scale_size_manual(values=c(0.15, 0.05), guide = "none")+
  scale_colour_manual(name="Mean SNP proportion", values=c("black", "red"))+
  scale_x_reverse()+
  theme_classic()+
  xlab("Relative time")+
  ylab("Townsend's phylogenetic informativeness")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  guides(linetype=guide_legend(title="Bait set"))+
  facet_grid(cols = vars(bait_set))

#Save the plot to a file.
ggsave(phylogenetic_informativeness_profiles_plot, filename = "8f.phylogenetic_informativeness_profiles_plot.png",
       width = 8,
       height = 6)

#Finally, let's calculate group means for quick & dirty comparison.
phylogenetic_informativeness_means<-data.frame(phylogenetic_informativeness_profiles[phylogenetic_informativeness_profiles$target_length>99,] %>%
  group_by(bait_set, mean_SNP_prop) %>%
  summarize(Mean = mean(PI),
            S.D. = sd(PI),
            n = n()/branching.points))

#Save this result as a table.
write.csv(phylogenetic_informativeness_means, file = "8g.phylogenetic_informativeness_means_table.csv")







