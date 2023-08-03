###This script collects and summarizes results from previous analyses by HybPhaser, ASTRAL-III, ASTRAL-Pro, and IQ-TREE.

#LOAD PACKAGES ----
library(ape)
library(treeio)
library(tidytree)
library(devtools); devtools::install_github("chutter/AstralPlane")
library(AstralPlane)
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2); library(ggtree)
library(ggimage)
library(ggnewscale)
library(cowplot)
library(colorRamps)
library(phytools)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(TreeDist)
library(phangorn)


#READ IN METADATA ----
metadata<-read.csv("2e.metadata.csv", header=T, stringsAsFactors=FALSE, fileEncoding="latin1")


#READ IN PHYLOGENETIC RESULTS ----

####TEMP NEW DATA CLOCKLIKE GENES
tree_ASTRALIII_superstrict_and_clocklike<-read.tree(file = "results_ASTRALIII_superstrict_and_clocklike/species_tree.tre")
tree_ASTRALIII_superstrict_and_clocklike<-root(tree_ASTRALIII_superstrict_and_clocklike, outgroup = "S1321")
tree_ASTRALIII_superstrict_and_clocklike_astral_annotations<-createAstralPlane(astral.tree = "results_ASTRALIII_superstrict_and_clocklike/species_tree_t2.tre",
                                                                               outgroups = "S1321",
                                                                               tip.length = 1)
tree_ASTRALIII_superstrict_and_clocklike_astral_annotations_piecharts<-nodepie(tree_ASTRALIII_superstrict_and_clocklike_astral_annotations@nodeData, cols=2:4)
tree_ASTRALIII_superstrict_and_clocklike_astral_annotations_piecharts<-lapply(tree_ASTRALIII_superstrict_and_clocklike_astral_annotations_piecharts, function(g) g+scale_fill_manual(values = c("darkblue","lightblue","grey")))

mean(tree_ASTRALIII_superstrict_and_clocklike_astral_annotations@nodeData$q1, na.rm = T)

## ASTRAL-III inclusive ----
tree_ASTRALIII_inclusive<-read.tree(file = "results_ASTRALIII_inclusive/species_tree.tre")
#Root the tree to the outgroup
tree_ASTRALIII_inclusive<-root(tree_ASTRALIII_inclusive, outgroup = "S1321")

#Import ASTRAL results to use node annotations with R package AstralPlane.
#Note that here we import the "-t 2" version from ASTRAL-III, in which all the annotations for node support have been logged.
tree_ASTRALIII_inclusive_astral_annotations<-createAstralPlane(astral.tree = "results_ASTRALIII_inclusive/species_tree_t2.tre",
                                                               outgroups = "S1321",
                                                               tip.length = 1)

#Use the nicely ordered dataframe from this object to create a list object with pie charts for plotting in ggtree later.
#Put pie chart data into a list that can be read by ggtree (see https://yulab-smu.top/treedata-book/chapter8.html).
tree_ASTRALIII_inclusive_astral_annotations_piecharts<-nodepie(tree_ASTRALIII_inclusive_astral_annotations@nodeData, cols=2:4)
tree_ASTRALIII_inclusive_astral_annotations_piecharts<-lapply(tree_ASTRALIII_inclusive_astral_annotations_piecharts, function(g) g+scale_fill_manual(values = c("darkblue","lightblue","grey")))

#Also read in the set of gene trees used by ASTRAL.
tree_ASTRALIII_inclusive_gene_trees<-read.tree(file = "results_ASTRALIII_inclusive/gene_trees_combined.tre")


## ASTRAL-III strict ----
tree_ASTRALIII_strict<-read.tree(file = "results_ASTRALIII_strict/species_tree.tre")
#Root the tree to the outgroup
tree_ASTRALIII_strict<-root(tree_ASTRALIII_strict, outgroup = "PAFTOL_019361")

#Import ASTRAL results to use node annotations with R package AstralPlane.
#Note that here we import the "-t 2" version from ASTRAL-III, in which all the annotations for node support have been logged.
tree_ASTRALIII_strict_astral_annotations<-createAstralPlane(astral.tree = "results_ASTRALIII_strict/species_tree_t2.tre",
                                                            outgroups = "PAFTOL_019361",
                                                            tip.length = 1)

#Use the nicely ordered dataframe from this object to create a list object with pie charts for plotting in ggtree later.
#Put pie chart data into a list that can be read by ggtree (see https://yulab-smu.top/treedata-book/chapter8.html).
tree_ASTRALIII_strict_astral_annotations_piecharts<-nodepie(tree_ASTRALIII_strict_astral_annotations@nodeData, cols=2:4)
tree_ASTRALIII_strict_astral_annotations_piecharts<-lapply(tree_ASTRALIII_strict_astral_annotations_piecharts, function(g) g+scale_fill_manual(values = c("darkblue","lightblue","grey")))

#Also read in the set of gene trees used by ASTRAL.
tree_ASTRALIII_strict_gene_tree<-read.tree(file = "results_ASTRALIII_strict/gene_trees_combined.tre")


## ASTRAL-III superstrict ----
tree_ASTRALIII_superstrict<-read.tree(file = "results_ASTRALIII_superstrict/species_tree.tre")
#Root the tree to the outgroup
tree_ASTRALIII_superstrict<-root(tree_ASTRALIII_superstrict, outgroup = "PAFTOL_019361")

#Import ASTRAL results to use node annotations with R package AstralPlane.
#Note that here we import the "-t 2" version from ASTRAL-III, in which all the annotations for node support have been logged.
tree_ASTRALIII_superstrict_astral_annotations<-createAstralPlane(astral.tree = "results_ASTRALIII_superstrict/species_tree_t2.tre",
                                                                 outgroups = "PAFTOL_019361",
                                                                 tip.length = 1)

#Use the nicely ordered dataframe from this object to create a list object with pie charts for plotting in ggtree later.
#Put pie chart data into a list that can be read by ggtree (see https://yulab-smu.top/treedata-book/chapter8.html).
tree_ASTRALIII_superstrict_astral_annotations_piecharts<-nodepie(tree_ASTRALIII_superstrict_astral_annotations@nodeData, cols=2:4)
tree_ASTRALIII_superstrict_astral_annotations_piecharts<-lapply(tree_ASTRALIII_superstrict_astral_annotations_piecharts, function(g) g+scale_fill_manual(values = c("darkblue","lightblue","grey")))

#Also read in the set of gene trees used by ASTRAL.
tree_ASTRALIII_superstrict_gene_trees<-read.tree(file = "results_ASTRALIII_superstrict/gene_trees_combined.tre")

## ASTRAL-III superstrict-by-tribe ----
tree_ASTRALIII_superstrict_by_tribe<-read.tree(file = "results_ASTRALIII_superstrict_by_tribe/species_tree.tre")
#Root the tree to the outgroup
tree_ASTRALIII_superstrict_by_tribe<-root(tree_ASTRALIII_superstrict_by_tribe, outgroup = "PAFTOL_019361")

#Import ASTRAL results to use node annotations with R package AstralPlane.
#Note that here we import the "-t 2" version from ASTRAL-III, in which all the annotations for node support have been logged.
tree_ASTRALIII_superstrict_by_tribe_astral_annotations<-createAstralPlane(astral.tree = "results_ASTRALIII_superstrict_by_tribe/species_tree_t2.tre",
                                                                          outgroups = "PAFTOL_019361",
                                                                          tip.length = 1)

#Use the nicely ordered dataframe from this object to create a list object with pie charts for plotting in ggtree later.
#Put pie chart data into a list that can be read by ggtree (see https://yulab-smu.top/treedata-book/chapter8.html).
tree_ASTRALIII_superstrict_by_tribe_astral_annotations_piecharts<-nodepie(tree_ASTRALIII_superstrict_by_tribe_astral_annotations@nodeData, cols=2:4)
tree_ASTRALIII_superstrict_by_tribe_astral_annotations_piecharts<-lapply(tree_ASTRALIII_superstrict_by_tribe_astral_annotations_piecharts, function(g) g+scale_fill_manual(values = c("darkblue","lightblue","grey")))

#Also read in the set of gene trees used by ASTRAL.
tree_ASTRALIII_superstrict_by_tribe_gene_trees<-read.tree(file = "results_ASTRALIII_superstrict_by_tribe/gene_trees_combined.tre")

## ASTRAL-Pro ----
tree_ASTRAL_Pro<-read.tree(file = "6i.ASTRAL_Pro_species_tree.tre")
#Root the tree to the outgroup
tree_ASTRAL_Pro<-root(tree_ASTRAL_Pro, outgroup = "S1321")

#Import ASTRAL results to use node annotations with R package AstralPlane.
#Note that here we import the "-t 2" version from ASTRAL-III, in which all the annotations for node support have been logged.
tree_ASTRAL_Pro_astral_annotations<-createAstralPlane(astral.tree = "6i.ASTRAL_Pro_species_tree.tre",
                                                      outgroups = "S1321",
                                                      tip.length = 1)

#Use the nicely ordered dataframe from this object to create a list object with pie charts for plotting in ggtree later.
#Put pie chart data into a list that can be read by ggtree (see https://yulab-smu.top/treedata-book/chapter8.html).
tree_ASTRAL_Pro_astral_annotations_piecharts<-nodepie(tree_ASTRAL_Pro_astral_annotations@nodeData, cols=2:4)
tree_ASTRAL_Pro_astral_annotations_piecharts<-lapply(tree_ASTRAL_Pro_astral_annotations_piecharts, function(g) g+scale_fill_manual(values = c("darkblue","lightblue","grey")))

#Also read in the set of gene trees used by ASTRAL.
tree_ASTRAL_Pro_gene_trees<-read.tree(file = "6f.ASTRAL_Pro_input_gene_trees.tre")


##ML nuclear ----

#Read in the ML IQ-TREE that was calibrated in treePL. This tree was already rooted to S1321 in IQ-TREE.
tree_IQ_TREE_supermatrix_calibrated<-read.beast(file = "results_treePL_calibration/output_treePL_dating/treePL_calibrated_phylogeny_summary.tre")

#Also read in the IQ-TREE version of this same tree with cGF and sCF concordance factors.
tree_IQ_TREE_supermatrix_concordance_factors<-read.tree(file = "results_treePL_calibration/output_iqtree_condordance_factors/iqtree_ML_tree_with_concordance_factors.tree")
#Add the concordance factors to the calibrated tree, such that all tree node annotations are held within a single tree object.
tree_IQ_TREE_supermatrix_calibrated@data$node_cf<-NA
tree_IQ_TREE_supermatrix_calibrated@data$gCF<-NA
tree_IQ_TREE_supermatrix_calibrated@data$sCF<-NA
tree_IQ_TREE_supermatrix_calibrated@data$concordance_factors<-NA
tree_IQ_TREE_supermatrix_calibrated@data$bootstrap<-NA
tree_IQ_TREE_supermatrix_calibrated@data$bootstrap_symbol<-NA
tree_IQ_TREE_supermatrix_calibrated@data$support<-NA
#And eight more columns to score gCF and sCF values in bins, which will be useful for plotting node labels later.
tree_IQ_TREE_supermatrix_calibrated@data$sCF_support<-NA
tree_IQ_TREE_supermatrix_calibrated@data$gCF_support<-NA
#Use phytools function matchNodes to find correspondence among node numbers in both trees.
matched_nodes<-data.frame(matchNodes(tree_IQ_TREE_supermatrix_calibrated@phylo, tree_IQ_TREE_supermatrix_concordance_factors, method = "descendants"))
tree_IQ_TREE_supermatrix_calibrated@data$node_cf<-matched_nodes$tr2[match(tree_IQ_TREE_supermatrix_calibrated@data$node, matched_nodes$tr1)]
#Loop through nodes to find concordance factors and store results.
for (i in 1:nrow(tree_IQ_TREE_supermatrix_calibrated@data)){
  node_cf<-tree_IQ_TREE_supermatrix_calibrated@data$node_cf[i]
  if(is.na(node_cf)==F){
    #Find concordance factors.
    concordance_factor<-tree_IQ_TREE_supermatrix_concordance_factors$node.label[node_cf-length(tree_IQ_TREE_supermatrix_concordance_factors$tip.label)]
    #Update notation.
    bootstrap<-as.numeric(strsplit(concordance_factor, "/")[[1]][1])
    bootstrap_symbol<-as.numeric(strsplit(concordance_factor, "/")[[1]][1])
    if(is.na(bootstrap_symbol)==F){if(round(bootstrap_symbol, 0)==100){bootstrap_symbol<-"*"}}
    gCF<-as.numeric(strsplit(concordance_factor, "/")[[1]][3])
    sCF<-as.numeric(strsplit(concordance_factor, "/")[[1]][4])
    #Store results.
    tree_IQ_TREE_supermatrix_calibrated@data$bootstrap[i]<-bootstrap
    tree_IQ_TREE_supermatrix_calibrated@data$bootstrap_symbol[i]<-bootstrap_symbol
    tree_IQ_TREE_supermatrix_calibrated@data$gCF[i]<-gCF
    tree_IQ_TREE_supermatrix_calibrated@data$sCF[i]<-sCF
    if(is.na(bootstrap)==F){if(round(bootstrap, 0)==100){bootstrap<-"*"}}
    #Define the support from these values to colour labels in tree plot later (rules based on suggestion by Lars Nauheimer).
    if(is.na(sCF)==F | is.na(gCF)==F){
      if(sCF>40){support="blue"}
      if(gCF>20){support="black"}
      if(sCF<40){support="orange"}
      if(sCF<35){support="red"}
      concordance_factor_combined<-paste(round(as.numeric(gCF), 0),
                                         round(as.numeric(sCF), 0),
                                         sep="/")
      #Store results.
      tree_IQ_TREE_supermatrix_calibrated@data$concordance_factors[i]<-concordance_factor_combined
      tree_IQ_TREE_supermatrix_calibrated@data$support[i]<-support
      #Also define support from gCF and sCF for plotting node support later using the bins defined above.
      if(sCF>50) {
        tree_IQ_TREE_supermatrix_calibrated@data$sCF_support[i]<-"black"
        } else if (sCF>40) {
        tree_IQ_TREE_supermatrix_calibrated@data$sCF_support[i]<-"blue"
        } else if(sCF>30) {
          tree_IQ_TREE_supermatrix_calibrated@data$sCF_support[i]<-"orange"
        } else {
            tree_IQ_TREE_supermatrix_calibrated@data$sCF_support[i]<-"red"
            }
      if(gCF>30) {
        tree_IQ_TREE_supermatrix_calibrated@data$gCF_support[i]<-"black"
      } else if (gCF>25) {
        tree_IQ_TREE_supermatrix_calibrated@data$gCF_support[i]<-"blue"
      } else if(gCF>20) {
        tree_IQ_TREE_supermatrix_calibrated@data$gCF_support[i]<-"orange"
      } else {
        tree_IQ_TREE_supermatrix_calibrated@data$gCF_support[i]<-"red"
      }
    }
  }
}

#Check from a plot that these categories of support are realistic and make sense.
plot(cbind(tree_IQ_TREE_supermatrix_calibrated_ingroups@data$gCF,tree_IQ_TREE_supermatrix_calibrated_ingroups@data$sCF))

#(Re)rooting the tree is not needed in this case, because the tree was already rooted to S1321 in IQ-TREE.


##ML nuclear excluding hybrids ----

#Read in the IQ-TREE version of this same tree (but excluding hybrids) with cGF and sCF concordance factors.
tree_IQ_TREE_supermatrix_excluding_hybrids_concordance_factors<-read.tree(file = "9b.IQ-TREE_supermatrix_supermatrix_excluding_hybrids_with_concordance_factors.tree")
#Add the concordance factors to the calibrated tree, such that all tree node annotations are held within a single tree object.


##ML plastome ----

#Read in the ML IQ-TREE that was calibrated in treePL. This tree was already rooted to S1321 in IQ-TREE.
tree_IQ_TREE_chloroplast_calibrated<-read.beast(file = "7.chloroplast_treePL_BS_analysis_TreeAnnotator_EDITED.tre")

#Also read in the IQ-TREE version of this same tree with BS and sCF concordance factors.
tree_IQ_TREE_chloroplast_concordance_factors<-read.tree(file = "7.chloroplast.concord.cf.EDITED.tree")
#Add the concordance factors to the calibrated tree, such that all tree node annotations are held within a single tree object.
tree_IQ_TREE_chloroplast_calibrated@data$node_cf<-NA
tree_IQ_TREE_chloroplast_calibrated@data$concordance_factors<-NA
tree_IQ_TREE_chloroplast_calibrated@data$bootstrap<-NA
tree_IQ_TREE_chloroplast_calibrated@data$bootstrap_symbol<-NA
tree_IQ_TREE_chloroplast_calibrated@data$sCF<-NA
tree_IQ_TREE_chloroplast_calibrated@data$sCF_support<-NA
#Use phytools function matchNodes to find correspondence among node numbers in both trees.
matched_nodes<-data.frame(matchNodes(tree_IQ_TREE_chloroplast_calibrated@phylo, tree_IQ_TREE_chloroplast_concordance_factors, method = "descendants"))
tree_IQ_TREE_chloroplast_calibrated@data$node_cf<-matched_nodes$tr2[match(tree_IQ_TREE_chloroplast_calibrated@data$node, matched_nodes$tr1)]
#Loop through nodes to find concordance factors and store results.
for (i in 1:nrow(tree_IQ_TREE_chloroplast_calibrated@data)){
  node_cf<-tree_IQ_TREE_chloroplast_calibrated@data$node_cf[i]
  if(is.na(node_cf)==F){
    #Find concordance factors.
    concordance_factor<-tree_IQ_TREE_chloroplast_concordance_factors$node.label[node_cf-length(tree_IQ_TREE_chloroplast_concordance_factors$tip.label)]
    #Update notation.
    sCF<-as.numeric(strsplit(concordance_factor, "/")[[1]][2])
    bootstrap<-as.numeric(strsplit(concordance_factor, "/")[[1]][1])
    bootstrap_symbol<-bootstrap
    if(is.na(bootstrap)==F){if(round(bootstrap, 0)==100){bootstrap_symbol<-"*"}}
    #Define the support from these values to colour labels in tree plot later (rules based on suggestion by Lars Nauheimer).
    if(is.na(sCF)==F){
      if(sCF>50){sCF_support="black"}
      if(sCF>40){sCF_support="blue"}
      if(sCF>30){sCF_support="orange"}
      if(sCF<=30){sCF_support="red"}
      concordance_factor_combined<-paste(round(as.numeric(sCF), 0),
                                         sep="/")
      #Store results.
      tree_IQ_TREE_chloroplast_calibrated@data$bootstrap[i]<-bootstrap
      tree_IQ_TREE_chloroplast_calibrated@data$bootstrap_symbol[i]<-bootstrap_symbol
      tree_IQ_TREE_chloroplast_calibrated@data$sCF[i]<-sCF
      tree_IQ_TREE_chloroplast_calibrated@data$concordance_factors[i]<-concordance_factor_combined
      tree_IQ_TREE_chloroplast_calibrated@data$sCF_support[i]<-sCF_support
    }
  }
}

# LOG ASTRAL NODE SUPPORT ----
#These are defined as "the local posterior Probability (PP) that the branch is in the species tree" (Sayyari et al 2016).
branch_support_compared<-rbind(cbind(tree_ASTRALIII_inclusive_astral_annotations@nodeData,sampling_routine = rep("ASTRAL-III inclusive", nrow(tree_ASTRALIII_inclusive_astral_annotations@nodeData))),
                               cbind(tree_ASTRALIII_strict_astral_annotations@nodeData,sampling_routine = rep("ASTRAL-III strict", nrow(tree_ASTRALIII_strict_astral_annotations@nodeData))),
                               cbind(tree_ASTRALIII_superstrict_astral_annotations@nodeData,sampling_routine = rep("ASTRAL-III superstrict", nrow(tree_ASTRALIII_superstrict_astral_annotations@nodeData))),
                               cbind(tree_ASTRALIII_superstrict_by_tribe_astral_annotations@nodeData,sampling_routine = rep("ASTRAL-III superstrict_by_tribe", nrow(tree_ASTRALIII_superstrict_by_tribe_astral_annotations@nodeData))),
                               cbind(tree_ASTRAL_Pro_astral_annotations@nodeData,sampling_routine = rep("ASTRAL-Pro", nrow(tree_ASTRAL_Pro_astral_annotations@nodeData))))


# ANNOTATIONS DATAFRAME ----
tip_labels<-unique(c(tree_ASTRALIII_inclusive$tip.label, 
                     tree_ASTRALIII_strict$tip.label, 
                     tree_ASTRALIII_superstrict$tip.label, 
                     tree_ASTRAL_Pro$tip.label, 
                     tree_IQ_TREE_supermatrix_concordance_factors$tip.label, 
                     tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label,
                     tree_IQ_TREE_supermatrix_excluding_hybrids_concordance_factors$tip.label))
# tip_labels<-tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label
annotations<-data.frame(Library_ID=tip_labels)

#Add taxonomic data to annotations dataframe.
annotations$tribe<-metadata$Tribe[match(annotations$Library_ID,metadata$Library_ID)]
annotations$genus<-metadata$Genus[match(annotations$Library_ID,metadata$Library_ID)]
list_of_Brassicaceae_ingroup_genera<-unique(annotations$genus[!grepl("Outgroup_", annotations$tribe)])
list_of_Brassicaceae_ingroup_tribes<-unique(annotations$tribe[!grepl("Outgroup_", annotations$tribe)])
annotations$species<-metadata$Name[match(annotations$Library_ID,metadata$Library_ID)]
annotations$Library_ID2<-metadata$Library_ID[match(annotations$Library_ID,metadata$Library_ID)]
#In case tip label data come from another source (such as in the case of the chloroplast tree), simply copy the Library_ID from that tip label.
annotations[is.na(annotations$Library_ID2),]$Library_ID2<-annotations[is.na(annotations$Library_ID2),]$Library_ID

#Add details of collection country for plot of map at the end.
annotations$country<-metadata$Country[match(annotations$Library_ID,metadata$Library_ID)]

##Sample numbers ASTRAL-III inclusive----
number_samples_ingroup_ASTRALIII_inclusive<-length(annotations$Library_ID[!grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_ASTRALIII_inclusive$tip.label])
number_samples_outgroup_ASTRALIII_inclusive<-length(annotations$Library_ID[grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_ASTRALIII_inclusive$tip.label])
number_species_ingroup_ASTRALIII_inclusive<-length(unique(annotations$species[annotations$Library_ID %in% annotations$Library_ID[!grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]]))-1 #Minus one for the two populations of Clausia aprica
number_species_outgroup_ASTRALIII_inclusive<-length(unique(annotations$species[annotations$Library_ID %in% annotations$Library_ID[grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]]))
number_genera_ingroup_ASTRALIII_inclusive<-length(unique(annotations$genus[annotations$Library_ID %in% annotations$Library_ID[!grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]]))
number_tribes_ingroup_ASTRALIII_inclusive<-length(sort(unique(annotations$tribe[annotations$Library_ID %in% annotations$Library_ID[!grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]])))-4 #Minus three for doubles in Brassica, Camelinea and Iberideae
number_families_outgroup_ASTRALIII_inclusive<-length(sort(unique(annotations$tribe[annotations$Library_ID %in% annotations$Library_ID[grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]])))
number_samples_new_in_this_study_ASTRALIII_inclusive<-sort(annotations$Library_ID[annotations$Library_ID %in% tree_ASTRALIII_inclusive$tip.label])

##Sample numbers chloroplast----
number_samples_ingroup_chloroplast<-length(annotations$Library_ID[!grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label])
number_samples_outgroup_chloroplast<-length(annotations$Library_ID[grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label])
number_species_ingroup_chloroplast<-length(unique(annotations$species[annotations$Library_ID %in% annotations$Library_ID[!grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label]]))-1 #Minus one for the two populations of Clausia aprica
number_species_outgroup_chloroplast<-length(unique(annotations$species[annotations$Library_ID %in% annotations$Library_ID[grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label]]))
number_genera_ingroup_chloroplast<-length(unique(annotations$genus[annotations$Library_ID %in% annotations$Library_ID[!grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label]]))
number_tribes_ingroup_chloroplast<-length(sort(unique(annotations$tribe[annotations$Library_ID %in% annotations$Library_ID[!grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label]])))-4 #Minus three for doubles in Brassica, Camelinea and Iberideae
number_families_outgroup_chloroplast<-length(sort(unique(annotations$tribe[annotations$Library_ID %in% annotations$Library_ID[grepl("Outgroup_", annotations$tribe) & annotations$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label]])))
number_samples_new_in_this_study_chloroplast<-sort(annotations$Library_ID[annotations$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label])

##Summary metrics for in publication----
(B764_samples_with_target_length<-sum(!is.na(as.numeric(metadata$bp_Nikolov2019[metadata$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]))))
(A353_samples_with_target_length<-sum(!is.na(as.numeric(metadata$bp_Angiosperms353[metadata$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]))))
(B764_mean_target_length<-mean(as.numeric(metadata$bp_Nikolov2019[metadata$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]), na.rm = T))
(B764_mean_target_length_prop<-mean(as.numeric(metadata$bpoftarget_Nikolov2019[metadata$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]), na.rm = T))
(B764_mean_target_length_prop<-mean(as.numeric(metadata$nloci_Nikolov2019[metadata$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]), na.rm = T))

(A353_mean_target_length<-mean(as.numeric(metadata$bp_Angiosperms353[metadata$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]), na.rm = T))
(A353_mean_target_length_prop<-mean(as.numeric(metadata$bpoftarget_Angiosperms353[metadata$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]), na.rm = T))
(A353_mean_target_length_prop<-mean(as.numeric(metadata$nloci_Angiosperms353[metadata$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]), na.rm = T))

#Register which samples were included in which analyses.
annotations$sample_in_inclusive_dataset<-""
annotations$sample_in_inclusive_dataset[annotations$Library_ID %in% tree_ASTRALIII_inclusive$tip.label]<-"Y"
annotations$sample_in_strict_dataset<-""
annotations$sample_in_strict_dataset[annotations$Library_ID %in% tree_ASTRALIII_strict$tip.label]<-"Y"
annotations$sample_in_superstrict_dataset<-""
annotations$sample_in_superstrict_dataset[annotations$Library_ID %in% tree_ASTRALIII_superstrict$tip.label]<-"Y"
annotations$sample_in_superstrict_by_tribe_dataset<-""
annotations$sample_in_superstrict_by_tribe_dataset[annotations$Library_ID %in% tree_ASTRALIII_superstrict_by_tribe$tip.label]<-"Y"
annotations$sample_in_ASTRAL_Pro_dataset<-""
annotations$sample_in_ASTRAL_Pro_dataset[annotations$Library_ID %in% tree_ASTRAL_Pro$tip.label]<-"Y"
annotations$tree_IQ_TREE_supermatrix_concordance_factors<-""
annotations$tree_IQ_TREE_supermatrix_concordance_factors[annotations$Library_ID %in% tree_IQ_TREE_supermatrix_concordance_factors$tip.label]<-"Y"
annotations$sample_in_chloroplast_dataset<-""
annotations$sample_in_chloroplast_dataset[annotations$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label]<-"Y"

#Further list collection details.
annotations$collection<-metadata[,17][match(annotations$Library_ID,metadata$Library_ID)]
annotations$collection_number<-metadata$voucher[match(annotations$Library_ID,metadata$Library_ID)]
annotations$collection_date<-metadata$Collection_year[match(annotations$Library_ID,metadata$Library_ID)]
annotations$collection_date[annotations$collection_date %in% c("NA","unknown")]<-""
annotations$collection_date[is.na(annotations$collection_date)]<-""
annotations$trait_state<-as.character(metadata$Sample_because[match(annotations$Library_ID,metadata$Library_ID)])

#Add taxon occupancy from ASTRAL log-files.
annotations$ASTRAL_taxon_occupancy_inclusive<-as.numeric(as.character(metadata$ASTRAL_taxon_occupancy_inclusive[match(annotations$Library_ID,metadata$Library_ID)]))
annotations$ASTRAL_taxon_occupancy_strict<-as.numeric(as.character(metadata$ASTRAL_taxon_occupancy_strict[match(annotations$Library_ID,metadata$Library_ID)]))
annotations$ASTRAL_taxon_occupancy_superstrict<-as.numeric(as.character(metadata$ASTRAL_taxon_occupancy_superstrict[match(annotations$Library_ID,metadata$Library_ID)]))
annotations$ASTRAL_taxon_occupancy_superstrict_by_tribe<-as.numeric(as.character(metadata$ASTRAL_taxon_occupancy_superstrict_by_tribe[match(annotations$Library_ID,metadata$Library_ID)]))
annotations$ASTRAL_Pro_taxon_occupancy<-as.numeric(as.character(metadata$ASTRAL_PRO_taxon_occupancy[match(annotations$Library_ID,metadata$Library_ID)]))

##ASTRAL sampling overview for publication ----

###KASPERH: need to still update this to include results from ASTRAL-Pro!!!
#Calculate the number of libraries, species, and tribes covered by this phylogenetic tree;
#make sure to exclude outgroup samples when counting species, genera, and tribes!
sampling_overview<-data.frame(sampling_routine=NA,
                              locus_sample_Prop_threshold=NA,
                              locus_sample_length_Prop_threshold=NA,
                              sample_target_length_threshold=NA,
                              sample_locus_Prop_threshold=NA,
                              locus_remove_for_all_samples_SNPs_theshold=NA,
                              sample_remove_outlier_loci=NA,
                              no_libraries=NA,
                              mean_ASTRAL_taxon_occupancy=NA,
                              no_species=NA,
                              no_genera=NA,
                              no_tribes=NA,
                              LPP_mean=NA,
                              LPP_median=NA,
                              q1_mean=NA,
                              q1_median=NA)
#Now fill this dataframe with relevant data.
k=0
for (i in c("inclusive","strict","superstrict","superstrict_by_tribe")){
  k=k+1
  #Add empty row if needed.
  if(sum(is.na(sampling_overview$sampling_routine))==0) {sampling_overview[nrow(sampling_overview)+1,] <- NA}
  #Fill out the data.
  sampling_overview$sampling_routine[k]<-i
  #Read config file for this routine.
  source(paste0("results_HybPhaser_allGenes/config_",i,".txt"))
  #Add the data from the config file to the dataframe.
  sampling_overview$locus_sample_Prop_threshold[k]<-remove_loci_with_less_than_this_propotion_of_samples_recovered
  sampling_overview$locus_sample_length_Prop_threshold[k]<-remove_loci_with_less_than_this_propotion_of_samples_recovered
  sampling_overview$sample_locus_Prop_threshold[k]<-remove_samples_with_less_than_this_propotion_of_loci_recovered
  sampling_overview$sample_target_length_threshold[k]<-remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered
  sampling_overview$locus_remove_for_all_samples_SNPs_theshold[k]<-remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs
  sampling_overview$sample_remove_outlier_loci[k]<-remove_outlier_loci_for_each_sample
  #Collect sampling statistics from ASTRAL-III output.
  sampling_overview$no_libraries[k]<-length(annotations[!grepl("^Outgroup",annotations$tribe) & annotations[,paste0("sample_in_",i,"_dataset")]=="Y",]$Library_ID)
  sampling_overview$mean_ASTRAL_taxon_occupancy[k]<-mean(as.numeric(annotations[,paste0("ASTRAL_taxon_occupancy_",i)]), na.rm = T)
  sampling_overview$no_species[k]<-length(unique(annotations[!grepl("^Outgroup",annotations$tribe) & annotations[,paste0("sample_in_",i,"_dataset")]=="Y",]$species))
  sampling_overview$no_genera[k]<-length(unique(annotations[!grepl("^Outgroup",annotations$tribe) & annotations[,paste0("sample_in_",i,"_dataset")]=="Y",]$genus))
  sampling_overview$no_tribes[k]<-length(unique(annotations[!grepl("^Outgroup",annotations$tribe) & annotations[,paste0("sample_in_",i,"_dataset")]=="Y",]$tribe))
  #Collect node support summary statistics from ASTRAL-III output.
  sampling_overview$LPP_mean[k]<-round(mean(branch_support_compared[which(branch_support_compared$sampling_routine==paste0("ASTRAL-III ",i)),]$pp1, na.rm = T),3)
  sampling_overview$LPP_median[k]<-round(median(branch_support_compared[which(branch_support_compared$sampling_routine==paste0("ASTRAL-III ",i)),]$pp1, na.rm = T),3)
  sampling_overview$q1_mean[k]<-round(mean(branch_support_compared[which(branch_support_compared$sampling_routine==paste0("ASTRAL-III ",i)),]$q1, na.rm = T),3)
  sampling_overview$q1_median[k]<-round(median(branch_support_compared[which(branch_support_compared$sampling_routine==paste0("ASTRAL-III ",i)),]$q1, na.rm = T),3)
}
sampling_overview<-t(sampling_overview)


##Plot comparison of ASTRAL quartet scores ----
plot_quartet_scores<-ggplot(data=branch_support_compared, aes(colour=sampling_routine, x=q1))+
  geom_density()+
  theme_classic()+
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme (plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
ggsave(file = "results_plots_phylogenies/plot_quartet_scores.pdf", plot_quartet_scores + theme(legend.position="none"),
       width = 3, height = 3)
ggsave(file = "results_plots_phylogenies/plot_quartet_scores_legend.pdf", cowplot::get_legend(plot_quartet_scores),
       width = 3, height = 1.5)


##Save table with sampling-overview for publication ----
write.table(sampling_overview, file = "11c.ASTRALIII_sample_summary.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#CONTINUE ANNOTATIONS DATAFRAME ----
#Let's remove the "Outgroup_" in the tribe names to save space when plotting later.
annotations<-data.frame(lapply(annotations, function(x) {gsub("Outgroup_","",x)}))

#Add annotation data for plotting later.
library(stringr)
annotations$type_specimen<-str_split_fixed(annotations$trait_state, ", ",2)[,2]

annotations$genus_type_species<-str_split_fixed(annotations$trait_state, ", ",2)[,1]
annotations$genus_type_species[!(annotations$genus_type_species %in% c("C","W","S","Ta","A","CW","C?","O", "N"))] <- "Y"
annotations$genus_type_species[annotations$genus_type_species %in% c("C","W","S","Ta","A","CW","C?","O", "N")] <- ""

annotations$calibration<-str_split_fixed(annotations$trait_state, ", ",2)[,1]
annotations$calibration[!(annotations$calibration %in% c("T","WT","Ta","S","A","W","O","T?"))] <- "Y"
annotations$calibration[annotations$calibration %in% c("T","WT","Ta","S","A","W","O","T?")] <- ""

annotations$woody_sample<-str_split_fixed(annotations$trait_state, ", ",2)[,1]
annotations$woody_sample[!(annotations$trait_state %in% c("WT","W"))] <- "H"
annotations$woody_sample[annotations$trait_state %in% c("WT","W")] <- "W"

annotations$old_sample<-""
annotations$old_sample[as.integer(as.character(annotations$collection_date))<1900 & !is.na(as.integer(as.character(annotations$collection_date)))]<-"Y"
annotations$old_sample[is.na(annotations$old_sample)]<-""

#Create multiple vectors for "species", dependent on the species being the genus type or not.
#This will assist plotting later.
annotations$species_genus_type<-NA
annotations[annotations$genus_type_species!="",]$species_genus_type<-annotations[annotations$genus_type_species!="",]$species
annotations$species_genus_type[is.na(annotations$species_genus_type)]<-""
annotations$not_species_genus_type<-NA
annotations[annotations$genus_type_species=="",]$not_species_genus_type<-annotations[annotations$genus_type_species=="",]$species
annotations$not_species_genus_type[is.na(annotations$not_species_genus_type)]<-""

#Add results from HybPhaser on AD and LH.
annotations$allele_divergence<-as.numeric(metadata$allele_divergence[match(annotations$Library_ID,metadata$Library_ID)])
annotations$locus_heterozygosity<-as.numeric(metadata$locus_heterozygosity[match(annotations$Library_ID,metadata$Library_ID)])
annotations$nloci<-as.numeric(metadata$nloci[match(annotations$Library_ID,metadata$Library_ID)])
annotations$HybPhaser_bp<-as.numeric(metadata$bp[match(annotations$Library_ID,metadata$Library_ID)])

#Add annotations from Lineage assignments by previous authors.
annotations$Lineages_Koch_Al_Shehbaz<-as.character(metadata$Lineages_Koch_Al_Shehbaz[match(annotations$Library_ID,metadata$Library_ID)])
annotations$Lineages_Franzke_et_al<-metadata$Lineages_Franzke_et_al[match(annotations$Library_ID,metadata$Library_ID)]
annotations$Lineages_Nikolov_et_al<-metadata$Lineages_Nikolov_et_al[match(annotations$Library_ID,metadata$Library_ID)]
annotations$Lineages_unplaced<-metadata$Lineages_unplaced[match(annotations$Library_ID,metadata$Library_ID)]

#Update rownames to be the same as tip labels from trees; this is needed for ggtree's gheatmap function.
rownames(annotations)<-annotations$Library_ID

#Save the annotations table for publication and use in following scripts.
write.table(annotations, file = "11e.Sampling_overview.csv", append = FALSE, sep = ", ", dec = ".",
            row.names = FALSE, col.names = TRUE)


#PLOT SMALL CLADOGRAMS OF TREES FOR COMPARISON ----

plot_tree_ASTRALIII_inclusive_cladogram<-ggtree(drop.tip(tree_ASTRALIII_inclusive, tree_ASTRALIII_inclusive$tip.label[!tree_ASTRALIII_inclusive$tip.label %in% c("S1116", "S1106sl", "S1239", "S1498", "S1093sl", "S1100", "PAFTOL_019361")]),
                                                aes(), branch.length='none', size=1.5) %<+% annotations
plot_tree_ASTRALIII_inclusive_cladogram<-gheatmap(plot_tree_ASTRALIII_inclusive_cladogram, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.6, colnames = FALSE)+
  scale_fill_manual(values=c("transparent","#FF4F33","lightblue","chartreuse3","pink","#FFC425","pink","red","#FFC425"), na.translate = F)+
  geom_tiplab(aes(label=c("I","Ae","Out","III","IV","V","II","","","","","","")), size=6, offset = 1.8, color = "black", hjust='center')+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))

plot_tree_ASTRALIII_strict_cladogram<-ggtree(drop.tip(tree_ASTRALIII_strict, tree_ASTRALIII_strict$tip.label[!tree_ASTRALIII_strict$tip.label %in% c("S1116", "S1106sl", "S1239", "S1498", "S1093sl", "S1100", "PAFTOL_019361")]),
                                             aes(), branch.length='none', size=1.5) %<+% annotations
plot_tree_ASTRALIII_strict_cladogram<-gheatmap(plot_tree_ASTRALIII_strict_cladogram, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.6, colnames = FALSE)+
  scale_fill_manual(values=c("transparent","#FF4F33","lightblue","chartreuse3","pink","#FFC425","pink","red","#FFC425"), na.translate = F)+
  geom_tiplab(aes(label=c("I","IV","Out","Ae","III","V","II","","","","","","")), size=6, offset = 1.8, color = "black", hjust='center')+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))

plot_tree_ASTRALIII_superstrict_cladogram<-ggtree(drop.tip(tree_ASTRALIII_superstrict, tree_ASTRALIII_superstrict$tip.label[!tree_ASTRALIII_superstrict$tip.label %in% c("S1116", "S1106sl", "S1239", "S1498", "S1093sl", "S1100", "PAFTOL_019361")]),
                                                  aes(), branch.length='none', size=1.5) %<+% annotations
plot_tree_ASTRALIII_superstrict_cladogram<-gheatmap(plot_tree_ASTRALIII_superstrict_cladogram, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.6, colnames = FALSE)+
  scale_fill_manual(values=c("transparent","#FF4F33","lightblue","chartreuse3","pink","#FFC425","pink","red","#FFC425"), na.translate = F)+
  geom_tiplab(aes(label=c("I","V","II","Out","Ae","III","IV","","","","","","")), size=6, offset = 1.8, color = "black", hjust='center')+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))

plot_tree_ASTRALIII_superstrict_by_tribe_cladogram<-ggtree(drop.tip(tree_ASTRALIII_superstrict_by_tribe, tree_ASTRALIII_superstrict_by_tribe$tip.label[!tree_ASTRALIII_superstrict_by_tribe$tip.label %in% c("S1116", "S1106sl", "S1239", "S1498", "S1093sl", "S1100", "PAFTOL_019361")]),
                                                           aes(), branch.length='none', size=1.5) %<+% annotations
plot_tree_ASTRALIII_superstrict_by_tribe_cladogram<-gheatmap(plot_tree_ASTRALIII_superstrict_by_tribe_cladogram, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.6, colnames = FALSE)+
  scale_fill_manual(values=c("transparent","#FF4F33","lightblue","chartreuse3","pink","#FFC425","pink","red","#FFC425"), na.translate = F)+
  geom_tiplab(aes(label=c("I","IV","Out","Ae","III","V","II","","","","","","")), size=6, offset = 1.8, color = "black", hjust='center')+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))

plot_tree_ASTRAL_Pro_cladogram<-ggtree(drop.tip(tree_ASTRAL_Pro, tree_ASTRAL_Pro$tip.label[!tree_ASTRAL_Pro$tip.label %in% c("S1116", "S1331", "S1239", "S1498", "S1579", "S1100", "PAFTOL_019361")]),
                                       aes(), branch.length='none', size=1.5) %<+% annotations
plot_tree_ASTRAL_Pro_cladogram<-gheatmap(plot_tree_ASTRAL_Pro_cladogram, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.6, colnames = FALSE)+
  scale_fill_manual(values=c("transparent","#FF4F33","#FFC425","chartreuse3","pink","lightblue","pink","red","#FFC425"), na.translate = F)+
  geom_tiplab(aes(label=c("II","V","III","Out","Ae","I","IV","","","","","","")), size=6, offset = 1.8, color = "black", hjust='center')+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
plot_tree_ASTRAL_Pro_cladogram<-flip(plot_tree_ASTRAL_Pro_cladogram, 2, 1)
#Note in this tree that the last tips needed to be flipped (which does not change the topology, only the view of the topology in a graph;
#this also required to reorder the colours in the heatmap!)

plot_tree_IQ_TREE_supermatrix_cladogram<-ggtree(drop.tip(tree_IQ_TREE_supermatrix_concordance_factors, tree_IQ_TREE_supermatrix_concordance_factors$tip.label[!tree_IQ_TREE_supermatrix_concordance_factors$tip.label %in% c("S1116", "S1106sl", "S1239", "S1498", "S1093sl", "S1100", "PAFTOL_019361")]),
                                                aes(), branch.length='none', size=1.5) %<+% annotations
plot_tree_IQ_TREE_supermatrix_cladogram<-gheatmap(plot_tree_IQ_TREE_supermatrix_cladogram, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.6, colnames = FALSE)+
  scale_fill_manual(values=c("transparent","#FF4F33","lightblue","chartreuse3","pink","#FFC425","pink","red","#FFC425"), na.translate = F)+
  geom_tiplab(aes(label=c("IV","I","V","II","III","Ae","Out","","","","","","")), size=6, offset = 1.8, color = "black", hjust='center')+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))

plot_tree_chloroplast_cladogram<-ggtree(drop.tip(tree_IQ_TREE_chloroplast_calibrated@phylo, tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label[!tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label %in% c("S1116", "S0399", "S1239", "S0748", "S1129", "S1100", "PAFTOL_019361")]),
                                        aes(), branch.length='none', size=1.5) %<+% annotations
plot_tree_chloroplast_cladogram<-gheatmap(plot_tree_chloroplast_cladogram, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.6, colnames = FALSE)+
  scale_fill_manual(values=c("transparent","#FF4F33","lightblue","chartreuse3","pink","#FFC425","pink","red","#FFC425"), na.translate = F)+
  geom_tiplab(aes(label=c("Out","Ae","IV","II*","V","III","I","","","","","","")), size=6, offset = 1.8, color = "black", hjust='center')+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))

#Plot the above cladograms side by side for easy comparison.
plot_trees_combined_cladograms<-plot_grid(plot_tree_IQ_TREE_supermatrix_cladogram,
                                          plot_tree_ASTRAL_Pro_cladogram,
                                          plot_tree_ASTRALIII_inclusive_cladogram,
                                          plot_tree_ASTRALIII_strict_cladogram,
                                          plot_tree_ASTRALIII_superstrict_cladogram,
                                          plot_tree_ASTRALIII_superstrict_by_tribe_cladogram,
                                          plot_tree_chloroplast_cladogram,
                                          nrow = 1,
                                          labels = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)"),
                                          label_size = 22,
                                          label_fontface = "bold")
ggsave(filename = "results_plots_phylogenies/plot_trees_combined_cladograms.png", plot_trees_combined_cladograms,
       width = 12,
       height = 2)
ggsave(filename = "results_plots_phylogenies/plot_trees_combined_cladograms.pdf", plot_trees_combined_cladograms,
       width = 12,
       height = 2)


#PLOT SMALL VERSIONS OF TREES FOR COMPARISON ----

plot_tree_ASTRALIII_inclusive_small<-ggtree(tree_ASTRALIII_inclusive,  aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")
plot_tree_ASTRALIII_inclusive_small<-gheatmap(plot_tree_ASTRALIII_inclusive_small, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.2, colnames = FALSE)+
#                           c("Anast.",         "basal",      "Bisc"     "Cochle" , "I",     "Iberis",    "II",        "III", "IV",     "Mega",        Teesd",     "V")
    scale_fill_manual(values=c("darkolivegreen","transparent","deeppink3","black","#FF4F33","darkblue","lightblue","chartreuse3","pink", "darkcyan", "chocolate", "#FFC425"), na.translate = F)+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_ASTRALIII_inclusive_small.pdf", plot_tree_ASTRALIII_inclusive_small,
       width = 6, height = 6)

plot_tree_ASTRALIII_strict_small<-ggtree(tree_ASTRALIII_strict,  aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")
plot_tree_ASTRALIII_strict_small<-gheatmap(plot_tree_ASTRALIII_strict_small, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.2, colnames = FALSE)+
  scale_fill_manual(values=c("darkolivegreen","transparent","deeppink3","black","#FF4F33","darkblue","lightblue","chartreuse3","pink", "darkcyan", "chocolate", "#FFC425"), na.translate = F)+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_ASTRALIII_strict_small.pdf", plot_tree_ASTRALIII_strict_small,
       width = 6, height = 6)

plot_tree_ASTRALIII_superstrict_small<-ggtree(tree_ASTRALIII_superstrict,  aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")
plot_tree_ASTRALIII_superstrict_small<-gheatmap(plot_tree_ASTRALIII_superstrict_small, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.2, colnames = FALSE)+
  scale_fill_manual(values=c("darkolivegreen","transparent","deeppink3","black","#FF4F33","darkblue","lightblue","chartreuse3","pink", "darkcyan", "chocolate", "#FFC425"), na.translate = F)+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_ASTRALIII_superstrict_small.pdf", plot_tree_ASTRALIII_superstrict_small,
       width = 6, height = 6)

plot_tree_ASTRALIII_superstrict_by_tribe_small<-ggtree(tree_ASTRALIII_superstrict_by_tribe,  aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")
plot_tree_ASTRALIII_superstrict_by_tribe_small<-gheatmap(plot_tree_ASTRALIII_superstrict_by_tribe_small, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.2, colnames = FALSE)+
  scale_fill_manual(values=c("darkolivegreen","transparent","deeppink3","black","#FF4F33","darkblue","lightblue","chartreuse3","pink", "darkcyan", "chocolate", "#FFC425"), na.translate = F)+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_ASTRALIII_superstrict_by_tribe_small.pdf", plot_tree_ASTRALIII_superstrict_by_tribe_small,
       width = 6, height = 6)

#For ASTRAL-Pro output, for which are listed the full node annotations, check how the LPP ("pp1") value was selected for annotation in this plot!
plot_tree_ASTRAL_Pro_small<-ggtree(tree_ASTRAL_Pro,  aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")
plot_tree_ASTRAL_Pro_small<-gheatmap(plot_tree_ASTRAL_Pro_small, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.2, colnames = FALSE)+
  scale_fill_manual(values=c("darkolivegreen","transparent","deeppink3","black","#FF4F33","darkblue","lightblue","chartreuse3","pink", "darkcyan", "chocolate", "#FFC425"), na.translate = F)+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_ASTRAL_Pro_small.pdf", plot_tree_ASTRAL_Pro_small,
       width = 6, height = 6)

plot_tree_IQ_TREE_supermatrix_small<-ggtree(tree_IQ_TREE_supermatrix_concordance_factors,  aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")
plot_tree_IQ_TREE_supermatrix_small<-gheatmap(plot_tree_IQ_TREE_supermatrix_small, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.2, colnames = FALSE)+
  scale_fill_manual(values=c("darkolivegreen","transparent","deeppink3","black","#FF4F33","darkblue","lightblue","chartreuse3","pink", "darkcyan", "chocolate", "#FFC425"), na.translate = F)+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_IQ_TREE_supermatrix_small.pdf", plot_tree_IQ_TREE_supermatrix_small,
       width = 6, height = 6)

#Take care to update the colour scale below for the main Lineages after final tree data comes in!
#Now Megacarpaeeae missing and thus one colour removed by hand...!
plot_tree_chloroplast_small<-ggtree(tree_IQ_TREE_chloroplast_calibrated@phylo,  aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")
plot_tree_chloroplast_small<-gheatmap(plot_tree_chloroplast_small, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.2, colnames = FALSE)+
  scale_fill_manual(values=c("darkolivegreen","transparent","deeppink3","black","#FF4F33","darkblue","lightblue","chartreuse3","pink", "darkcyan", "chocolate", "#FFC425"), na.translate = F)+
  theme (legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_chloroplast_small.pdf", plot_tree_chloroplast_small,
       width = 6, height = 6)

#Plot one of the legends for inclusion in publication plot later.
ggsave(filename = "results_plots_phylogenies/plot_tree_legend_small.pdf", cowplot::get_legend(plot_tree_ASTRALIII_inclusive_small + theme(legend.position= "right")),
       width = 3, height = 3)



#PLOT FULL-SIZE VERSIONS OF TREES FOR PUBLICATION ----

####TEMP NEW PLOT CLOCKLIKE GENES
plot_tree_ASTRALIII_superstrict_and_clocklike_large<-ggtree(tree_ASTRALIII_superstrict_and_clocklike,  aes(color="grey"), size=0.3, branch.length = "none") %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")+
  geom_inset(tree_ASTRALIII_superstrict_and_clocklike_astral_annotations_piecharts, width = 0.05, height = 0.05, hjust = 0.175, vjust = 0.15)+
  #Highlight the calibration nodes.
  geom_tiplab(aes(label=Library_ID2), size=1, offset = 8.1, color = "black")+
  geom_tiplab(aes(label=tribe), size=1, offset = 2.7, color = "black")+
  geom_tiplab(aes(label=species), fontface='italic', size=1, offset = 4.5, color = "black")
plot_tree_ASTRALIII_superstrict_and_clocklike_large<-plot_tree_ASTRALIII_superstrict_and_clocklike_large+new_scale_fill()
plot_tree_ASTRALIII_superstrict_and_clocklike_large<-gheatmap(plot_tree_ASTRALIII_superstrict_and_clocklike_large, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.05, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("Nikolov","jumpy"), colnames_offset_y = 2)+
  # scale_fill_manual(values=c("darkslategrey","transparent","darkmagenta","blue","#FF4F33","lightblue","chartreuse3","pink","red","#FFC425"), na.translate = F)
  scale_fill_manual(values=c("darkslategrey","transparent","darkmagenta","#FF4F33","lightblue","chartreuse3","pink","red","#FFC425","red","#FFC425"), na.translate = F)
# plot_tree_ASTRALIII_superstrict_and_clocklike_large<-plot_tree_ASTRALIII_superstrict_and_clocklike_large+new_scale_fill()
# plot_tree_ASTRALIII_superstrict_and_clocklike_large<-gheatmap(plot_tree_ASTRALIII_superstrict_and_clocklike_large, data = annotations[,c("locus_heterozygosity"), drop=FALSE], width=.01, offset = 1.8, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("LH"), colnames_offset_y = 2)+
#   scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
# plot_tree_ASTRALIII_superstrict_and_clocklike_large<-plot_tree_ASTRALIII_superstrict_and_clocklike_large+new_scale_fill()
# plot_tree_ASTRALIII_superstrict_and_clocklike_large<-gheatmap(plot_tree_ASTRALIII_superstrict_and_clocklike_large, data = annotations[,c("allele_divergence"), drop=FALSE], width=.01, offset = 2.1, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("AD"), colnames_offset_y = 2) +
#   scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
# theme(legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_ASTRALIII_superstrict_and_clocklike_large.pdf", plot_tree_ASTRALIII_superstrict_and_clocklike_large,
       width = 12, height = 18)



##ML nuclear ----

#Find all tribes_n in the tree for proper annotation below.
tribes_n<-unique(metadata[metadata$Library_ID %in% tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label, ]$Tribe)
tribes_n<-tribes_n[!grepl("Outgroup_", tribes_n)]

plot_tree_IQ_TREE_supermatrix_large<-ggtree(tree_IQ_TREE_supermatrix_calibrated,  aes(color="black"), size=0.25) %<+% annotations +
  scale_color_manual(values = c("black","#00BFC4"), na.value = "black")+
  #Highlight the calibration nodes.
  geom_point2(aes(subset=(node %in% c(MRCA(tree_IQ_TREE_supermatrix_calibrated, "PAFTOL_014331", "S1408")))), shape=23, size=0.2, fill='red', stroke=NA)+
  # geom_point2(aes(subset=(node %in% c(MRCA(tree_IQ_TREE_supermatrix_calibrated, "[ml]", "HYZL"), MRCA(tree_IQ_TREE_supermatrix_calibrated, "PAFTOL_014331", "S1408"), MRCA(tree_IQ_TREE_supermatrix_calibrated, "PAFTOL_019361", "S1408"), MRCA(tree_IQ_TREE_supermatrix_calibrated, "S0416", "S1376")))), shape=23, size=0.5, fill='red', stroke=NA)+
  geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=0.5)
plot_tree_IQ_TREE_supermatrix_large<-plot_tree_IQ_TREE_supermatrix_large+new_scale_color()
plot_tree_IQ_TREE_supermatrix_large<-plot_tree_IQ_TREE_supermatrix_large +
  geom_nodelab(aes(x=branch, label=concordance_factors, color=support), vjust=1.4, size=1)+
  scale_colour_manual(values=sort(unique(tree_IQ_TREE_supermatrix_calibrated@data$support)[c(2,3,4,5)]))+
  geom_nodelab(aes(x=branch, label=bootstrap), vjust=-.6, size=1, color = "black")
#Loop through the tribes_n to add annotation layers to the ggtree object.
for (i in 1:length(tribes_n)){
  tribes_n[i]
  samples_in_tribe<-length(which(tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label %in% annotations[annotations$tribe==tribes_n[i],]$Library_ID2))
  node_tribe<-getMRCA(tree_IQ_TREE_supermatrix_calibrated@phylo, which(tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label %in% annotations[annotations$tribe==tribes_n[i],]$Library_ID2))
  plot_tree_IQ_TREE_supermatrix_large<-plot_tree_IQ_TREE_supermatrix_large+
    geom_cladelab(node=
                    if(samples_in_tribe==1)
                    {which(tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label %in% annotations[annotations$tribe==tribes_n[i],]$Library_ID2)}
                  else
                  {node_tribe}, 
                  label=tribes_n[i], align=TRUE, offset = 11, barsize=2, fontsize=1.5)
}
plot_tree_IQ_TREE_supermatrix_large<-plot_tree_IQ_TREE_supermatrix_large+
  # geom_tiplab(aes(label=tribe), size=1, offset = 14, color = "black")+
  # geom_tiplab(aes(label=Library_ID2), size=1, offset = 7, color = "black") +
  geom_tiplab(aes(label=species_genus_type), fontface='bold.italic', size=1, offset = 0.1, color = "black") +
  geom_tiplab(aes(label=not_species_genus_type), fontface='italic', size=1, offset = 0.1, color = "black") +
  theme_tree2()
plot_tree_IQ_TREE_supermatrix_large<-plot_tree_IQ_TREE_supermatrix_large+new_scale_fill()
plot_tree_IQ_TREE_supermatrix_large<-gheatmap(plot_tree_IQ_TREE_supermatrix_large, data = annotations[,c("Lineages_Nikolov_et_al"), drop=F], width=.01, offset=16, colnames=F)+
  scale_fill_manual(values=c("transparent", "orange", "lightblue","green","pink","yellow"), na.translate = F)
# plot_tree_IQ_TREE_supermatrix_large<-plot_tree_IQ_TREE_supermatrix_large+new_scale_fill()
# plot_tree_IQ_TREE_supermatrix_large<-gheatmap(plot_tree_IQ_TREE_supermatrix_large, data = annotations[,c("locus_heterozygosity"), drop=FALSE], width=.01, offset = 5.2, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("LH"), colnames_offset_y = 2)+
#   scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
# plot_tree_IQ_TREE_supermatrix_large<-plot_tree_IQ_TREE_supermatrix_large+new_scale_fill()
# plot_tree_IQ_TREE_supermatrix_large<-gheatmap(plot_tree_IQ_TREE_supermatrix_large, data = annotations[,c("allele_divergence"), drop=FALSE], width=.01, offset = 6, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("AD"), colnames_offset_y = 2) +
#   scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
#   theme_tree2(legend.position= "right", plot.margin = unit(c(0,0,0,0), "cm")) 
# plot_tree_IQ_TREE_supermatrix_large<-revts(plot_tree_IQ_TREE_supermatrix_large)

ggsave(filename = "results_plots_phylogenies/plot_tree_IQ_TREE_supermatrix_large.pdf", plot_tree_IQ_TREE_supermatrix_large + theme(legend.position="none"),
       width = 8, height = 21)


##ML nuclear circular ----

#Remove outgroups from the tree that we don't want to print now. This helps reduce branch lengths of most ancestral branches.
tree_IQ_TREE_supermatrix_calibrated_ingroups<-drop.tip(tree_IQ_TREE_supermatrix_calibrated, which(tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label %in% annotations$Library_ID2[!annotations$genus %in% c(list_of_Brassicaceae_ingroup_genera, 'Cleome')]))

plot_tree_IQ_TREE_supermatrix_large_circular<-ggtree(tree_IQ_TREE_supermatrix_calibrated_ingroups,  aes(color="black"), size=0.25, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("black","#00BFC4"), na.value = "black")+
  geom_nodelab(aes(x=branch, label=round(height_median,1)), vjust=-.6, size=1, color = "darkgreen")+
  geom_nodelab(aes(x=branch, label=sCF), vjust=-.6, size=1, color = "cyan3")+
  geom_nodelab(aes(x=branch, label=gCF), vjust=-.6, size=1, color = "darksalmon")
plot_tree_IQ_TREE_supermatrix_large_circular<-plot_tree_IQ_TREE_supermatrix_large_circular+
  new_scale_color()+
  geom_point2(aes(fill=sCF_support), shape=21 , stroke=0.5, colour="deeppink2", size=1.5)+
  geom_point2(aes(fill=gCF_support), shape=21 , stroke=0.5, colour="firebrick2", size=1)+
  scale_fill_manual(values = c("black","darkgrey","grey", "lightgrey"), na.value = "darkseagreen2")
#Loop through the tribes_n to add annotation layers to the ggtree object.
for (i in 1:length(tribes_n)){
  tribes_n[i]
  samples_in_tribe<-length(which(tree_IQ_TREE_supermatrix_calibrated_ingroups@phylo$tip.label %in% annotations[annotations$tribe==tribes_n[i],]$Library_ID2))
  if(samples_in_tribe>0){
    node_tribe<-getMRCA(tree_IQ_TREE_supermatrix_calibrated_ingroups@phylo, which(tree_IQ_TREE_supermatrix_calibrated_ingroups@phylo$tip.label %in% annotations[annotations$tribe==tribes_n[i],]$Library_ID2))
    plot_tree_IQ_TREE_supermatrix_large_circular<-plot_tree_IQ_TREE_supermatrix_large_circular+
      geom_cladelab(node=
                      if(samples_in_tribe==1)
                      {which(tree_IQ_TREE_supermatrix_calibrated_ingroups@phylo$tip.label %in% annotations[annotations$tribe==tribes_n[i],]$Library_ID2)}
                    else
                    {node_tribe}, 
                    label=tribes_n[i], align=TRUE, offset = 11, barsize=2, fontsize=1.5)
  }
}
plot_tree_IQ_TREE_supermatrix_large_circular<-plot_tree_IQ_TREE_supermatrix_large_circular+
  # geom_tiplab(aes(label=tribe), size=1, offset = 14, color = "black")+
  # geom_tiplab(aes(label=Library_ID2), size=1, offset = 7, color = "black") +
  geom_tiplab(aes(label=species_genus_type), fontface='bold.italic', size=1, offset = 0.3, color = "black") +
  geom_tiplab(aes(label=not_species_genus_type), fontface='italic', size=1, offset = 0.3, color = "black") +
  theme_tree2()
plot_tree_IQ_TREE_supermatrix_large_circular<-plot_tree_IQ_TREE_supermatrix_large_circular+new_scale_fill()
plot_tree_IQ_TREE_supermatrix_large_circular<-gheatmap(plot_tree_IQ_TREE_supermatrix_large_circular, data = annotations[,c("Lineages_Nikolov_et_al"), drop=F], width=.03, offset=14, colnames=F)+
  scale_fill_manual(values=c("transparent", "orange", "lightblue","green","pink","yellow"), na.translate = F)
# plot_tree_IQ_TREE_supermatrix_large_circular<-plot_tree_IQ_TREE_supermatrix_large_circular+new_scale_fill()
# plot_tree_IQ_TREE_supermatrix_large_circular<-gheatmap(plot_tree_IQ_TREE_supermatrix_large_circular, data = annotations[,c("locus_heterozygosity"), drop=FALSE], width=.01, offset = 5.2, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("LH"), colnames_offset_y = 2)+
#   scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
# plot_tree_IQ_TREE_supermatrix_large_circular<-plot_tree_IQ_TREE_supermatrix_large_circular+new_scale_fill()
# plot_tree_IQ_TREE_supermatrix_large_circular<-gheatmap(plot_tree_IQ_TREE_supermatrix_large_circular, data = annotations[,c("allele_divergence"), drop=FALSE], width=.01, offset = 6, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("AD"), colnames_offset_y = 2) +
#   scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
#   theme_tree2(legend.position= "right", plot.margin = unit(c(0,0,0,0), "cm")) 
# plot_tree_IQ_TREE_supermatrix_large_circular<-revts(plot_tree_IQ_TREE_supermatrix_large_circular)

ggsave(filename = "results_plots_phylogenies/plot_tree_IQ_TREE_supermatrix_large_circular_for_print.pdf", plot_tree_IQ_TREE_supermatrix_large_circular + theme(legend.position="none"),
       width = 10, height = 10)


##ML nuclear excluding hybrids ----

#Create a plot of the phylogeny with branch lengths.
plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large<-ggtree(tree_IQ_TREE_supermatrix_excluding_hybrids_concordance_factors,  aes(color="black"), size=0.25) %<+% annotations +
  scale_color_manual(values = c("black","#00BFC4"), na.value = "black")+
  geom_nodelab(aes(x=branch, label=c(rep("", length(tree_IQ_TREE_supermatrix_excluding_hybrids_concordance_factors$tip.label)),tree_IQ_TREE_supermatrix_excluding_hybrids_concordance_factors$node.label)), vjust=-.3, size=1.2, color = "black")

#And a second plot with the annotations.
#We do this because ggtree cannot properly add annotations to a tree that has branch lengths, so we need to be creative here.
plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large_annotations<-ggtree(tree_IQ_TREE_supermatrix_excluding_hybrids_concordance_factors,  aes(color="black"), size=0.25, branch.length = "none") %<+% annotations +
  scale_color_manual(values = c("grey","grey"), na.value = "black") +
  geom_tiplab(aes(label=tribe), size=1, offset = 1, color = "black") +
  geom_tiplab(aes(label=Library_ID2),  size=1, offset = 6, color = "black")
plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large_annotations<-gheatmap(plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large_annotations, data = annotations[,c("Lineages_Nikolov_et_al"), drop=F], offset=0.5, colnames=F)+
  scale_fill_manual(values=c("transparent", "orange", "lightblue","green","pink","yellow"), na.translate = F)
plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large_annotations<-plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large_annotations +
  geom_tiplab(aes(label=species_genus_type), fontface='bold.italic', size=1, offset = 3, color = "black") +
  geom_tiplab(aes(label=not_species_genus_type), fontface='italic', size=1, offset = 3, color = "black")
  
#Combine the above two plots.
plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large_combined<-
  cowplot::plot_grid(plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large,
                     plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large_annotations,
                     nrow = 1)

#And save for publication.
ggsave(filename = "results_plots_phylogenies/plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large.pdf", plot_tree_IQ_TREE_supermatrix_excluding_hybrids_large_combined + theme(legend.position="none"),
       width = 18, height = 10)


  


##ML plastome ----

#Copy the annotations object in order for us to make small adjustments relating to polyphyly in the plastome phylogeny.
annotations_cp<-annotations
annotations_cp$tribe[annotations_cp$Library_ID=="SRR8528386"]<-"Brassiceae"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637787"]<-"Sisymbrieae II"
annotations_cp$tribe[annotations_cp$Library_ID=="S0753"]<-"Cremolobeae III"
annotations_cp$tribe[annotations_cp$Library_ID=="S1095"]<-"Schizopetaleae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637821"]<-"Schizopetaleae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637768"]<-"Cremolobeae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637749"]<-"Cremolobeae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637760"]<-"Cremolobeae IIIb"
annotations_cp$tribe[annotations_cp$Library_ID=="S0399"]<-"Eudemeae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637822"]<-"Conringieae II"
annotations_cp$tribe[annotations_cp$Library_ID=="S1371"]<-"Hemilophieae trib. nov. II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637702"]<-"Hemilophieae trib. nov. II"

#Add outgroup family names as in metadata dataframe.
annotations_cp$tribe2<-metadata$Tribe[match(annotations_cp$Library_ID, metadata$Library_ID)]

#Find all tribes in the tree for proper annotation below.
tribes_cp<-unique(annotations_cp[annotations_cp$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label & !grepl("Outgroup_", annotations_cp$tribe2), ]$tribe)

plot_tree_IQ_TREE_chloroplast_large<-ggtree(tree_IQ_TREE_chloroplast_calibrated,  aes(color="black"), size=0.25) %<+% annotations_cp +
  scale_color_manual(values = c("black","#00BFC4"), na.value = "black")+
  #Highlight the calibration nodes.
  # geom_point2(aes(subset=(node %in% c(MRCA(tree_IQ_TREE_chloroplast_calibrated, "PAFTOL_014331", "S1408")))), shape=23, size=0.2, fill='red', stroke=NA)+
  geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=0.5)
plot_tree_IQ_TREE_chloroplast_large<-plot_tree_IQ_TREE_chloroplast_large+new_scale_color()
plot_tree_IQ_TREE_chloroplast_large<-plot_tree_IQ_TREE_chloroplast_large +
  geom_nodelab(aes(x=branch, label=concordance_factors, color=support), vjust=1.4, size=1)+
  scale_colour_manual(values=sort(unique(tree_IQ_TREE_chloroplast_calibrated@data$support)[c(2,3,5,4)]))+
  geom_nodelab(aes(x=branch, label=bootstrap), vjust=-.6, size=1, color = "black")
#Loop through the tribes_cp to add annotation layers to the ggtree object.
for (i in 1:length(tribes_cp)){
  tribes_cp[i]
  samples_in_tribe<-length(which(tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label %in% annotations_cp[annotations_cp$tribe==tribes_cp[i],]$Library_ID2))
  node_tribe<-getMRCA(tree_IQ_TREE_chloroplast_calibrated@phylo, which(tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label %in% annotations_cp[annotations_cp$tribe==tribes_cp[i],]$Library_ID2))
  plot_tree_IQ_TREE_chloroplast_large<-plot_tree_IQ_TREE_chloroplast_large+
    geom_cladelab(node=
                    if(samples_in_tribe==1)
                    {which(tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label %in% annotations_cp[annotations_cp$tribe==tribes_cp[i],]$Library_ID2)}
                  else
                  {node_tribe}, 
                  label=tribes_cp[i], align=TRUE, offset = 14, barsize=2, fontsize=1.5)
}
plot_tree_IQ_TREE_chloroplast_large<-plot_tree_IQ_TREE_chloroplast_large+
  geom_tiplab(aes(label=tribe), size=1, offset = 11, color = "black")+
  geom_tiplab(aes(label=Library_ID2), size=1, offset = 5.5, color = "black") +
  geom_tiplab(aes(label=collection_number), size=1, offset = 8, color = "black") +
  geom_tiplab(aes(label=species_genus_type), fontface='bold.italic', size=1, offset = 0.1, color = "black") +
  geom_tiplab(aes(label=not_species_genus_type), fontface='italic', size=1, offset = 0.1, color = "black") +
  theme_tree2()
plot_tree_IQ_TREE_chloroplast_large<-plot_tree_IQ_TREE_chloroplast_large+new_scale_fill()
plot_tree_IQ_TREE_chloroplast_large<-gheatmap(plot_tree_IQ_TREE_chloroplast_large, data = annotations_cp[,c("Lineages_Nikolov_et_al"), drop=F], width=.01, offset=13, colnames=F)+
  scale_fill_manual(values=c("transparent", "orange", "lightblue","green","pink","yellow"), na.translate = F)
# plot_tree_IQ_TREE_chloroplast_large<-plot_tree_IQ_TREE_chloroplast_large+new_scale_fill()
# plot_tree_IQ_TREE_chloroplast_large<-gheatmap(plot_tree_IQ_TREE_chloroplast_large, data = annotations_cp[,c("locus_heterozygosity"), drop=FALSE], width=.01, offset = 5.2, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("LH"), colnames_offset_y = 2)+
#   scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
# plot_tree_IQ_TREE_chloroplast_large<-plot_tree_IQ_TREE_chloroplast_large+new_scale_fill()
# plot_tree_IQ_TREE_chloroplast_large<-gheatmap(plot_tree_IQ_TREE_chloroplast_large, data = annotations_cp[,c("allele_divergence"), drop=FALSE], width=.01, offset = 6, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("AD"), colnames_offset_y = 2) +
#   scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
#   theme_tree2(legend.position= "right", plot.margin = unit(c(0,0,0,0), "cm")) 
# plot_tree_IQ_TREE_chloroplast_large<-revts(plot_tree_IQ_TREE_chloroplast_large)

ggsave(filename = "results_plots_phylogenies/plot_tree_IQ_TREE_chloroplast_large.pdf", plot_tree_IQ_TREE_chloroplast_large + theme(legend.position="right"),
       width = 18, height = 25)


##ML plastome circular ----

#Copy the annotations object in order for us to make small adjustments relating to polyphyly in the plastome phylogeny.
annotations_cp<-annotations
annotations_cp$tribe[annotations_cp$Library_ID=="SRR8528386"]<-"Brassiceae"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637787"]<-"Sisymbrieae II"
annotations_cp$tribe[annotations_cp$Library_ID=="S0753"]<-"Cremolobeae III"
annotations_cp$tribe[annotations_cp$Library_ID=="S1095"]<-"Schizopetaleae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637821"]<-"Schizopetaleae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637768"]<-"Cremolobeae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637749"]<-"Cremolobeae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637760"]<-"Cremolobeae IIIb"
annotations_cp$tribe[annotations_cp$Library_ID=="S0399"]<-"Eudemeae II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637822"]<-"Conringieae II"
annotations_cp$tribe[annotations_cp$Library_ID=="S1371"]<-"Hemilophieae trib. nov. II"
annotations_cp$tribe[annotations_cp$Library_ID=="MK637702"]<-"Hemilophieae trib. nov. II"

#Add outgroup family names as in metadata dataframe.
annotations_cp$tribe2<-metadata$Tribe[match(annotations_cp$Library_ID, metadata$Library_ID)]

#Remove outgroups from the tree that we don't want to print now. This helps reduce branch lengths of most ancestral branches.
tree_IQ_TREE_chloroplast_calibrated_ingroups<-drop.tip(tree_IQ_TREE_chloroplast_calibrated, which(tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label %in% annotations_cp$Library_ID2[!annotations_cp$genus %in% c(list_of_Brassicaceae_ingroup_genera, 'Cleome')]))

#Find all tribes in the tree for proper annotation below.
tribes_cp<-unique(annotations_cp[annotations_cp$Library_ID %in% tree_IQ_TREE_chloroplast_calibrated_ingroups@phylo$tip.label & !grepl("Outgroup_", annotations_cp$tribe2), ]$tribe)

plot_tree_IQ_TREE_chloroplast_large_circular<-ggtree(tree_IQ_TREE_chloroplast_calibrated_ingroups, aes(color="black"), size=0.25, layout='circular') %<+% annotations_cp +
  scale_color_manual(values = c("black","#00BFC4"), na.value = "black")+
  geom_nodelab(aes(x=branch, label=round(height_median,1)), vjust=-.6, size=1, color = "darkgreen")+
  geom_nodelab(aes(x=branch, label=sCF), vjust=-.6, size=1, color = "cyan3")+
  geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=0.5)
plot_tree_IQ_TREE_chloroplast_large_circular<-plot_tree_IQ_TREE_chloroplast_large_circular+
  new_scale_color()+
  geom_point2(aes(fill=sCF_support), shape=21, colour = "deeppink", size=1.5)+
  scale_fill_manual(values = c("black","darkgrey","grey", "lightgrey"), na.value = "darkseagreen2")

plot_tree_IQ_TREE_chloroplast_large_circular<-plot_tree_IQ_TREE_chloroplast_large_circular+
  new_scale_color()
#Loop through the tribes_cp to add annotation layers to the ggtree object.
for (i in 1:length(tribes_cp)){
  tribes_cp[i]
  samples_in_tribe<-length(which(tree_IQ_TREE_chloroplast_calibrated_ingroups@phylo$tip.label %in% annotations_cp[annotations_cp$tribe==tribes_cp[i],]$Library_ID2))
  if(samples_in_tribe>1){
    node_tribe<-getMRCA(tree_IQ_TREE_chloroplast_calibrated_ingroups@phylo, which(tree_IQ_TREE_chloroplast_calibrated_ingroups@phylo$tip.label %in% annotations_cp[annotations_cp$tribe==tribes_cp[i],]$Library_ID2))
    plot_tree_IQ_TREE_chloroplast_large_circular<-plot_tree_IQ_TREE_chloroplast_large_circular+
      geom_cladelab(node=node_tribe, label=tribes_cp[i], align=TRUE, offset = 14, barsize=2, fontsize=1.5,textcolor=col_vector[i], barcolor=col_vector[i] )
  }
}
plot_tree_IQ_TREE_chloroplast_large_circular<-plot_tree_IQ_TREE_chloroplast_large_circular+
  geom_tiplab(aes(label=species_genus_type), fontface='bold.italic', size=1, offset = 0.4, color = "black") +
  geom_tiplab(aes(label=not_species_genus_type), fontface='italic', size=1, offset = 0.4, color = "black") +
  theme_tree2()

plot_tree_IQ_TREE_chloroplast_large_circular<-plot_tree_IQ_TREE_chloroplast_large_circular+
  new_scale_fill()
plot_tree_IQ_TREE_chloroplast_large_circular<-gheatmap(plot_tree_IQ_TREE_chloroplast_large_circular, data = annotations_cp[,c("Lineages_Nikolov_et_al"), drop=F], width=.01, offset=13, colnames=F)+
  scale_fill_manual(values=c("transparent", "orange", "lightblue","green","pink","yellow"), na.translate = F)

ggsave(filename = "results_plots_phylogenies/plot_tree_IQ_TREE_chloroplast_large_circular_for_print.pdf", plot_tree_IQ_TREE_chloroplast_large_circular + theme(legend.position="none"),
       width = 10, height = 10)


##ASTRAL-III inclusive ----

plot_tree_ASTRALIII_inclusive_large<-ggtree(tree_ASTRALIII_inclusive,  aes(color="grey"), size=0.3, branch.length = "none") %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")+
  geom_inset(tree_ASTRALIII_inclusive_astral_annotations_piecharts, width = 0.05, height = 0.05, hjust = 0.175, vjust = 0.15)+
  #Highlight the calibration nodes.
  geom_tiplab(aes(label=ASTRAL_taxon_occupancy_inclusive), size=0.6, offset = 7.6, color = "black")+
  geom_tiplab(aes(label=Library_ID2), size=1, offset = 8.1, color = "black")+
  geom_tiplab(aes(label=tribe), size=1, offset = 2.7, color = "black")+
  geom_tiplab(aes(label=species), fontface='italic', size=1, offset = 4.5, color = "black")
plot_tree_ASTRALIII_inclusive_large<-plot_tree_ASTRALIII_inclusive_large+new_scale_fill()
plot_tree_ASTRALIII_inclusive_large<-gheatmap(plot_tree_ASTRALIII_inclusive_large, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.05, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("Nikolov","Franzke","Koch","jumpy"), colnames_offset_y = 2)+
  scale_fill_manual(values=c("darkslategrey","transparent","darkmagenta","blue","#FF4F33","lightblue","chartreuse3","pink","red","#FFC425"), na.translate = F)
plot_tree_ASTRALIII_inclusive_large<-plot_tree_ASTRALIII_inclusive_large+new_scale_fill()
plot_tree_ASTRALIII_inclusive_large<-gheatmap(plot_tree_ASTRALIII_inclusive_large, data = annotations[,c("locus_heterozygosity"), drop=FALSE], width=.01, offset = 1.8, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("LH"), colnames_offset_y = 2)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
plot_tree_ASTRALIII_inclusive_large<-plot_tree_ASTRALIII_inclusive_large+new_scale_fill()
plot_tree_ASTRALIII_inclusive_large<-gheatmap(plot_tree_ASTRALIII_inclusive_large, data = annotations[,c("allele_divergence"), drop=FALSE], width=.01, offset = 2.1, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("AD"), colnames_offset_y = 2) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
  theme(legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_ASTRALIII_inclusive_large.pdf", plot_tree_ASTRALIII_inclusive_large,
       width = 12, height = 18)

##ASTRAL-III strict ----

plot_tree_ASTRALIII_strict_large<-ggtree(tree_ASTRALIII_strict,  aes(color="grey"), size=0.3, branch.length = "none") %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")+
  geom_inset(tree_ASTRALIII_strict_astral_annotations_piecharts, width = 0.05, height = 0.05, hjust = 0.175, vjust = 0.15)+
  #Highlight the calibration nodes.
  geom_tiplab(aes(label=ASTRAL_taxon_occupancy_strict), size=0.6, offset = 7.6, color = "black")+
  geom_tiplab(aes(label=Library_ID2), size=1, offset = 8.1, color = "black")+
  geom_tiplab(aes(label=tribe), size=1, offset = 2.7, color = "black")+
  geom_tiplab(aes(label=species), fontface='italic', size=1, offset = 4.5, color = "black")
plot_tree_ASTRALIII_strict_large<-plot_tree_ASTRALIII_strict_large+new_scale_fill()
plot_tree_ASTRALIII_strict_large<-gheatmap(plot_tree_ASTRALIII_strict_large, data = annotations[,c("Lineages_Nikolov_et_al","Lineages_unplaced")], width=.05, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("Nikolov","Franzke","Koch","jumpy"), colnames_offset_y = 2)+
  scale_fill_manual(values=c("darkslategrey","transparent","darkmagenta","blue","#FF4F33","lightblue","chartreuse3","pink","red","#FFC425"), na.translate = F)
plot_tree_ASTRALIII_strict_large<-plot_tree_ASTRALIII_strict_large+new_scale_fill()
plot_tree_ASTRALIII_strict_large<-gheatmap(plot_tree_ASTRALIII_strict_large, data = annotations[,c("locus_heterozygosity"), drop=FALSE], width=.01, offset = 1.8, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("LH"), colnames_offset_y = 2)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
plot_tree_ASTRALIII_strict_large<-plot_tree_ASTRALIII_strict_large+new_scale_fill()
plot_tree_ASTRALIII_strict_large<-gheatmap(plot_tree_ASTRALIII_strict_large, data = annotations[,c("allele_divergence"), drop=FALSE], width=.01, offset = 2.1, colnames_position = "top", colnames_angle = 45, font.size = 2, custom_column_labels = c("AD"), colnames_offset_y = 2) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
  theme(legend.position= "none", plot.margin = unit(c(0,0,0,0), "cm"))
ggsave(filename = "results_plots_phylogenies/plot_tree_ASTRALIII_strict_large.pdf", plot_tree_ASTRALIII_strict_large,
       width = 12, height = 18)


#ASTRAL ROUTINE BINARY HEATMAPS ----

##Inclusive routine ----
#Read in the names of all gene tree files.
tree_ASTRALIII_inclusive_gene_trees_files<-list.files("results_iqtree_final_inclusive", pattern="*.treefile")
#Sort by origin of bait kit (using simply length of the file name).
tree_ASTRALIII_inclusive_gene_trees_files<-tree_ASTRALIII_inclusive_gene_trees_files[order(nchar(tree_ASTRALIII_inclusive_gene_trees_files))]
#Create a dataframe to be filled with binary data for presence/absence for each sample in each gene tree.
tree_ASTRALIII_inclusive_presence<-data.frame(matrix(nrow=length(tree_ASTRALIII_inclusive$tip.label), ncol=length(tree_ASTRALIII_inclusive_gene_trees_files)), 0)
colnames(tree_ASTRALIII_inclusive_presence)<-lapply(strsplit(tree_ASTRALIII_inclusive_gene_trees_files, "_"), "[", 1)
rownames(tree_ASTRALIII_inclusive_presence)<-tree_ASTRALIII_inclusive$tip.label
#Loop through the gene trees and tip labels to score presence of the samples in each gene tree.
for (i in 1:length(tree_ASTRALIII_inclusive_gene_trees_files)){
  gene_tree_i<-read.tree(paste0("results_iqtree_final_inclusive/", tree_ASTRALIII_inclusive_gene_trees_files[i]))
  for (j in 1:length(gene_tree_i$tip.label)){
    gene_tree_i_tip_label_j<-gene_tree_i$tip.label[j]
    #And put a "1" in the table in the correct cell.
    tree_ASTRALIII_inclusive_presence[which(rownames(tree_ASTRALIII_inclusive_presence) == gene_tree_i_tip_label_j), i] <- 1
  }
}
#Reformat the dataframe into long format.
tree_ASTRALIII_inclusive_presence<-tree_ASTRALIII_inclusive_presence[,which(!is.na(colnames(tree_ASTRALIII_inclusive_presence))==T)] #Remove any columns with a colname that does not exist (NA)
tree_ASTRALIII_inclusive_presence$samples<-rownames(tree_ASTRALIII_inclusive_presence)
tree_ASTRALIII_inclusive_presence_long<-tree_ASTRALIII_inclusive_presence %>% gather(genes, value, colnames(tree_ASTRALIII_inclusive_presence)[1:length(colnames(tree_ASTRALIII_inclusive_presence))-1])


##Strict routine ----
#Read in the names of all gene tree files.
tree_ASTRALIII_strict_gene_trees_files<-list.files("results_iqtree_final_strict", pattern="*.treefile")
#Sort by origin of bait kit (using simply length of the file name).
tree_ASTRALIII_strict_gene_trees_files<-tree_ASTRALIII_strict_gene_trees_files[order(nchar(tree_ASTRALIII_strict_gene_trees_files))]
#Create a dataframe to be filled with binary data for presence/absence for each sample in each gene tree.
tree_ASTRALIII_strict_presence<-data.frame(matrix(nrow=length(tree_ASTRALIII_strict$tip.label), ncol=length(tree_ASTRALIII_strict_gene_trees_files)), 0)
colnames(tree_ASTRALIII_strict_presence)<-lapply(strsplit(tree_ASTRALIII_strict_gene_trees_files, "_"), "[", 1)
rownames(tree_ASTRALIII_strict_presence)<-tree_ASTRALIII_strict$tip.label
#Loop through the gene trees and tip labels to score presence of the samples in each gene tree.
for (i in 1:length(tree_ASTRALIII_strict_gene_trees_files)){
  gene_tree_i<-read.tree(paste0("results_iqtree_final_strict/", tree_ASTRALIII_strict_gene_trees_files[i]))
  for (j in 1:length(gene_tree_i$tip.label)){
    gene_tree_i_tip_label_j<-gene_tree_i$tip.label[j]
    #And put a "1" in the table in the correct cell.
    tree_ASTRALIII_strict_presence[which(rownames(tree_ASTRALIII_strict_presence) == gene_tree_i_tip_label_j), i] <- 1
  }
}
#Reformat the dataframe into long format.
tree_ASTRALIII_strict_presence<-tree_ASTRALIII_strict_presence[,which(!is.na(colnames(tree_ASTRALIII_strict_presence))==T)] #Remove any columns with a colname that does not exist (NA)
tree_ASTRALIII_strict_presence$samples<-rownames(tree_ASTRALIII_strict_presence)
tree_ASTRALIII_strict_presence_long<-tree_ASTRALIII_strict_presence %>% gather(genes, value, colnames(tree_ASTRALIII_strict_presence)[1:length(colnames(tree_ASTRALIII_strict_presence))-1])


##Superstrict routine ----
#Read in the names of all gene tree files.
tree_ASTRALIII_superstrict_gene_trees_files<-list.files("results_iqtree_final_superstrict", pattern="*.treefile")
#Sort by origin of bait kit (using simply length of the file name).
tree_ASTRALIII_superstrict_gene_trees_files<-tree_ASTRALIII_superstrict_gene_trees_files[order(nchar(tree_ASTRALIII_superstrict_gene_trees_files))]
#Create a dataframe to be filled with binary data for presence/absence for each sample in each gene tree.
tree_ASTRALIII_superstrict_presence<-data.frame(matrix(nrow=length(tree_ASTRALIII_superstrict$tip.label), ncol=length(tree_ASTRALIII_superstrict_gene_trees_files)), 0)
colnames(tree_ASTRALIII_superstrict_presence)<-lapply(strsplit(tree_ASTRALIII_superstrict_gene_trees_files, "_"), "[", 1)
rownames(tree_ASTRALIII_superstrict_presence)<-tree_ASTRALIII_superstrict$tip.label
#Loop through the gene trees and tip labels to score presence of the samples in each gene tree.
for (i in 1:length(tree_ASTRALIII_superstrict_gene_trees_files)){
  gene_tree_i<-read.tree(paste0("results_iqtree_final_superstrict/", tree_ASTRALIII_superstrict_gene_trees_files[i]))
  for (j in 1:length(gene_tree_i$tip.label)){
    gene_tree_i_tip_label_j<-gene_tree_i$tip.label[j]
    #And put a "1" in the table in the correct cell.
    tree_ASTRALIII_superstrict_presence[which(rownames(tree_ASTRALIII_superstrict_presence) == gene_tree_i_tip_label_j), i] <- 1
  }
}
#Reformat the dataframe into long format.
tree_ASTRALIII_superstrict_presence<-tree_ASTRALIII_superstrict_presence[,which(!is.na(colnames(tree_ASTRALIII_superstrict_presence))==T)] #Remove any columns with a colname that does not exist (NA)
tree_ASTRALIII_superstrict_presence$samples<-rownames(tree_ASTRALIII_superstrict_presence)
tree_ASTRALIII_superstrict_presence_long<-tree_ASTRALIII_superstrict_presence %>% gather(genes, value, colnames(tree_ASTRALIII_superstrict_presence)[1:length(colnames(tree_ASTRALIII_superstrict_presence))-1])


##Superstrict by tribe routine ----
#Read in the names of all gene tree files.
tree_ASTRALIII_superstrict_by_tribe_gene_trees_files<-list.files("results_iqtree_final_superstrict_by_tribe", pattern="*.treefile")
#Sort by origin of bait kit (using simply length of the file name).
tree_ASTRALIII_superstrict_by_tribe_gene_trees_files<-tree_ASTRALIII_superstrict_by_tribe_gene_trees_files[order(nchar(tree_ASTRALIII_superstrict_by_tribe_gene_trees_files))]
#Create a dataframe to be filled with binary data for presence/absence for each sample in each gene tree.
tree_ASTRALIII_superstrict_by_tribe_presence<-data.frame(matrix(nrow=length(tree_ASTRALIII_superstrict_by_tribe$tip.label), ncol=length(tree_ASTRALIII_superstrict_by_tribe_gene_trees_files)), 0)
colnames(tree_ASTRALIII_superstrict_by_tribe_presence)<-lapply(strsplit(tree_ASTRALIII_superstrict_by_tribe_gene_trees_files, "_"), "[", 1)
rownames(tree_ASTRALIII_superstrict_by_tribe_presence)<-tree_ASTRALIII_superstrict_by_tribe$tip.label
#Loop through the gene trees and tip labels to score presence of the samples in each gene tree.
for (i in 1:length(tree_ASTRALIII_superstrict_by_tribe_gene_trees_files)){
  gene_tree_i<-read.tree(paste0("results_iqtree_final_superstrict_by_tribe/", tree_ASTRALIII_superstrict_by_tribe_gene_trees_files[i]))
  for (j in 1:length(gene_tree_i$tip.label)){
    gene_tree_i_tip_label_j<-gene_tree_i$tip.label[j]
    #And put a "1" in the table in the correct cell.
    tree_ASTRALIII_superstrict_by_tribe_presence[which(rownames(tree_ASTRALIII_superstrict_by_tribe_presence) == gene_tree_i_tip_label_j), i] <- 1
  }
}
#Reformat the dataframe into long format.
tree_ASTRALIII_superstrict_by_tribe_presence<-tree_ASTRALIII_superstrict_by_tribe_presence[,which(!is.na(colnames(tree_ASTRALIII_superstrict_by_tribe_presence))==T)] #Remove any columns with a colname that does not exist (NA)
tree_ASTRALIII_superstrict_by_tribe_presence$samples<-rownames(tree_ASTRALIII_superstrict_by_tribe_presence)
tree_ASTRALIII_superstrict_by_tribe_presence_long<-tree_ASTRALIII_superstrict_by_tribe_presence %>% gather(genes, value, colnames(tree_ASTRALIII_superstrict_by_tribe_presence)[1:length(colnames(tree_ASTRALIII_superstrict_by_tribe_presence))-1])

##Pro routine ----
#Note that many samples occur multiple times in a gene tree; these are the different alleles saved by the HybPiper paralog script.
#So specifically for the ASTRAL_Pro output we'll score number of occurrence instead of presence/absence.
#Read in the names of all gene tree files.
tree_ASTRAL_Pro_gene_trees_files<-list.files("results_ASTRAL_Pro/results_iqtree_gene_trees/", pattern="*.treefile")
#Sort by origin of bait kit (using simply length of the file name).
tree_ASTRAL_Pro_gene_trees_files<-tree_ASTRAL_Pro_gene_trees_files[order(nchar(tree_ASTRAL_Pro_gene_trees_files))]
#Create a dataframe to be filled with binary data for presence/absence for each sample in each gene tree.
tree_ASTRAL_Pro_presence<-data.frame(matrix(nrow=length(tree_ASTRAL_Pro$tip.label), ncol=length(tree_ASTRAL_Pro_gene_trees_files)), 0)
colnames(tree_ASTRAL_Pro_presence)<-lapply(strsplit(tree_ASTRAL_Pro_gene_trees_files, "_"), "[", 2)
rownames(tree_ASTRAL_Pro_presence)<-tree_ASTRAL_Pro$tip.label
#Loop through the gene trees and tip labels to score presence of the samples in each gene tree.
for (i in 1:length(tree_ASTRAL_Pro_gene_trees_files)){
  gene_tree_i<-read.tree(paste0("results_ASTRAL_Pro/results_iqtree_gene_trees/", tree_ASTRAL_Pro_gene_trees_files[i]))
  #Specifically for the ASTRAL-Pro output, change remove ".*" from the tip labels, i.e. the indication of different allele versions for a sample in the tree.
  gene_tree_i$tip.label<-sapply(strsplit(gene_tree_i$tip.label, "[.]"), "[", 1)
  gene_tree_i_table<-table(gene_tree_i$tip.label)
  for (j in 1:length(gene_tree_i_table)){
    #And log the count data from the gene tree label table in the overall heatmap table.
    tree_ASTRAL_Pro_presence[which(rownames(tree_ASTRAL_Pro_presence) == names(gene_tree_i_table)[j]), i] <- as.numeric(gene_tree_i_table[j])
  }
}
#We need to remove one gene for which an unlikely high number of alleles was found for some samples, A353 gene 6128.
tree_ASTRAL_Pro_presence<-tree_ASTRAL_Pro_presence[,which(!colnames(tree_ASTRAL_Pro_presence)=="6128")]
#Reformat the dataframe into long format.
tree_ASTRAL_Pro_presence$samples<-rownames(tree_ASTRAL_Pro_presence)
tree_ASTRAL_Pro_presence_long<-tree_ASTRAL_Pro_presence %>% gather(genes, value, colnames(tree_ASTRAL_Pro_presence)[1:length(colnames(tree_ASTRAL_Pro_presence))-1])


##Combine & plot results ----
#Combine results.
tree_ASTRAL_all_routines_presence_long<-rbind(cbind(tree_ASTRALIII_inclusive_presence_long, routine="ASTRAL-III inclusive"),
                                              cbind(tree_ASTRALIII_strict_presence_long, routine="ASTRAL-III strict"),
                                              cbind(tree_ASTRALIII_superstrict_presence_long, routine="ASTRAL-III superstrict"),
                                              cbind(tree_ASTRALIII_superstrict_by_tribe_presence_long, routine="ASTRAL-III superstrict by tribe"),
                                              cbind(tree_ASTRAL_Pro_presence_long, routine="ASTRAL-Pro"))
#Change all NAs to zeros (but remember that a zero can mean "no data available", or "data removed by HybPhaser").
tree_ASTRAL_all_routines_presence_long$value[is.na(tree_ASTRAL_all_routines_presence_long$value)]<-0

#Add some (meta)data for ordering the heatmap plot.
tree_ASTRAL_all_routines_presence_long$presence<-as.integer(tree_ASTRAL_all_routines_presence_long$value) #This prevents ggplot to introduce a continuous scale.
#Add tribe and main lineage information so we can order samples roughly by their systematic position.
tree_ASTRAL_all_routines_presence_long$tribe<-metadata$Tribe[match(tree_ASTRAL_all_routines_presence_long$samples, metadata$Library_ID)]
tree_ASTRAL_all_routines_presence_long$Lineages_Nikolov_et_al<-metadata$Lineages_Nikolov_et_al[match(tree_ASTRAL_all_routines_presence_long$samples, metadata$Library_ID)]
tree_ASTRAL_all_routines_presence_long$Lineages_Nikolov_et_al[which(tree_ASTRAL_all_routines_presence_long$Lineages_Nikolov_et_al=="")]<-"unplaced"
tree_ASTRAL_all_routines_presence_long$Lineages_Nikolov_et_al[which(tree_ASTRAL_all_routines_presence_long$Lineages_Nikolov_et_al=="basal")]<-"Aethionemeae"
tree_ASTRAL_all_routines_presence_long$Lineages_Nikolov_et_al[which(grepl("Outgroup_", tree_ASTRAL_all_routines_presence_long$tribe))]<-"Outgroup"
tree_ASTRAL_all_routines_presence_long$sample_label<-paste0(tree_ASTRAL_all_routines_presence_long$Lineages_Nikolov_et_al, " - ", tree_ASTRAL_all_routines_presence_long$tribe, " - ", tree_ASTRAL_all_routines_presence_long$samples)
#Add a distinction for the origin of the genes based on bait kit.
tree_ASTRAL_all_routines_presence_long$bait_kit<-NA
tree_ASTRAL_all_routines_presence_long$bait_kit[nchar(tree_ASTRAL_all_routines_presence_long$genes)<5]<-"A353"
tree_ASTRAL_all_routines_presence_long$bait_kit[nchar(tree_ASTRAL_all_routines_presence_long$genes)>5]<-"B764"


#Add order by (1) main lineage cf. Nikolov et al (2019) (y-axis) and (2) tribe (also y-axis).
#First arrange the dataframe in the desired order.
tree_ASTRAL_all_routines_presence_long<-tree_ASTRAL_all_routines_presence_long %>%
  arrange(factor(tribe, levels = sort(unique(tree_ASTRAL_all_routines_presence_long$tribe)))) %>%
  arrange(factor(Lineages_Nikolov_et_al, levels = c("Outgroup", "Aethionemeae", "I", "II", "III", "IV", "V", "unplaced")))
#Add a column to give a subsequent number that can be used to order in the plot.
tree_ASTRAL_all_routines_presence_long$y_axis_plot_order<-1:nrow(tree_ASTRAL_all_routines_presence_long)

tree_ASTRAL_gene_heatmaps<-ggplot(tree_ASTRAL_all_routines_presence_long, 
                                  aes(x=fct_reorder(genes, presence, .fun='sum'), 
                                      y=fct_reorder(sample_label, desc(y_axis_plot_order))))+
  geom_tile(aes(fill=as.character(presence)))+
  facet_grid(rows = vars(routine))+
  scale_fill_manual(values=c("transparent", colorRamps::blue2red(9)))+
  theme(axis.text.y = element_text(size=1), 
        axis.text.x = element_text(angle = 90, hjust=1, size=1),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  xlab("Genes ordered by sum of presence across all samples and routines")+
  ylab("Samples ordered by main lineage (sensu Nikolov, et al. 2019) > tribe")+
  guides(fill=guide_legend("Presence or number of alleles"))

ggsave(file = "results_plots_phylogenies/tree_ASTRAL_gene_heatmaps12.pdf", tree_ASTRAL_gene_heatmaps,
       width = 12, height = 25)

ggsave(file = "results_plots_phylogenies/tree_ASTRAL_gene_heatmaps12.png", tree_ASTRAL_gene_heatmaps,
       width = 12, height = 25, dpi = 600)





#quick test for qhy empty colums
library(dplyr)
totals<-tree_ASTRAL_all_routines_presence_long %>% group_by(genes) %>% summarise(sum = sum(value)) %>%as.data.frame




###TEMP READ IN PAFTOL TREE AND PRUNE ALL EXCEPT BRASSICALES
tree_PAFTOL<-read.tree(file="10.PAFTOL_treeoflife.all_support_values.2.0.tree")
#Collect tip labels and create a vector of tips to drop.
tip_labels_PAFTOL_all<-tree_PAFTOL$tip.label
tip_labels_PAFTOL_droptips<-tip_labels_PAFTOL_all[!grepl("Brassicaceae|Cleome_albescens", tip_labels_PAFTOL_all)]
#Prune the tree to only save Brassicales
tree_PAFTOL_Brassicaceae<-drop.tip(tree_PAFTOL, tip_labels_PAFTOL_droptips)
#Root the tree.
tree_PAFTOL_Brassicaceae<-root(tree_PAFTOL_Brassicaceae, outgroup = "Brassicales_Cleomaceae_Cleome_albescens_ERR7621937")
#Create a dataframe with annotations.
annotations_PAFTOL<-data.frame(sample=tree_PAFTOL_Brassicaceae$tip.label)
library(tidyr)
annotations_PAFTOL<-cbind(annotations_PAFTOL,
                          separate(annotations_PAFTOL, sample, into = c("Order", "Family", "Genus", "Species", "Sample"), "_"))
#Update some synonyms to the currently accepted names.
annotations_PAFTOL$Genus[annotations_PAFTOL$Genus=="Hugueninia"]<-"Descurainia"
annotations_PAFTOL$Genus[annotations_PAFTOL$Genus=="Cuphonotus"]<-"Cuphonotus (Lemphoria)"
annotations_PAFTOL$Genus[annotations_PAFTOL$Genus=="Zuvanda"]<-"Plagioloba"
annotations_PAFTOL$Genus[annotations_PAFTOL$Genus=="Ochthodium"]<-"Sisymbrium"
#Define the final annotation label for plotting.
annotations_PAFTOL$tiplabel<-paste0(annotations_PAFTOL$Genus, " ", annotations_PAFTOL$Species)
#Add lineage designation, for now taken from the large dataframe annotations for the phylogenies Produced from own data.
annotations_PAFTOL$Lineages_Koch_Al_Shehbaz<-as.character(metadata$Lineages_Koch_Al_Shehbaz[match(annotations_PAFTOL$Genus,metadata$Genus)])
annotations_PAFTOL$Lineages_Franzke_et_al<-as.character(metadata$Lineages_Franzke_et_al[match(annotations_PAFTOL$Genus,metadata$Genus)])
annotations_PAFTOL$Lineages_Nikolov_et_al<-as.character(metadata$Lineages_Nikolov_et_al[match(annotations_PAFTOL$Genus,metadata$Genus)])
#Similarly add tribes for Brassicaceae species.
annotations_PAFTOL$Tribes<-as.character(annotations$tribe[match(annotations_PAFTOL$Genus,annotations$genus)])
annotations_PAFTOL$Tribes[which(grepl("Outgroup_", annotations_PAFTOL$Tribes))]<-""
#Remove all NA from table.
#annotations_PAFTOL<-na.omit(annotations_PAFTOL)
#Create a ggtree for plotting.
#Plot the 1000nuclearGenes species tree.
PAFTOL_species_tree_v2_circular<-ggtree(tree_PAFTOL_Brassicaceae,  aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations_PAFTOL +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")+
  ggnewscale::new_scale_fill()+
  geom_tippoint(aes(x=x+0.5,fill=Lineages_Nikolov_et_al), shape=22, size=2, stroke=0, colour="transparent")+
  scale_fill_manual(values = c("white","white","orange","lightblue","lightgreen","pink","yellow"), na.translate = F)+
  ggnewscale::new_scale_fill()+
  geom_tippoint(aes(x=x+1.2, fill=Lineages_Franzke_et_al), shape=22, size=2, stroke=0, colour="transparent")+
  scale_fill_manual(values = c("white","white","blue","orange","lightblue","lightgreen"), na.translate = F)+
  ggnewscale::new_scale_fill()+
  geom_tippoint(aes(x=x+1.9, fill=Lineages_Koch_Al_Shehbaz), shape=22, size=2, stroke=0, colour="transparent")+
  scale_fill_manual(values = c("white","white","orange","lightblue","lightgreen"), na.translate = F)+
  geom_tiplab(aes(label=tiplabel), size=0.6, offset = 8, color = "black")+
  geom_tiplab(aes(label=Tribes), size=0.6, offset = 2.5, color = "black")+
  #geom_text(aes(label=c(rep("",length(tree_PAFTOL_Brassicaceae$tip.label)),tree_PAFTOL_Brassicaceae$node.label)), size=0.6, color = "black", hjust=2, vjust=-0.7)+
  theme(legend.position = "none")
ggsave(filename="10.PAFTOL_species_tree_v2_circular.pdf", PAFTOL_species_tree_v2_circular,
       width = 9, height = 9)




##TEMP READ IN ASTRAL-Pro MICROLIPEDIEAE TREE
tree_astral_Pro_microlepidieae<-read.tree(file = "6i.ASTRAL_Pro_species_tree.tre")
tree_astral_Pro_microlepidieae<-root(tree_astral_Pro_microlepidieae, outgroup = "S0435")
tip_labels<-tree_astral_Pro_microlepidieae$tip.label
annotations<-data.frame(Library_ID=tip_labels)
#Add taxonomic data to annotations dataframe.
annotations$Library_ID2<-metadata$Library_ID[match(annotations$Library_ID,metadata$Library_ID)]
annotations$tribe<-metadata$Tribe[match(annotations$Library_ID,metadata$Library_ID)]
annotations$genus<-metadata$Genus[match(annotations$Library_ID,metadata$Library_ID)]
list_of_Brassicaceae_ingroup_genera<-unique(annotations$genus[!grepl("Outgroup_", annotations$tribe)])
annotations$species<-metadata$Name[match(annotations$Library_ID,metadata$Library_ID)]
ggtree(tree_astral_Pro_microlepidieae,  aes(color="grey"), branch.length='none', size=0.3) %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")+
  geom_tiplab(aes(label=Library_ID2), size=2, offset = 11.2, color = "black")+
  #geom_tiplab(aes(label=collection), size=0.6, offset = 13.2, color = "black")+
  #geom_tiplab(aes(label=ASTRAL_taxon_occupancy_strict), size=0.6, offset = 9.7, color = "black")+
  #geom_tiplab(aes(label="#", color = as.numeric(locus_heterozygosity)), fontface='bold', size=3, offset = 0.5)+
  #geom_tiplab(aes(label=collection_number), size=0.6, offset = 14.7, color = "black")+
  geom_tiplab(aes(label=tribe), size=2, offset = 0.5, color = "black")+
  geom_tiplab(aes(label=species), fontface='italic', size=2, offset = 4, color = "black")+
  geom_text(aes(label=c(rep("",length(tree_astral_Pro_microlepidieae$tip.label)),tree_astral_Pro_microlepidieae$node.label)), size=3, color = "black")
####END TEMP

###TEMP Read in the results from the iq-tree concatenated supermatrix analysis ran by Lars Nauheimer
tree_supermatrix_iqtree_superstrict<-read.tree(file = "../2022-02-17_Brassicaceae_nuclear_concatenated_iqtree_by_Lars_Nauheimer/concord_superstrict_tribe_1.cf.tree")
#Root the tree to the outgroup
tree_supermatrix_iqtree_superstrict<-root(tree_supermatrix_iqtree_superstrict, outgroup = "PAFTOL_019361")
#Create a ggtree for plotting.
#Plot the 1000nuclearGenes species tree.
IQ_TREE_superstrict_circular<-ggtree(tree_supermatrix_iqtree_superstrict,  aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")+
  ggnewscale::new_scale_fill()+
  geom_tippoint(aes(x=x+0.5,fill=Lineages_Nikolov_et_al), shape=22, size=2, stroke=0, colour="transparent")+
  scale_fill_manual(values = c("white","white","orange","lightblue","lightgreen","pink","yellow"))+
  ggnewscale::new_scale_fill()+
  geom_tippoint(aes(x=x+1.2, fill=Lineages_Franzke_et_al), shape=22, size=2, stroke=0, colour="transparent")+
  scale_fill_manual(values = c("white","white","blue","orange","lightblue","lightgreen"))+
  ggnewscale::new_scale_fill()+
  geom_tippoint(aes(x=x+1.9, fill=Lineages_Koch_Al_Shehbaz), shape=22, size=2, stroke=0, colour="transparent")+
  scale_fill_manual(values = c("white","white","orange","lightblue","lightgreen"))+
  geom_tiplab(aes(label=Library_ID2), size=0.6, offset = 11.2, color = "black")+
  geom_tiplab(aes(label=tribe), size=0.6, offset = 2.5, color = "black")+
  geom_tiplab(aes(label=species), fontface='italic', size=0.6, offset = 6, color = "black")+
  geom_text(aes(label=c(rep("",length(tree_supermatrix_iqtree_superstrict$tip.label)),tree_supermatrix_iqtree_superstrict$node.label)), size=0.6, color = "black", hjust=0, vjust=0)+
  theme(legend.position = "none")
ggsave(filename="9.IQ_TREE_superstrict_circular.pdf", IQ_TREE_superstrict_circular,
       width = 9, height = 9)




####TEMP SCRIPT BELOW -- CAN BE DELETED LATER
L1<-c("S1460", "S1235sl", "S0316sl", "S0809", "S0282", "S1227", "S1369", "S0690", "S0293", "S0400", "S0268", "S0301", "S_670", "S0763", "S1204", "S0670", "S0718", "S0641", "S1318-2", "S0721sl", "S1381", "S1211", "SRR8528391", "S0904", "S0671", "S0734sl", "S0733sl", "S0298", "S1276", "S0593", "S1541sl", "S0427", "S0672", "S0893", "S0376", "SRR8528343", "S0756", "S1407", "S0558", "S0412sl", "S0742", "S0762-2", "S0055", "S1279", "S1370", "S0405", "S0643", "S1278", "S1181-2", "S0751", "S1561", "S0591", "S0547", "S1375", "S0509", "ERR2789774", "S0306", "S0285", "S1116", "S0728sl", "SRR8528342", "S0896", "S0947sl", "S0941", "S0741", "S1400", "S0435", "S0949", "S0895", "S0946", "S1424", "SRR8528350", "S0612", "S0294", "S0903", "S0300", "S0764", "S1052sl", "S0383", "S0275", "S0689", "S0722", "S0717", "SRR8528357", "SRR8528359", "S1180")
not_L1<-tree_ASTRALIII_inclusive$tip.label[!(tree_ASTRALIII_inclusive$tip.label %in% L1)]


#TEMP plot Lineage 1 for Martin Lysak.
tree_ASTRALIII_inclusive_L1<-drop.tip(tree_ASTRALIII_inclusive ,not_L1)
tree_ASTRALIII_inclusive_L1<-root.phylo(tree_ASTRALIII_inclusive_L1, "S1180")
L1_plot_MartinLysak<-ggtree(tree_ASTRALIII_inclusive_L1, aes(color="grey"), branch.length='none', size=0.3, layout = 'circular') %<+% annotations +
  scale_color_manual(values = c("darkgrey","#00BFC4"), na.value = "darkgrey")+
  #geom_tiplab(aes(label=tribe), size=2, offset = 3.7, color = "black")+
  geom_tiplab(aes(label=species), fontface='italic', size=2, offset = 0.4, color = "black")+
  geom_tiplab(aes(label=species), fontface='italic', size=2, offset = 10, color = "transparent")+
  geom_text(aes(label=c(rep("",length(tree_ASTRALIII_inclusive_L1$tip.label)),tree_ASTRALIII_inclusive_L1$node.label)), size=2, color = "black", hjust=1.3, vjust=-0.5)+
  theme(legend.position= "none")
ggsave(filename = "L1_plot_MartinLysak.pdf", L1_plot_MartinLysak, width = 8, height = 8)


#STUDY TRIBAL SUPPORT AND AGE----

#Setup dataframe to store results in.

#Use the vector of Brassicaceae ingroup tribes created above for nuclear dataset.
# list_of_Brassicaceae_ingroup_tribes_n<-sort(list_of_Brassicaceae_ingroup_tribes_n[!is.na(list_of_Brassicaceae_ingroup_tribes_n)])
list_of_Brassicaceae_ingroup_tribes_n<-tribes_n
   
#List any tribes to exclude.
list_of_Brassicaceae_ingroup_tribes_n_to_exclude<-c("Brassiceae II", "Camelineae IIIb")

#Use the vector of Brassicaceae ingroup tribes created above for plastome dataset.
# list_of_Brassicaceae_ingroup_tribes_cp<-sort(list_of_Brassicaceae_ingroup_tribes_cp[!is.na(list_of_Brassicaceae_ingroup_tribes_cp)])
list_of_Brassicaceae_ingroup_tribes_cp<-tribes_cp

#List any tribes to exclude.
list_of_Brassicaceae_ingroup_tribes_cp_to_exclude<-c("Brassiceae II", "Schizopetaleae II", "Camelineae IIIb", "Cremolobeae IIIb", "Cremolobeae III", "Cremolobeae II", "Sisymbrieae II", "Eudemeae II", "Conringieae II", "Hemilophieae trib. nov. II")

#Create combined vector.
tribes_combined<-sort(unique(c(list_of_Brassicaceae_ingroup_tribes_n[!list_of_Brassicaceae_ingroup_tribes_n %in% list_of_Brassicaceae_ingroup_tribes_n_to_exclude],
                   list_of_Brassicaceae_ingroup_tribes_cp[!list_of_Brassicaceae_ingroup_tribes_cp %in% list_of_Brassicaceae_ingroup_tribes_cp_to_exclude])))

#Create a dataframe to store results.
# support_tribal<-data.frame(clade=list_of_Brassicaceae_ingroup_tribes_n[!list_of_Brassicaceae_ingroup_tribes_n %in% list_of_Brassicaceae_ingroup_tribes_n_to_exclude],
support_tribal<-data.frame(clade=tribes_combined,
                             main_lineage=NA,
                             nuclear_BS=NA,
                             nuclear_gCF=NA,
                             nuclear_sCF=NA,
                             nuclear_LPP=NA,
                             nuclear_Q1=NA,
                             plastome_BS=NA,
                             plastome_sCF=NA,
                             nuclear_crown_age_median=NA,
                             nuclear_crown_age_95CI_lower=NA,
                             nuclear_crown_age_95CI_upper=NA,
                             nuclear_stem_age_median=NA,
                             nuclear_stem_age_95CI_lower=NA,
                             nuclear_stem_age_95CI_upper=NA,
                             plastome_crown_age_median=NA,
                             plastome_crown_age_95CI_lower=NA,
                             plastome_crown_age_95CI_upper=NA,
                             plastome_stem_age_median=NA,
                             plastome_stem_age_95CI_lower=NA,
                             plastome_stem_age_95CI_upper=NA)

#Loop through the tribes in the phylogenies to find and store values of interest.
#Note that for quartet scores from the coalescent ASTRAL analysis we use the 'strict' routine.
for (t in 1:length(tribes_combined)){
  tribe<-tribes_combined[t]
  #List main lineage following Nikolov et al. (2019).
  support_tribal$main_lineage[t]<-annotations$Lineages_Nikolov_et_al[which(annotations$tribe==tribe)[1]]
  
  #Note we first need to check if the tribe is in the phylogeny at all.
  if (length(which(annotations$tribe[match(tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label, annotations$Library_ID)]==tribe))>0){
    #To find nodal support from the ML supermatrix approach, first find the node with the MRCA.
    #Note that in cases where a tribe is represented by a single sample, we cannot know these values, except for the stem age.
    #We then need to take the parent node of the mrca to get the stem age.
    tribal_node_ML<-MRCA(tree_IQ_TREE_supermatrix_calibrated, which(annotations$tribe[match(tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label, annotations$Library_ID)]==tribe))
    if (length(which(annotations$tribe[match(tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label, annotations$Library_ID)]==tribe))==1){
      support_tribal$nuclear_BS[t]<-"-"
      support_tribal$nuclear_gCF[t]<-"-"
      support_tribal$nuclear_sCF[t]<-"-"
      support_tribal$nuclear_crown_age_median[t]<-"-"
      support_tribal$nuclear_crown_age_95CI_lower[t]<-"-"
      support_tribal$nuclear_crown_age_95CI_upper[t]<-"-"
      support_tribal$nuclear_stem_age_median[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, tribal_node_ML),]$height_median
      support_tribal$nuclear_stem_age_95CI_lower[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, tribal_node_ML),]$height_0.95_HPD[[1]][1]
      support_tribal$nuclear_stem_age_95CI_upper[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, tribal_node_ML),]$height_0.95_HPD[[1]][2]
      
    } else {
      #Find and save nodal support.
      support_tribal$nuclear_BS[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==tribal_node_ML,]$bootstrap
      support_tribal$nuclear_gCF[t]<-strsplit(tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==tribal_node_ML,]$concordance_factors, "/", 1)[[1]][1]
      support_tribal$nuclear_sCF[t]<-strsplit(tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==tribal_node_ML,]$concordance_factors, "/", 1)[[1]][2]
      #Find and save time calibration data.
      support_tribal$nuclear_crown_age_median[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==tribal_node_ML,]$height_median
      support_tribal$nuclear_crown_age_95CI_lower[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==tribal_node_ML,]$height_0.95_HPD[[1]][1]
      support_tribal$nuclear_crown_age_95CI_upper[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==tribal_node_ML,]$height_0.95_HPD[[1]][2]
      #Repeat this for the stem age, so focusing on the parent of the MRCA.
      support_tribal$nuclear_stem_age_median[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, tribal_node_ML),]$height_median
      support_tribal$nuclear_stem_age_95CI_lower[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, tribal_node_ML),]$height_0.95_HPD[[1]][1]
      support_tribal$nuclear_stem_age_95CI_upper[t]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, tribal_node_ML),]$height_0.95_HPD[[1]][2]
    }
  } else {
    support_tribal$nuclear_BS[t]<-"-"
    support_tribal$nuclear_gCF[t]<-"-"
    support_tribal$nuclear_sCF[t]<-"-"
    support_tribal$nuclear_crown_age_median[t]<-"-"
    support_tribal$nuclear_crown_age_95CI_lower[t]<-"-"
    support_tribal$nuclear_crown_age_95CI_upper[t]<-"-"
    support_tribal$nuclear_stem_age_median[t]<-"-"
    support_tribal$nuclear_stem_age_95CI_lower[t]<-"-"
    support_tribal$nuclear_stem_age_95CI_upper[t]<-"-"
  }
  
  #Repeat for the coalescent analysis.
  #Note we first need to check if the tribe is in the phylogeny at all.
  if (length(which(annotations$tribe[match(tree_ASTRALIII_strict_astral_annotations@phylo$tip.label, annotations$Library_ID)]==tribe))>0){
    tribal_node_coal<-MRCA(tree_ASTRALIII_strict_astral_annotations@phylo, which(annotations$tribe[match(tree_ASTRALIII_strict_astral_annotations@phylo$tip.label, annotations$Library_ID)]==tribe))
    #Find and save nodal support.
    #Note that here, in the case that a tribe is represented by a single sample only, we cannot know the tribal support.
    if(length(which(annotations$tribe[match(tree_ASTRALIII_strict_astral_annotations@phylo$tip.label, annotations$Library_ID)]==tribe))==1){
      support_tribal$nuclear_LPP[t]<-"-"
      support_tribal$nuclear_Q1[t]<-"-"
    } else {
      support_tribal$nuclear_LPP[t]<-tree_ASTRALIII_strict_astral_annotations@nodeData[tree_ASTRALIII_strict_astral_annotations@nodeData$node==tribal_node_coal,]$pp1
      support_tribal$nuclear_Q1[t]<-tree_ASTRALIII_strict_astral_annotations@nodeData[tree_ASTRALIII_strict_astral_annotations@nodeData$node==tribal_node_coal,]$q1
    }
  } else {
    support_tribal$nuclear_LPP[t]<-"-"
    support_tribal$nuclear_Q1[t]<-"-"
  }
  
  #Repeat for the plastome phylogeny.
  #Note we first need to check if the tribe is in the phylogeny at all.
  if (length(which(annotations_cp$tribe[match(tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label, annotations_cp$Library_ID)]==tribe))>0){
    #To find nodal support from the ML supermatrix approach, first find the node with the MRCA.
    #Note that in cases where a tribe is represented by a single sample, we cannot know these values, except for the stem age.
    #We then need to take the parent node of the mrca to get the stem age.
    tribal_node_plastome<-MRCA(tree_IQ_TREE_chloroplast_calibrated, which(annotations_cp$tribe[match(tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label, annotations_cp$Library_ID)]==tribe))
    if (length(which(annotations_cp$tribe[match(tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label, annotations_cp$Library_ID)]==tribe))==1){
      support_tribal$plastome_BS[t]<-"-"
      support_tribal$plastome_sCF[t]<-"-"
      support_tribal$plastome_crown_age_median[t]<-"-"
      support_tribal$plastome_crown_age_95CI_lower[t]<-"-"
      support_tribal$plastome_crown_age_95CI_upper[t]<-"-"
      support_tribal$plastome_stem_age_median[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, tribal_node_plastome),]$height_median
      support_tribal$plastome_stem_age_95CI_lower[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, tribal_node_plastome),]$height_0.95_HPD[[1]][1]
      support_tribal$plastome_stem_age_95CI_upper[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, tribal_node_plastome),]$height_0.95_HPD[[1]][2]
      
    } else {
      #Find and save nodal support.
      support_tribal$plastome_BS[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==tribal_node_plastome,]$bootstrap
      support_tribal$plastome_sCF[t]<-strsplit(tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==tribal_node_plastome,]$concordance_factors, "/", 1)[[1]][1]
      #Find and save time calibration data.
      support_tribal$plastome_crown_age_median[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==tribal_node_plastome,]$height_median
      support_tribal$plastome_crown_age_95CI_lower[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==tribal_node_plastome,]$height_0.95_HPD[[1]][1]
      support_tribal$plastome_crown_age_95CI_upper[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==tribal_node_plastome,]$height_0.95_HPD[[1]][2]
      #Repeat this for the stem age, so focusing on the parent of the MRCA.
      support_tribal$plastome_stem_age_median[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, tribal_node_plastome),]$height_median
      support_tribal$plastome_stem_age_95CI_lower[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, tribal_node_plastome),]$height_0.95_HPD[[1]][1]
      support_tribal$plastome_stem_age_95CI_upper[t]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, tribal_node_plastome),]$height_0.95_HPD[[1]][2]
    }
  } else {
    support_tribal$plastome_BS[t]<-"-"
    support_tribal$plastome_sCF[t]<-"-"
    support_tribal$plastome_crown_age_median[t]<-"-"
    support_tribal$plastome_crown_age_95CI_lower[t]<-"-"
    support_tribal$plastome_crown_age_95CI_upper[t]<-"-"
    support_tribal$plastome_stem_age_median[t]<-"-"
    support_tribal$plastome_stem_age_95CI_lower[t]<-"-"
    support_tribal$plastome_stem_age_95CI_upper[t]<-"-"
  }
  
}


#Repeat for the main lineages.
list_main_lineages<-c("basal", "I", "II", "III", "IV", "V")

#Create a dataframe to store results.
support_main_lineages<-data.frame(clade=list_main_lineages,
                                  main_lineage=NA,
                                  nuclear_BS=NA,
                                  nuclear_gCF=NA,
                                  nuclear_sCF=NA,
                                  nuclear_LPP=NA,
                                  nuclear_Q1=NA,
                                  plastome_BS=NA,
                                  plastome_sCF=NA,
                                  nuclear_crown_age_median=NA,
                                  nuclear_crown_age_95CI_lower=NA,
                                  nuclear_crown_age_95CI_upper=NA,
                                  nuclear_stem_age_median=NA,
                                  nuclear_stem_age_95CI_lower=NA,
                                  nuclear_stem_age_95CI_upper=NA,
                                  plastome_crown_age_median=NA,
                                  plastome_crown_age_95CI_lower=NA,
                                  plastome_crown_age_95CI_upper=NA,
                                  plastome_stem_age_median=NA,
                                  plastome_stem_age_95CI_lower=NA,
                                  plastome_stem_age_95CI_upper=NA)

for (ml in 1:length(list_main_lineages)){
  main_lineage<-list_main_lineages[ml]
  support_main_lineages$main_lineage[ml]<-list_main_lineages[ml]
  
  main_lineage_node_ML<-MRCA(tree_IQ_TREE_supermatrix_calibrated, which(annotations$Lineages_Nikolov_et_al[match(tree_IQ_TREE_supermatrix_calibrated@phylo$tip.label, annotations$Library_ID)]==main_lineage))
  support_main_lineages$nuclear_BS[ml]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==main_lineage_node_ML,]$bootstrap
  support_main_lineages$nuclear_gCF[ml]<-strsplit(tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==main_lineage_node_ML,]$concordance_factors, "/", 1)[[1]][1]
  support_main_lineages$nuclear_sCF[ml]<-strsplit(tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==main_lineage_node_ML,]$concordance_factors, "/", 1)[[1]][2]
  support_main_lineages$nuclear_crown_age_median[ml]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==main_lineage_node_ML,]$height_median
  support_main_lineages$nuclear_crown_age_95CI_lower[ml]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==main_lineage_node_ML,]$height_0.95_HPD[[1]][1]
  support_main_lineages$nuclear_crown_age_95CI_upper[ml]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==main_lineage_node_ML,]$height_0.95_HPD[[1]][2]
  support_main_lineages$nuclear_stem_age_median[ml]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, main_lineage_node_ML),]$height_median
  support_main_lineages$nuclear_stem_age_95CI_lower[ml]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, main_lineage_node_ML),]$height_0.95_HPD[[1]][1]
  support_main_lineages$nuclear_stem_age_95CI_upper[ml]<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, main_lineage_node_ML),]$height_0.95_HPD[[1]][2]
  
  main_lineage_node_coal<-MRCA(tree_ASTRALIII_strict_astral_annotations@phylo, which(annotations$Lineages_Nikolov_et_al[match(tree_ASTRALIII_strict_astral_annotations@phylo$tip.label, annotations$Library_ID)]==main_lineage))
  
  support_main_lineages$nuclear_LPP[ml]<-tree_ASTRALIII_strict_astral_annotations@nodeData[tree_ASTRALIII_strict_astral_annotations@nodeData$node==main_lineage_node_coal,]$pp1
  support_main_lineages$nuclear_Q1[ml]<-tree_ASTRALIII_strict_astral_annotations@nodeData[tree_ASTRALIII_strict_astral_annotations@nodeData$node==main_lineage_node_coal,]$q1
  
  #For the plastome phylogeny we'll only consider main lineages basal, I and II, as the others are polyphyletic.
  if(main_lineage %in% c("basal", "I", "III")){
    main_lineage_node_plastome<-MRCA(tree_IQ_TREE_chloroplast_calibrated, which(annotations_cp$Lineages_Nikolov_et_al[match(tree_IQ_TREE_chloroplast_calibrated@phylo$tip.label, annotations_cp$Library_ID)]==main_lineage))
    support_main_lineages$plastome_BS[ml]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==main_lineage_node_plastome,]$bootstrap
    support_main_lineages$plastome_sCF[ml]<-strsplit(tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==main_lineage_node_plastome,]$concordance_factors, "/", 1)[[1]][1]
    support_main_lineages$plastome_crown_age_median[ml]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==main_lineage_node_plastome,]$height_median
    support_main_lineages$plastome_crown_age_95CI_lower[ml]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==main_lineage_node_plastome,]$height_0.95_HPD[[1]][1]
    support_main_lineages$plastome_crown_age_95CI_upper[ml]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==main_lineage_node_plastome,]$height_0.95_HPD[[1]][2]
    support_main_lineages$plastome_stem_age_median[ml]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, main_lineage_node_plastome),]$height_median
    support_main_lineages$plastome_stem_age_95CI_lower[ml]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, main_lineage_node_plastome),]$height_0.95_HPD[[1]][1]
    support_main_lineages$plastome_stem_age_95CI_upper[ml]<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, main_lineage_node_plastome),]$height_0.95_HPD[[1]][2]
  } else {
    support_main_lineages$plastome_BS[ml]<-"-"
    support_main_lineages$plastome_sCF[ml]<-"-"
    support_main_lineages$plastome_crown_age_median[ml]<-"-"
    support_main_lineages$plastome_crown_age_95CI_lower[ml]<-"-"
    support_main_lineages$plastome_crown_age_95CI_upper[ml]<-"-"
    support_main_lineages$plastome_stem_age_median[ml]<-"-"
    support_main_lineages$plastome_stem_age_95CI_lower[ml]<-"-"
    support_main_lineages$plastome_stem_age_95CI_upper[ml]<-"-"
  }
}

#Combine tables for a single overview.
support_tribal<-rbind(support_main_lineages,
                      support_tribal)

#Sort table.
support_tribal<-support_tribal[with(support_tribal, order(main_lineage, clade)),]


##Save table with nodal support for publication.
write.table(support_tribal, file = "11d.Nodal_support_and_ages_tribes_and_main_lineages.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)





#CALIBRATION AGE CHECKS----

##ML nuclear----

#†Akania sp.
node_Akaniaceae_crown_n<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("S1384", "HYZL"))
(age_median_Akaniaceae_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Akaniaceae_crown_n,]$height_median)
(age_95CI_Akaniaceae_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Akaniaceae_crown_n,]$height_0.95_HPD)
(age_median_Akaniaceae_stem_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, node_Akaniaceae_crown_n),]$height_median)
(age_95CI_Akaniaceae_stem_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, node_Akaniaceae_crown_n),]$height_0.95_HPD)

#†Palaeocleome lakensis.
node_Brassicaceae_stem_n<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("PAFTOL_019361", "S1408"))
(age_median_Brassicaceae_stem_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Brassicaceae_stem_n,]$height_median)
(age_95CI_Brassicaceae_stem_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Brassicaceae_stem_n,]$height_0.95_HPD)

#†Capparidoxylon holleisii.
node_Capparis_crown_n<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("S1388", "S1389"))
(age_median_Capparis_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Capparis_crown_n,]$height_median)
(age_95CI_Capparis_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Capparis_crown_n,]$height_0.95_HPD)
(age_median_Capparis_stem_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, node_Capparis_crown_n),]$height_median)
(age_95CI_Capparis_stem_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, node_Capparis_crown_n),]$height_0.95_HPD)

#†Thlaspi primaevum.
node_Thlaspi<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("S1376", "S1419"))
(age_median_Thlaspi_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Thlaspi,]$height_median)
(age_95CI_Thlaspi_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Thlaspi,]$height_0.95_HPD)

#Ricotia.
node_Ricotia_crown_n<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("S0746", "S1074", "S1077", "SRR8528399"))
(age_median_Ricotia_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Ricotia_crown_n,]$height_median)
(age_95CI_Ricotia_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Ricotia_crown_n,]$height_0.95_HPD)

#Arabis.
node_Arabis_crown_n<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("SRR11470320", "S1083"))
(age_median_Arabis_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Arabis_crown_n,]$height_median)
(age_95CI_Arabis_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Arabis_crown_n,]$height_0.95_HPD)

#Pachycladon.
node_Pachycladon_crown_n<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("S0435", "S1400"))
(age_median_Pachycladon_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Pachycladon_crown_n,]$height_median)
(age_95CI_Pachycladon_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Pachycladon_crown_n,]$height_0.95_HPD)

#Clausia aprica.
node_Clausia_crown_n<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("S1187sl", "S1188"))
(age_median_Clausia_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Clausia_crown_n,]$height_median)
(age_95CI_Clausia_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Clausia_crown_n,]$height_0.95_HPD)

#Lepidium.
node_Lepidium_crown_n<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("S0547", "S1561"))
(age_median_Lepidium_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Lepidium_crown_n,]$height_median)
(age_95CI_Lepidium_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Lepidium_crown_n,]$height_0.95_HPD)

#Other
node_Lepidium_crown_n<-MRCA(tree_IQ_TREE_supermatrix_calibrated, c("S0547", "S1223"))
(age_median_Lepidium_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Lepidium_crown_n,]$height_median)
(age_95CI_Lepidium_crown_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==node_Lepidium_crown_n,]$height_0.95_HPD)
(age_median_Capparis_stem_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, node_Lepidium_crown_n),]$height_median)
(age_95CI_Capparis_stem_n<-tree_IQ_TREE_supermatrix_calibrated@data[tree_IQ_TREE_supermatrix_calibrated@data$node==getParent(tree_IQ_TREE_supermatrix_calibrated@phylo, node_Lepidium_crown_n),]$height_0.95_HPD)


##ML plastome----

#†Akania sp.
node_Akaniaceae_stem_cp<-MRCA(tree_IQ_TREE_chloroplast_calibrated, c("S1384", "MK637814"))
(age_median_Akaniaceae_stem_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Akaniaceae_stem_cp,]$height_median)
(age_95CI_Akaniaceae_stem_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Akaniaceae_stem_cp,]$height_0.95_HPD)

#†Palaeocleome lakensis.
node_Brassicaceae_stem_cp<-MRCA(tree_IQ_TREE_chloroplast_calibrated, c("MK637687", "S1408"))
(age_median_Brassicaceae_stem_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Brassicaceae_stem_cp,]$height_median)
(age_95CI_Brassicaceae_stem_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Brassicaceae_stem_cp,]$height_0.95_HPD)

#†Capparidoxylon holleisii.
node_Capparis_crown_cp<-MRCA(tree_IQ_TREE_chloroplast_calibrated, c("S1388", "S1389"))
(age_median_Capparis_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Capparis_crown_cp,]$height_median)
(age_95CI_Capparis_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Capparis_crown_cp,]$height_0.95_HPD)
(age_median_Capparis_stem_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, node_Capparis_crown_cp),]$height_median)
(age_95CI_Capparis_stem_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==getParent(tree_IQ_TREE_chloroplast_calibrated@phylo, node_Capparis_crown_cp),]$height_0.95_HPD)

#†Thlaspi primaevum.
node_Thlaspi<-MRCA(tree_IQ_TREE_chloroplast_calibrated, c("S1376", "S1419"))
(age_median_Thlaspi_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Thlaspi,]$height_median)
(age_95CI_Thlaspi_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Thlaspi,]$height_0.95_HPD)

#Ricotia.
node_Ricotia_crown_cp<-MRCA(tree_IQ_TREE_chloroplast_calibrated, c("S0746", "S1074", "S1077"))
(age_median_Ricotia_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Ricotia_crown_cp,]$height_median)
(age_95CI_Ricotia_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Ricotia_crown_cp,]$height_0.95_HPD)

#Arabis.
node_Arabis_crown_cp<-MRCA(tree_IQ_TREE_chloroplast_calibrated, c("SRR11470320", "S1083"))
(age_median_Arabis_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Arabis_crown_cp,]$height_median)
(age_95CI_Arabis_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Arabis_crown_cp,]$height_0.95_HPD)

#Pachycladon.
node_Pachycladon_crown_cp<-MRCA(tree_IQ_TREE_chloroplast_calibrated, c("S0435", "S1400"))
(age_median_Pachycladon_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Pachycladon_crown_cp,]$height_median)
(age_95CI_Pachycladon_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Pachycladon_crown_cp,]$height_0.95_HPD)

#Clausia aprica.
#note: not both populations (east and west) sampled in phylogeny, as was done in nuclear phylogeny

#Lepidium.
node_Lepidium_crown_cp<-MRCA(tree_IQ_TREE_chloroplast_calibrated, c("NC_009273", "S1561"))
(age_median_Lepidium_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Lepidium_crown_cp,]$height_median)
(age_95CI_Lepidium_crown_cp<-tree_IQ_TREE_chloroplast_calibrated@data[tree_IQ_TREE_chloroplast_calibrated@data$node==node_Lepidium_crown_cp,]$height_0.95_HPD)


#PLOT GENE RECOVERY VERSUS AGE ---- 

gene_recovery<-rbind(
  data.frame(routine="inclusive",
                          collection_year=annotations$collection_date,
                          genes_recovered=annotations$ASTRAL_taxon_occupancy_inclusive),
  data.frame(routine="strict",
             collection_year=annotations$collection_date,
             genes_recovered=annotations$ASTRAL_taxon_occupancy_strict),
  data.frame(routine="superstrict",
             collection_year=annotations$collection_date,
             genes_recovered=annotations$ASTRAL_taxon_occupancy_superstrict),
  data.frame(routine="superstrict by tribe",
             collection_year=annotations$collection_date,
             genes_recovered=annotations$ASTRAL_taxon_occupancy_superstrict_by_tribe)
)

#Plot results.
ggplot(data=gene_recovery[gene_recovery$collection_year!="" & !is.na(gene_recovery$genes_recovered),], aes(x=as.numeric(collection_year), y=genes_recovered)) +
  geom_point() +
  facet_wrap(.~routine)+
  geom_smooth(method='lm', formula= y~x)+
  xlab("Collection year") +
  ylab("Number of genes recovered")+
  theme_classic()



#PLOT COPHYLOGENY ---- 

##Load metadata ----
metadata<-read.csv("2e.metadata.csv", header=T, stringsAsFactors=FALSE, fileEncoding="latin1")

##Species trees ----

#Load tree files.

#Load nuclear species tree.
tree_n<-read.beast(file = "10f.treePL_supermatrix_approach_calibrated.tree")
tree_n<-tree_n@phylo

#Load chloroplast species tree.
#note the following tree was calibrated
tree_cp<-read.beast(file = "7.chloroplast_treePL_BS_analysis_TreeAnnotator_EDITED.tre")
tree_cp<-tree_cp@phylo

##Prune trees to tribal level.

#For the nuclear species tree.

#Get a vector of all tribe names corresponding to the tip labels.
tree_n_tribes<-metadata$Tribe[match(tree_n$tip.label, metadata$Library_ID)]
#First prune all tips that are outgroups, except Cleomaceae, the sister family to the Brassicaceae.
tips_to_drop_outgroups_n<-tree_n_tribes[grepl("Outgroup_", tree_n_tribes) & !grepl("Outgroup_Cleomaceae", tree_n_tribes)]
tips_to_drop_outgroups_n<-which(metadata$Tribe[match(tree_n$tip.label, metadata$Library_ID)] %in% tips_to_drop_outgroups_n)
#Now drop outgroup tips.
tree_n<-drop.tip(tree_n, tip = tips_to_drop_outgroups_n)
#Also drop tips if from tribe "Brassiceae II", which is not a formal tribe but used here for polyphyly of the Brassiceae.
tree_n<-drop.tip(tree_n, tip = which(metadata$Tribe[match(tree_n$tip.label, metadata$Library_ID)] == "Brassiceae II"))
#Now we will subset to a single, randomly chosen representative for each tribe.
tips_unique_tribes_to_keep_n<-unique(ave(seq_along(metadata$Tribe[match(tree_n$tip.label, metadata$Library_ID)] ), metadata$Tribe[match(tree_n$tip.label, metadata$Library_ID)] , FUN = function(x) if(length(x) > 1) head(sample(x), 1) else x))
tips_unique_tribes_to_drop_n<-which(!(1:length(tree_n$tip.label) %in% tips_unique_tribes_to_keep_n))
tree_n<-drop.tip(tree_n, tip = tips_unique_tribes_to_drop_n)
#Finally, change tip label names to tribe names.
tree_n$tip.label<-metadata$Tribe[match(tree_n$tip.label, metadata$Library_ID)]

#For the chloroplast species tree.
#Get a vector of all tribe names corresponding to the tip labels.
tree_cp_tribes<-metadata$Tribe[match(tree_cp$tip.label, metadata$Library_ID)]
#First prune all tips that are outgroups, except Cleomaceae, the sister family to the Brassicaceae.
tips_to_drop_outgroups_cp<-tree_cp_tribes[grepl("Outgroup_", tree_cp_tribes) & !grepl("Outgroup_Cleomaceae", tree_cp_tribes) | is.na(tree_cp_tribes)]
tips_to_drop_outgroups_cp<-which(metadata$Tribe[match(tree_cp$tip.label, metadata$Library_ID)] %in% tips_to_drop_outgroups_cp)
#Now drop outgroup tips.
tree_cp<-drop.tip(tree_cp, tip = tips_to_drop_outgroups_cp)
#Also drop tips if from tribe "Brassiceae II", which is not a formal tribe but used here for polyphyly of the Brassiceae.
tree_cp<-drop.tip(tree_cp, tip = which(metadata$Tribe[match(tree_cp$tip.label, metadata$Library_ID)] == "Brassiceae II"))
#Now we will subset to a single, randomly chosen representative for each tribe.
tips_unique_tribes_to_keep_cp<-unique(ave(seq_along(metadata$Tribe[match(tree_cp$tip.label, metadata$Library_ID)] ), metadata$Tribe[match(tree_cp$tip.label, metadata$Library_ID)] , FUN = function(x) if(length(x) > 1) head(sample(x), 1) else x))
tips_unique_tribes_to_drop_cp<-which(!(1:length(tree_cp$tip.label) %in% tips_unique_tribes_to_keep_cp))
tree_cp<-drop.tip(tree_cp, tip = tips_unique_tribes_to_drop_cp)
#Finally, change tip label names to tribe names.
tree_cp$tip.label<-metadata$Tribe[match(tree_cp$tip.label, metadata$Library_ID)]

#Both trees need to have exactly the same tip labels, so we need to further prune them.
tribes_in_both_trees<-intersect(tree_n$tip.label, tree_cp$tip.label)
tree_n<-drop.tip(tree_n, tip = which(!tree_n$tip.label %in% tribes_in_both_trees))
tree_cp<-drop.tip(tree_cp, tip = which(!tree_cp$tip.label %in% tribes_in_both_trees))

#Remove any spaces in tip label names with underscores for the below script to work properly.
tree_n$tip.label[which(grepl(" trib. nov.", tree_n$tip.label))]<-paste0(sapply(strsplit(tree_n$tip.label, " "), "[", 1), "_",sapply(strsplit(tree_n$tip.label, " "), "[", 2), "_", sapply(strsplit(tree_n$tip.label, " "), "[", 3))[which(grepl(" trib. nov.", tree_n$tip.label))]
tree_n$tip.label[which(grepl(" ", tree_n$tip.label))]<-paste0(sapply(strsplit(tree_n$tip.label, " "), "[", 1), "_",sapply(strsplit(tree_n$tip.label, " "), "[", 2))[which(grepl(" ", tree_n$tip.label))]
tree_cp$tip.label[which(grepl(" trib. nov.", tree_cp$tip.label))]<-paste0(sapply(strsplit(tree_cp$tip.label, " "), "[", 1), "_",sapply(strsplit(tree_cp$tip.label, " "), "[", 2), "_", sapply(strsplit(tree_cp$tip.label, " "), "[", 3))[which(grepl(" trib. nov.", tree_cp$tip.label))]
tree_cp$tip.label[which(grepl(" ", tree_cp$tip.label))]<-paste0(sapply(strsplit(tree_cp$tip.label, " "), "[", 1), "_",sapply(strsplit(tree_cp$tip.label, " "), "[", 2))[which(grepl(" ", tree_cp$tip.label))]

#Check that exactly the same tip labels are now in both tree.
#The following line should result only in TRUE values.
tree_n$tip.label %in% tree_cp$tip.label

# tree_n_untangle<-untangle(tree_n, "reorder")

#Reorder tips to be consitent with the descending order in FigTree and ggtree used above.
tip.order.n<-c("Outgroup_Cleomaceae","Aethionemeae","Dontostemoneae","Chorisporeae","Shehbazieae","Hesperideae","Euclidieae","Anchonieae","Buniadeae","Arabideae","Stevenieae","Alysseae","Asperuginoideae_trib._nov.","Cardamineae","Lepidieae","Descurainieae","Smelowskieae","Yinshanieae","Erysimeae","Malcolmieae","Physarieae","Arabidopsideae_trib._nov.","Oreophytoneae","Alyssopsideae","Turritideae","Camelineae_I","Microlepidieae","Hemilophieae_trib._nov.","Crucihimalayeae","Boechereae","Halimolobeae","Megacarpaeeae","Biscutelleae","Anastaticeae","Iberideae_I","Notothlaspideae","Heliophileae","Chamireae","Iberideae_II","Subularieae","Asteae","Eudemeae","Cremolobeae","Schizopetaleae","Aphragmeae","Kernereae","Cochlearieae","Conringieae","Coluteocarpeae","Plagiolobeae","Calepineae","Thlaspideae","Eutremeae","Schrenkielleae_trib._nov.","Fourraeeae","Brassiceae","Isatideae","Thelypodieae","Sisymbrieae")
tree_n<-ape::rotateConstr(tree_n, tip.order.n)

tip.order.cp<-c("Outgroup_Cleomaceae","Aethionemeae","Descurainieae","Smelowskieae","Yinshanieae","Lepidieae","Cardamineae","Stevenieae","Physarieae","Halimolobeae","Boechereae","Hemilophieae_trib._nov.","Crucihimalayeae","Microlepidieae","Alyssopsideae","Erysimeae","Malcolmieae","Turritideae","Oreophytoneae","Camelineae_I","Arabidopsideae_trib._nov.","Dontostemoneae","Shehbazieae","Chorisporeae","Buniadeae","Hesperideae","Anchonieae","Euclidieae","Asperuginoideae_trib._nov.","Biscutelleae","Cochlearieae","Iberideae_II","Iberideae_I","Anastaticeae","Megacarpaeeae","Aphragmeae","Conringieae","Coluteocarpeae","Plagiolobeae","Notothlaspideae","Chamireae","Heliophileae","Schizopetaleae","Kernereae","Asteae","Subularieae","Cremolobeae","Eudemeae","Alysseae","Arabideae","Thlaspideae","Calepineae","Eutremeae","Schrenkielleae_trib._nov.","Fourraeeae","Isatideae","Sisymbrieae","Thelypodieae","Brassiceae")
tree_cp<-ape::rotateConstr(tree_cp, tip.order.cp)

#Create the cophylogeny using phytools.
brassitol_cophylo<-cophylo(tree_n,
                           tree_cp,
                           rotate = F)

#Open the pdf printer.
pdf(file = "12.brassitol_cophylogeny.pdf",
    width = 14,
    height = 10)

#Plot the cophylogeny.
plot(brassitol_cophylo,
     link.type="curved",
     link.lwd=4,
     link.lty="solid",
     link.col=make.transparent("red", 0.25))

#Close the pdf printer.
dev.off()



#PLOT BRANCH LENGTH AND SUPPORT VERSUS AGE ---- 

##Load metadata
metadata<-read.csv("2e.metadata.csv", header=T, stringsAsFactors=FALSE, fileEncoding="latin1")

#Load again tree files.

#Load nuclear species tree.
#note the following tree was not calibrated
# tree_n<-read.tree(file = "9b.IQ-TREE_supermatrix_approach.tree")
tree_n<-tree_IQ_TREE_supermatrix_calibrated
tree_n@data$mean_edge_height<-tree_n@data$height-0.5*tree_n@data$length

#Load chloroplast species tree.
#note the following tree was calibrated
tree_cp<-tree_IQ_TREE_chloroplast_calibrated
tree_cp@data$mean_edge_height<-tree_cp@data$height-0.5*tree_cp@data$length

#Create a dataframe that lists all edges from the two BrassiToLs with their height (i.e., age) and support.
brassitol_edges_sCF<-data.frame(tree=c(rep("nuclear", length(tree_n@data$height)),
                                   rep("plastome", length(tree_cp@data$height))),
                                metric=c(rep("sCF", length(tree_n@data$height)),
                                       rep("sCF", length(tree_cp@data$height))),
                                edge_height_mean=c(tree_n@data$mean_edge_height,
                                                   tree_cp@data$mean_edge_height),
                                edge_length=c(tree_n@data$length,
                                              tree_cp@data$length),
                                support=c(as.numeric(tree_n@data$sCF),
                                          as.numeric(tree_cp@data$sCF)))
brassitol_edges_gCF<-data.frame(tree=c(rep("nuclear", length(tree_n@data$height))),
                                metric=c(rep("gCF", length(tree_n@data$height))),
                                edge_height_mean=tree_n@data$mean_edge_height,
                                edge_length=tree_n@data$length,
                                support=as.numeric(tree_n@data$gCF))
brassitol_edges<-rbind(brassitol_edges_sCF,
                       brassitol_edges_gCF)

#Remove any edges with a negative length; these are the terminal branches for which mean height could not be calculated in the above
#way, and for which there will be no node support values.
brassitol_edges<-brassitol_edges[brassitol_edges$edge_height_mean>0, ]

#Highlight edges associated with the onset of the Brassicaceae family and it's supertribes.
#This is between 30 and 14 mya for the nuclear dataset.
brassitol_edges$family_onset<-NA
brassitol_edges$family_onset[brassitol_edges$tree=="nuclear" & 
                               brassitol_edges$edge_height_mean>14 & 
                               brassitol_edges$edge_height_mean<30]<- 
  "family_onset"

#This is between 25 and 8 mya for the plastome dataset.
brassitol_edges$family_onset[brassitol_edges$tree=="plastome" & 
                               brassitol_edges$edge_height_mean>8 & 
                               brassitol_edges$edge_height_mean<25]<- 
  "family_onset"

#Name the other edges, too.
brassitol_edges$family_onset[is.na(brassitol_edges$family_onset)]<-"other"


#Plot results in an easy to interpret way.
brassitol_edge_height_support<-ggplot(data=brassitol_edges, aes(x=-1*edge_height_mean, y=support))+
  geom_point(size=2, shape = 21, colour = "black", fill = "white")+
  geom_point(data=brassitol_edges[brassitol_edges$family_onset=="family_onset",], size=2, shape = 21, colour = "black", fill = "black")+
  facet_grid(rows = vars(metric),
             cols = vars(tree),
             scales='free')+
  ylim(c(0,100))+
  stat_summary(fun.data=mean_cl_normal, size = 0.4, color = "transparent") +
  geom_smooth(method='lm', formula= y~x, color="darkgreen", size=0.4)+
  theme_bw()

brassitol_edge_length_support<-ggplot(data=brassitol_edges, aes(x=edge_length, y=support))+
  geom_point(size=2, shape = 21, colour = "black", fill = "white")+
  geom_point(data=brassitol_edges[brassitol_edges$family_onset=="family_onset",], size=2, shape = 21, colour = "black", fill = "black")+
  facet_grid(rows = vars(metric),
             cols = vars(tree),
             scales='free')+
  ylim(c(0,100))+
  stat_summary(fun.data=mean_cl_normal, size = 0.4, color = "transparent") +
  geom_smooth(method='lm', formula= y~x, color="darkgreen", size=0.4)+
  theme_bw()

brassitol_edge_support<-plot_grid(brassitol_edge_height_support,
          brassitol_edge_length_support,
          nrow=2,
          ncol=1,
          labels= c("A", "B"))

#Save plots for publication.
ggsave(filename = "13.brassitol_edge_support.pdf", 
       brassitol_edge_support + theme(legend.position="none"),
       width = 8, height = 10)


#TREE SIMILARITY NUCLEAR VS. PLASTOME ---- 

##Load metadata
metadata<-read.csv("2e.metadata.csv", header=T, stringsAsFactors=FALSE, fileEncoding="latin1")

#Load again tree files.
tree_n<-tree_IQ_TREE_supermatrix_calibrated@phylo
tree_cp<-tree_IQ_TREE_chloroplast_calibrated@phylo

#Get all tribes present in the two BrassiToLs.
tribes<-
  sort(unique(metadata$Tribe[match(unique(c(tree_n$tip.label, tree_cp$tip.label)),
                     metadata$Library_ID)]))
#Remove 'tribes' from outgroups and tribes with numbers because polyphyletic.
tribes<-
  tribes[!grepl("Outgroup_", tribes)]
tribes<-
  tribes[!tribes %in% c("Brassiceae I", "Brassiceae II", "Camelineae I", "Camelineae III", "Iberideae I", "Iberideae II")]

#Create a dataframe to store results.
brassitol_tree_distances_n_vs_cp<-
  data.frame(Tribe=tribes,
             tips_n=NA,
             tips_cp=NA,
             tips_shared=NA,
             dist_RF_classic=NA,
             dist_RF_gen=NA)

#Loop through the tribes.

for (t in 1:length(tribes)){
  #Select the t-th tribe.
  tribe<-
    tribes[t]

  #Subset both trees to the tips that belong to a tribe.
  tree_n_tribe<-
    keep.tip(tree_n,
           which(metadata$Tribe[match(tree_n$tip.label, metadata$Library_ID)]==tribe)
    )
  tree_cp_tribe<-
    keep.tip(tree_cp,
             which(metadata$Tribe[match(tree_cp$tip.label, metadata$Library_ID)]==tribe)
    )
  
  #Store number of tips.
  tree_n_tribe_no_tips<-
    length(tree_n_tribe$tip.label)
  tree_cp_tribe_no_tips<-
    length(tree_cp_tribe$tip.label)
  
  #Rename tip labels to the genus.
  tree_n_tribe$tip.label<-
    metadata$Name[match(tree_n_tribe$tip.label, metadata$Library_ID)]
  tree_cp_tribe$tip.label<-
    metadata$Name[match(tree_cp_tribe$tip.label, metadata$Library_ID)]
  
  #Find the species present in both BrassiToLs.
  species_matching<-
    intersect(tree_n_tribe$tip.label, tree_cp_tribe$tip.label)
  
  #Only continue if number of shared species is > 3.
  if(length(species_matching)!=0){
    if(length(species_matching)>3){
    #Subset both BrassiToLs to the matching species.
    tree_n_tribe<-
      keep.tip(tree_n_tribe,
               which(tree_n_tribe$tip.label %in% species_matching))
    tree_cp_tribe<-
      keep.tip(tree_cp_tribe,
               which(tree_cp_tribe$tip.label %in% species_matching))
    
    #Remove duplicates by removing a random tip from each duplicate.
    tree_n_tribe<-
      drop.tip(tree_n_tribe,
               which(duplicated(tree_n_tribe$tip.label)))
    tree_cp_tribe<-
      drop.tip(tree_cp_tribe,
               which(duplicated(tree_cp_tribe$tip.label)))
    
    #Calculate two useful distance metrics.
    dist_RF_classic<-TreeDist::RobinsonFoulds(tree_n_tribe,
                                              tree_cp_tribe)
    dist_RF_gen<-TreeDist::TreeDistance(tree_n_tribe,
                                        tree_cp_tribe)
    }
  }
  
  #Store results in dataframe.
  brassitol_tree_distances_n_vs_cp$Tribe[t]<-
    tribe
  brassitol_tree_distances_n_vs_cp$tips_n[t]<-
    tree_n_tribe_no_tips
  brassitol_tree_distances_n_vs_cp$tips_cp[t]<-
    tree_cp_tribe_no_tips
  brassitol_tree_distances_n_vs_cp$tips_shared[t]<-
    length(species_matching)
  if(length(species_matching)!=0){
    if(length(species_matching)>3){
    brassitol_tree_distances_n_vs_cp$dist_RF_classic[t]<-
      dist_RF_classic
    brassitol_tree_distances_n_vs_cp$dist_RF_gen[t]<-
      dist_RF_gen
    }
  }
}

#Save table with results for publication.
write.table(brassitol_tree_distances_n_vs_cp,
            file = "14.brassitol_tree_distances_n_vs_cp.csv")

#Calculate some summary statistics for publication.
mean(brassitol_tree_distances_n_vs_cp$dist_RF_gen,
     na.rm = T)
median(brassitol_tree_distances_n_vs_cp$dist_RF_gen,
     na.rm = T)
min(brassitol_tree_distances_n_vs_cp$dist_RF_gen,
    na.rm = T)
max(brassitol_tree_distances_n_vs_cp$dist_RF_gen,
    na.rm = T)
sum(!(is.na(brassitol_tree_distances_n_vs_cp$dist_RF_gen)))

