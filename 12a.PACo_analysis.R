#With this script we perform a PACo analysis to test for differences between the nuclear and chloroplast phylogenies.

#From step 2 onwards, much of the scripts was borrowed from the "rumbling orchids" script of 
#Pérez-Escobar et al. (2016) Syst. Biol.


#STEP 0: LOAD LIBRARIES ----
library(ape)
library(treeio)
library(cluster)
library(gplots)
library(phytools)
library(vegan)
library(ggplot2)
library(ggtree)
library(dplyr)
library(cowplot)
library(svglite)
library(reshape2)


#STEP 1: LOAD DATA ----

##Load metadata ----
metadata<-read.csv("2e.metadata.csv", header=T, stringsAsFactors=FALSE, fileEncoding="latin1")

##Species trees ----

###Load tree files ----

#Load nuclear species tree.
#note the following tree was not calibrated
tree_n<-read.tree(file = "9b.IQ-TREE_supermatrix_approach.tree")

#Load chloroplast species tree.
#note the following tree was calibrated
tree_cp<-read.tree(file = "results_PACo_analysis/input_species_trees/iqtree_species_n.tree")
# tree_cp<-tree_cp@phylo


###Prune trees to tribal level ----

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


##Bootstrap species trees ----
#note: we follow the same routine for the bootstrapped trees below, but now leave out the comments

###Nuclear bootstrap trees ----

#Create a trees object to store pruned bootstrap trees in.
tree_n_bs_pruned<-NA
class(tree_n_bs_pruned)<-"multiPhylo"

tree_n_bs<-read.tree(file = "results_PACo_analysis/input_bootstrap_trees/iqtree_bootstrap_nuclear.trees")
# tree_n_bs<-tree_n_bs[1:25]
for (j in 1:length(tree_n_bs)){
  tree_n_bs_j<-tree_n_bs[[j]]
  #now continue the pruning etc as above for species trees
  tree_n_tribes_j<-metadata$Tribe[match(tree_n_bs_j$tip.label, metadata$Library_ID)]
  tips_to_drop_outgroups_n_j<-tree_n_tribes_j[grepl("Outgroup_", tree_n_tribes_j) & !grepl("Outgroup_Cleomaceae", tree_n_tribes_j) | is.na(tree_n_tribes_j)]
  tips_to_drop_outgroups_n_j<-which(metadata$Tribe[match(tree_n_bs_j$tip.label, metadata$Library_ID)] %in% tips_to_drop_outgroups_n_j)
  tree_n_bs_j<-drop.tip(tree_n_bs_j, tip = tips_to_drop_outgroups_n_j)
  tree_n_bs_j<-drop.tip(tree_n_bs_j, tip = which(metadata$Tribe[match(tree_n_bs_j$tip.label, metadata$Library_ID)] == "Brassiceae II"))
  tips_unique_tribes_to_keep_n_j<-unique(ave(seq_along(metadata$Tribe[match(tree_n_bs_j$tip.label, metadata$Library_ID)] ), metadata$Tribe[match(tree_n_bs_j$tip.label, metadata$Library_ID)] , FUN = function(x) if(length(x) > 1) head(sample(x), 1) else x))
  tips_unique_tribes_to_drop_n_j<-which(!(1:length(tree_n_bs_j$tip.label) %in% tips_unique_tribes_to_keep_n_j))
  tree_n_bs_j<-drop.tip(tree_n_bs_j, tip = tips_unique_tribes_to_drop_n_j)
  tree_n_bs_j$tip.label<-metadata$Tribe[match(tree_n_bs_j$tip.label, metadata$Library_ID)]
  #again prune to the tribes that are found in both species trees
  tree_n_bs_j<-drop.tip(tree_n_bs_j, tip = which(!tree_n_bs_j$tip.label %in% tribes_in_both_trees))
  #and again replace any spaces with underscores
  tree_n_bs_j$tip.label[which(grepl(" trib. nov.", tree_n_bs_j$tip.label))]<-paste0(sapply(strsplit(tree_n_bs_j$tip.label, " "), "[", 1), "_",sapply(strsplit(tree_n_bs_j$tip.label, " "), "[", 2), "_", sapply(strsplit(tree_n_bs_j$tip.label, " "), "[", 3))[which(grepl(" trib. nov.", tree_n_bs_j$tip.label))]
  tree_n_bs_j$tip.label[which(grepl(" ", tree_n_bs_j$tip.label))]<-paste0(sapply(strsplit(tree_n_bs_j$tip.label, " "), "[", 1), "_",sapply(strsplit(tree_n_bs_j$tip.label, " "), "[", 2))[which(grepl(" ", tree_n_bs_j$tip.label))]
  #store pruned bootstrap tree
  tree_n_bs_pruned<-c(tree_n_bs_pruned, tree_n_bs_j)
}
#Remove first NA tree from multiPhylo object.
tree_n_bs_pruned<-tree_n_bs_pruned[2:length(tree_n_bs_pruned)]


###Chloroplast bootstrap trees ----

#Create a trees object to store pruned bootstrap trees in.
tree_cp_bs_pruned<-NA
class(tree_cp_bs_pruned)<-"multiPhylo"

tree_cp_bs<-read.tree(file = "results_PACo_analysis/input_bootstrap_trees/iqtree_bootstrap_plastome.trees")
# tree_cp_bs<-tree_cp_bs[1:25]
for (j in 1:length(tree_cp_bs)){
  #first rename some samples to correct names in metadata file
  tree_cp_bs_j<-tree_cp_bs[[j]]
  tree_cp_bs_j$tip.label[which(grepl("_2", tree_cp_bs_j$tip.label))]<-paste0(substring(tree_cp_bs_j$tip.label[which(grepl("_2", tree_cp_bs_j$tip.label))],1,5), "-", 2)
  tree_cp_bs_j$tip.label[tree_cp_bs_j$tip.label=="KY126841.1_Arabis_stelleri_japonica"]<-"KY126841.1"
  tree_cp_bs_j$tip.label[tree_cp_bs_j$tip.label=="doi_10_5061_dryad_kc5b8_BM"]<-"doi:10.5061/dryad.kc5b8_BM"
  tree_cp_bs_j$tip.label[tree_cp_bs_j$tip.label=="doi_10_5061_dryad_kc5b8_Ssp"]<-"doi:10.5061/dryad.kc5b8_Ssp"
  tree_cp_bs_j$tip.label[tree_cp_bs_j$tip.label=="doi_10_5061_dryad_kc5b8_RO"]<-"doi:10.5061/dryad.kc5b8_RO"
  tree_cp_bs_j$tip.label[tree_cp_bs_j$tip.label=="MK637676_Clausia_aprica_western"]<-"MK637676"
  tree_cp_bs_j$tip.label[tree_cp_bs_j$tip.label=="LN877375_Clausia_aprica_eastern"]<-"LN877375"
  tree_cp_bs_j$tip.label[tree_cp_bs_j$tip.label=="doi_10_5061_dryad_kc5b8_RO"]<-"doi:10.5061/dryad.kc5b8_RO"
  #now continue the pruning etc as above for species trees
  tree_cp_tribes_j<-metadata$Tribe[match(tree_cp_bs_j$tip.label, metadata$Library_ID)]
  tips_to_drop_outgroups_cp_j<-tree_cp_tribes_j[grepl("Outgroup_", tree_cp_tribes_j) & !grepl("Outgroup_Cleomaceae", tree_cp_tribes_j) | is.na(tree_cp_tribes_j)]
  tips_to_drop_outgroups_cp_j<-which(metadata$Tribe[match(tree_cp_bs_j$tip.label, metadata$Library_ID)] %in% tips_to_drop_outgroups_cp_j)
  tree_cp_bs_j<-drop.tip(tree_cp_bs_j, tip = tips_to_drop_outgroups_cp_j)
  tree_cp_bs_j<-drop.tip(tree_cp_bs_j, tip = which(metadata$Tribe[match(tree_cp_bs_j$tip.label, metadata$Library_ID)] == "Brassiceae II"))
  tips_unique_tribes_to_keep_cp_j<-unique(ave(seq_along(metadata$Tribe[match(tree_cp_bs_j$tip.label, metadata$Library_ID)] ), metadata$Tribe[match(tree_cp_bs_j$tip.label, metadata$Library_ID)] , FUN = function(x) if(length(x) > 1) head(sample(x), 1) else x))
  tips_unique_tribes_to_drop_cp_j<-which(!(1:length(tree_cp_bs_j$tip.label) %in% tips_unique_tribes_to_keep_cp_j))
  tree_cp_bs_j<-drop.tip(tree_cp_bs_j, tip = tips_unique_tribes_to_drop_cp_j)
  tree_cp_bs_j$tip.label<-metadata$Tribe[match(tree_cp_bs_j$tip.label, metadata$Library_ID)]
  #again prune to the tribes that are found in both species trees
  tree_cp_bs_j<-drop.tip(tree_cp_bs_j, tip = which(!tree_cp_bs_j$tip.label %in% tribes_in_both_trees))
  #and again replace any spaces with underscores
  tree_cp_bs_j$tip.label[which(grepl(" trib. nov.", tree_cp_bs_j$tip.label))]<-paste0(sapply(strsplit(tree_cp_bs_j$tip.label, " "), "[", 1), "_",sapply(strsplit(tree_cp_bs_j$tip.label, " "), "[", 2), "_", sapply(strsplit(tree_cp_bs_j$tip.label, " "), "[", 3))[which(grepl(" trib. nov.", tree_cp_bs_j$tip.label))]
  tree_cp_bs_j$tip.label[which(grepl(" ", tree_cp_bs_j$tip.label))]<-paste0(sapply(strsplit(tree_cp_bs_j$tip.label, " "), "[", 1), "_",sapply(strsplit(tree_cp_bs_j$tip.label, " "), "[", 2))[which(grepl(" ", tree_cp_bs_j$tip.label))]
  #store pruned bootstrap tree
  tree_cp_bs_pruned<-c(tree_cp_bs_pruned, tree_cp_bs_j)
}
#Remove first NA tree from multiPhylo object.
tree_cp_bs_pruned<-tree_cp_bs_pruned[2:length(tree_cp_bs_pruned)]


#STEP 2: SET UP PACO ANALYSIS ----

##Load functions ----

# Load functions PACo and D.wrapper
# correspondences between the trees
# Distance transformations executed by using deVienne et al. 2011 method 
PACo.dV <- function (H.dist, P.dist, HP.bin) { 
  HP.bin <- which(HP.bin > 0, arr.in=TRUE)
  H.PCo <- pcoa(sqrt(H.dist), correction="none")$vectors 
  P.PCo <- pcoa(sqrt(P.dist), correction="none")$vectors 
  H.PCo <- H.PCo[HP.bin[,1],] 
  P.PCo <- P.PCo[HP.bin[,2],] 
  list (H.PCo = H.PCo, P.PCo = P.PCo)
}

# D.wrapper reads 1 nexus Bayesian tree & executes PACo and Parafit
D.wrapper <- function(n) {
  DH.add <- cophenetic(treeH[[n]]) 
  DP.add <- cophenetic(treeP[[n]]) 
  DH.top <- cophenetic(compute.brlen(treeH[[n]], 1))   
  DP.top <- cophenetic(compute.brlen(treeP[[n]], 1))
  DH.add <- DH.add[rownames(NCP),rownames(NCP)]
  DP.add <- DP.add[colnames(NCP), colnames(NCP)]
  DH.top <- DH.top[rownames(NCP),rownames(NCP)]
  DP.top <- DP.top[colnames(NCP), colnames(NCP)]
  
  PACo.add <- PACo.dV(DH.add, DP.add, HP)
  Proc.add <- procrustes(PACo.add$H.PCo, PACo.add$P.PCo) 
  add.res <- residuals(Proc.add)
  HostX <- Proc.add$X
  ParY <- Proc.add$Yrot
  colnamesPACo <- paste(rownames(HostX),rownames(ParY), sep="_")
  
  PACo.top <- PACo.dV(DH.top, DP.top, HP)
  Proc.top <- procrustes(PACo.top$H.PCo, PACo.top$P.PCo) 
  top.res <- residuals(Proc.top)
  
  PF.add <- parafit(sqrt(DH.add), sqrt(DP.add), HP, nperm=1, test.links=TRUE, silent=TRUE)
  PFL2.add <- c(PF.add$link.table[,5])
  
  PF.top <- parafit(sqrt(DH.top), sqrt(DP.top), HP, nperm=1, test.links=TRUE, silent=TRUE)
  PFL2.top <- c(PF.top$link.table[,5])
  
  # 2.1 Save PACo and Parafit outputs for further use 
  write (add.res, file="results_PACo_analysis/PACo_output_text_files/PACo_res_add.txt", ncolumns = NLinks , append=TRUE, sep="\t")
  write (top.res, file="results_PACo_analysis/PACo_output_text_files/PACo_res_top.txt", ncolumns = NLinks , append=TRUE, sep="\t")
  write (PFL2.add, file="results_PACo_analysis/PACo_output_text_files/PFL2_add.txt", ncolumns = NLinks , append=TRUE, sep="\t")
  write (PFL2.top, file="results_PACo_analysis/PACo_output_text_files/PFL2_top.txt", ncolumns = NLinks , append=TRUE, sep="\t")
  write (colnamesPACo, "results_PACo_analysis/PACo_output_text_files/colnamesPACo.txt", ncolumns=NLinks, sep="\t")
}

##Prepare tree and matrix ----

#Create PACo-ready tree objects..
NTree  <-  tree_n
CPTree  <- tree_cp

NTaxa  <- sort(NTree$tip.label)
CPTaxa <- sort(CPTree$tip.label)
NCP <- as.matrix(table(NTaxa, CPTaxa))


#STEP 3: RUN PACO ----

## 3.1. Compute patristic distances from consensus trees ----

N.D <- cophenetic (NTree)
CP.D <- cophenetic (CPTree)

## 3.2. Sort N and CP taxa in distance matrices to match the N-CP matrix: ----

N.D <- N.D[rownames(NCP),rownames(NCP)]
CP.D <- CP.D [colnames(NCP), colnames(NCP)]


## 3.3. Apply PACo function ----

PACo.fit <- PACo.dV(N.D, CP.D, NCP)
NCP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo) 
NLinks = sum(NCP) 

## 3.4. Test of phylogenetic congruence between N and CP trees ----

m2.obs <- NCP.proc$ss 
N.perm = 1000 
P.value = 0
set.seed(5) 
for (n in c(1:N.perm))
{
  if (NLinks <= nrow(NCP) | NLinks <= ncol(NCP))    
  {  flag2 <- TRUE 
  while (flag2 == TRUE)  { 
    NCP.perm <- t(apply(NCP,1,sample))
    if(any(colSums(NCP.perm) == NLinks)) flag2 <- TRUE else flag2 <- FALSE
  }  
  } else { NCP.perm <- t(apply(NCP,1,sample))} 
  PACo.perm <- PACo.dV(N.D, CP.D, NCP.perm)
  m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss 
  if (m2.perm <= m2.obs)
  {P.value = P.value + 1} 
}
P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\n", "P-value = ", P.value, " based on ", N.perm," permutations.")

#STEP 4: LOAD BOOTSTRAP TREES FOR ID OF OUTLIERS ----

treeH <- tree_n_bs_pruned
treeP <- tree_cp_bs_pruned

# 4.1. drop burn-in sets:
#note: this is not relevant in our analysis in which we use bootstrap replicates, not Bayesian treespace search results
# treeH <- treeH[19000: length(treeH)]
# treeP <- treeP[19000: length(treeP)]

NLinks = sum(NCP)
HP <- NCP


#STEP 5: RUN PACO AND PARAFIT ANALYSES ----

for (k in 1:length(treeH)){
  D.wrapper(k)
}

#note: original line of code below does not work for our bootstrap tree objects
# lapply(1:length(treeH), D.wrapper)


#STEP 6: PLOT PACO SQ. RES. AND PFL2 STATS ----

## 6.1. Read PACo and Parafit output files ----

colnamesPACo <- read.table(file="results_PACo_analysis/PACo_output_text_files/colnamesPACo.txt", header=TRUE)
colnamesPACo <- colnames(colnamesPACo)

#note: since we study tribe-by-tribe correspondence, we update the colnames to simpy the tribe
colnamesPACo<-sapply(colnamesPACo, function(x) substring(x, 1, nchar(x)/2))
names(colnamesPACo)<-NULL
#We can now also replace underscores with spaces again; run twice to replace also second underscores.
colnamesPACo<-sub("_", " ", colnamesPACo)
colnamesPACo<-sub("_", " ", colnamesPACo)

pac.add <- read.table(file="results_PACo_analysis/PACo_output_text_files/PACo_res_add.txt", header=FALSE, col.names=colnamesPACo)
pac.top <- read.table(file="results_PACo_analysis/PACo_output_text_files/PACo_res_top.txt", header=FALSE, col.names=colnamesPACo)
pf2.add <- read.table(file="results_PACo_analysis/PACo_output_text_files/PFL2_add.txt", header=FALSE, col.names=colnamesPACo)
pf2.top <- read.table(file="results_PACo_analysis/PACo_output_text_files/PFL2_top.txt", header=FALSE, col.names=colnamesPACo)

## 6.2. Normalize PACo sq residuals to total m2 ----

m2A <- apply(pac.add, 1, sum)
pac.norm.add <- pac.add/m2A 

m2T <- apply(pac.top, 1, sum)
pac.norm.top <- pac.top/m2T  

## 6.3. Plotting parameters ----

op <- par(oma=c(3,2,1,1))
par (mfrow=c(2,1),mar = c(4,4,1,1)) 

## 6.4. Plot PACo results using phylograms ----

mA <- apply(pac.norm.add, 2, median) ### OBS! ##### MEDIAN
uCI.A <- apply(pac.norm.add, 2, quantile, probs = 0.975)
lCI.A <- apply(pac.norm.add, 2, quantile, probs = 0.025)
cols <- c("lightgreen", "mistyrose")[(mA > 1/NLinks) + 1] 
barplot2(mA, main = "PAco squared residuals - additive trees", xlab="Tribe", ylab="Normalized PACo sqr. residuals",
         cex.axis=0.5, col=cols, border="lightgrey", names.arg=colnamesPACo, las=2, cex.names=0.5, plot.ci=T, ci.l=lCI.A,
         ci.u=uCI.A, ci.color="blue")
abline(h=1/NLinks, col="red") 

#The above plot was suggested by the authors of the original script.
#We add the following ggplot2-based plot for our Brassicaceae study.

#Create a dataframe that contains data from all tribes.
PACo.results.phylogram.df<-as.data.frame(pac.norm.add)

#Update tribe names that were changed as a result of the above script.
colnames(PACo.results.phylogram.df)
colnames(PACo.results.phylogram.df)[colnames(PACo.results.phylogram.df)=="Outgroup.Cleomaceae"]<-"Cleomaceae (outgroup)"
colnames(PACo.results.phylogram.df)[colnames(PACo.results.phylogram.df)=="Arabidopsideae.trib..nov."]<-"Arabidopsideae trib. nov."
colnames(PACo.results.phylogram.df)[colnames(PACo.results.phylogram.df)=="Asperuginoideae.trib..nov."]<-"Asperuginoideae trib. nov."
colnames(PACo.results.phylogram.df)[colnames(PACo.results.phylogram.df)=="Hemilophieae.trib..nov."]<-"Hemilophieae trib. nov."
colnames(PACo.results.phylogram.df)[colnames(PACo.results.phylogram.df)=="Schrenkielleae.trib..nov."]<-"Schrenkielleae trib. nov."

#Change to long format for use in ggplot2.
PACo.results.phylogram.df$bootstrap<-1:nrow(PACo.results.phylogram.df)
PACo.results.phylogram.df.long<-melt(PACo.results.phylogram.df, id.vars = c("bootstrap"), variable.name = "tribe")

#Add information on main lineage.
PACo.results.phylogram.df.long$main_lineage<-metadata$Lineages_Nikolov_et_al[match(PACo.results.phylogram.df.long$tribe, metadata$Tribe)]
#There are still NAs in the dataframe, for family Cleomaceae.
unique(PACo.results.phylogram.df.long$main_lineage)
PACo.results.phylogram.df.long$main_lineage[which(is.na(PACo.results.phylogram.df.long$main_lineage))]<-""

#Add a value for ordering x-axis in ggplot later.
PACo.results.phylogram.df.long$main_lineage2<-NA
PACo.results.phylogram.df.long$main_lineage2[PACo.results.phylogram.df.long$main_lineage %in% c("basal", "")]<-0
PACo.results.phylogram.df.long$main_lineage2[PACo.results.phylogram.df.long$main_lineage=="I"]<-1
PACo.results.phylogram.df.long$main_lineage2[PACo.results.phylogram.df.long$main_lineage=="II"]<-2
PACo.results.phylogram.df.long$main_lineage2[PACo.results.phylogram.df.long$main_lineage=="III"]<-3
PACo.results.phylogram.df.long$main_lineage2[PACo.results.phylogram.df.long$main_lineage=="IV"]<-4
PACo.results.phylogram.df.long$main_lineage2[PACo.results.phylogram.df.long$main_lineage=="V"]<-5

#Add tribe name for plotting.
PACo.results.phylogram.df.long$tribe2<-paste0(PACo.results.phylogram.df.long$main_lineage2,"_",PACo.results.phylogram.df.long$tribe)

#Add information on outlier Y/N.
outliers_by_tribe<-as.data.frame(cbind(mA, cols))
rownames(outliers_by_tribe)[rownames(outliers_by_tribe)=="Outgroup.Cleomaceae"]<-"Cleomaceae (outgroup)"
rownames(outliers_by_tribe)[rownames(outliers_by_tribe)=="Arabidopsideae.trib..nov."]<-"Arabidopsideae trib. nov."
rownames(outliers_by_tribe)[rownames(outliers_by_tribe)=="Asperuginoideae.trib..nov."]<-"Asperuginoideae trib. nov."
rownames(outliers_by_tribe)[rownames(outliers_by_tribe)=="Hemilophieae.trib..nov."]<-"Hemilophieae trib. nov."
rownames(outliers_by_tribe)[rownames(outliers_by_tribe)=="Schrenkielleae.trib..nov."]<-"Schrenkielleae trib. nov."
PACo.results.phylogram.df.long$congruence<-outliers_by_tribe$cols[match(PACo.results.phylogram.df.long$tribe, rownames(outliers_by_tribe))]

#Create the ggplot.
ggplot(data=PACo.results.phylogram.df.long, aes(x=tribe2, y=value, fill=congruence)) +
  geom_boxplot(outlier.size=0.2)+
  scale_fill_manual(values=c("lightgreen", "pink"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), legend.position="none")+
  ylab("Normalized PACo sqr. residuals")+
  xlab("Tribe")+
  geom_hline(yintercept=1/NLinks, color="red", size=1)

#Save the plot for publication.
ggsave(last_plot(), 
       file = "12b.PACo_squared_residuals.svg",
       width = 10,
       height = 6)


## 6.5. Plot PACo results using unit branch length trees ----

mA <- apply(pac.norm.top, 2, median) ### OBS! ##### MEDIAN
uCI.A <- apply(pac.norm.top, 2, quantile, probs = 0.975)
lCI.A <- apply(pac.norm.top, 2, quantile, probs = 0.025)
cols <- c("lightgreen", "mistyrose")[(mA > 1/NLinks) + 1] 
barplot2(mA, main = "PAco squared residuals - unit branch length trees", xlab="Tribe", ylab="Normalized PACo sqr. residuals",
         cex.axis=0.5, col=cols, border="lightgrey", names.arg=colnamesPACo, las=2, cex.names=0.5, plot.ci=T, ci.l=lCI.A,
         ci.u=uCI.A, ci.color="blue")
abline(h=1/NLinks, col="red") 


## 6.6. Plot ParaFit2 results using additive trees ----

mA <- apply(pf2.add, 2, median) ### OBS! ##### MEDIAN
uCI.A <- apply(pf2.add, 2, quantile, probs = 0.975)
lCI.A <- apply(pf2.add, 2, quantile, probs = 0.025)
cols <- c("lightgreen", "mistyrose")[(mA > 0) + 1] 
barplot2(mA, main = "pfl2 statistic - additive trees", xlab="Tribe", ylab="Normalized PACo sqr. residuals",
         cex.axis=0.5, col=cols, border="lightgrey", names.arg=colnamesPACo, las=2, cex.names=0.5, plot.ci=T, ci.l=lCI.A,
         ci.u=uCI.A, ci.color="blue")
abline(h=0, col="red") 


## 6.7. Plot ParaFit2 using unit branch length trees ----

mA <- apply(pf2.top, 2,  median) ### OBS! ##### MEDIAN
uCI.A <- apply(pf2.top, 2, quantile, probs = 0.975)
lCI.A <- apply(pf2.top, 2, quantile, probs = 0.025)
cols <- c("lightgreen", "mistyrose")[(mA > 0) + 1] 
barplot2(mA, main = "pfl2 statistic - unit branch length trees", xlab="Tribe", ylab="Normalized PACo sqr. residuals",
         cex.axis=0.5, col=cols, border="lightgrey", names.arg=colnamesPACo, las=2, cex.names=0.5, plot.ci=T, ci.l=lCI.A,
         ci.u=uCI.A, ci.color="blue")
abline(h=0, col="red") 



#STEP 7: Validate terminal classification using Parition Around Medoids Algorithm (PAM) ----
## 7.1 normalize PACo sq residuals to total m2 ----

sum.pac.add <- apply(pac.add, 1, sum) # Additive trees
pac.add <- pac.add/sum.pac.add - 1/NLinks

sum.pac.top <- apply(pac.top, 1, sum) # Unit branch length trees
pac.top <- pac.top/sum.pac.top - 1/NLinks

## 7.4. incongruence metrics ----

im.paco.add <- apply(pac.add, 2, median) 
im.paco.top <- apply(pac.top, 2, median) 

im.pf2.add <- apply(pf2.add, 2, median)
im.pf2.top <- apply(pf2.top, 2, median)

## 7.5. Combine 2 metrics in one plot ----

par (mfrow=c(2,1))  
plot(im.paco.add,im.pf2.add)
plot(im.paco.top,im.pf2.top)

## 7.6. standardize metrics ----

x.paco.add <- mean(im.paco.add) ; x.pf2.add <- mean(im.pf2.add)
sd.paco.add<- sd(im.paco.add)  ; sd.pf2.add <- sd(im.pf2.add)
im.paco.stadd <- (x.paco.add - im.paco.add)/sd.paco.add
im.pf2.stadd <- (x.pf2.add - im.pf2.add)/sd.pf2.add
metrics.stadd <- data.frame(im.paco.stadd, im.pf2.stadd)

x.paco.top <- mean(im.paco.top) ; x.pf2.top <- mean(im.pf2.top)
sd.paco.top <- sd(im.paco.top)  ; sd.pf2.top <- sd(im.pf2.top)
im.paco.sttop <- (x.paco.top - im.paco.top)/sd.paco.top
im.pf2.sttop <- (x.pf2.top - im.pf2.top)/sd.pf2.top
metrics.sttop <- data.frame(im.paco.sttop, im.pf2.sttop)

## 7.7. execute clustering analyisis (PAM algorithm) ----

nclust = 2 # set no. clusters. Sometimes nclust=3 works better

## 7.8. Apply PAM with PACo and Parafit statistics using additive trees ----

K.PAM <- pam(metrics.stadd, nclust, diss=FALSE)
plot(im.paco.add,im.pf2.add, col=c("red","blue")[K.PAM$clustering])
title(main=list("PACo-Parafit - additive trees", cex=0.8))
SPaPf.add <- silhouette(K.PAM) # measures strength of separation
cat(summary(SPaPf.add)$avg.width)
SPaPf.add  <- summary(SPaPf.add)$avg.width
cat("\n")

## 7.9. Apply PAM with PACO statistic only using additive trees ----

K.PAM <- pam(metrics.stadd[1], nclust, diss=FALSE)
plot(im.paco.add,im.pf2.add, col=c("red","blue")[K.PAM$clustering])
title(main=list("PACo + additive trees", cex=0.8))
SPa.add <- silhouette(K.PAM) # measures strength of separation
cat(summary(SPa.add)$avg.width)
SPa.add <- summary(SPa.add)$avg.width
cat("\n")

## 7.10. Apply PAM with PACo and Parafit statistics using unit branch length trees ----

K.PAM <- pam(metrics.sttop, nclust, diss=FALSE)
plot(im.paco.top,im.pf2.top, col=c("red","blue")[K.PAM$clustering])
title(main=list("PACo-pf2 - unit branch length trees", cex=0.8))
SPaPf.top <- silhouette(K.PAM) # measures strength of separation
cat(summary(SPaPf.top)$avg.width)
SPaPf.top  <- summary(SPaPf.top)$avg.width
cat("\n")

## 7.11. Apply PAM with PACO statistic only using unit branch length trees ----

K.PAM <- pam(metrics.sttop[1], nclust, diss=FALSE)
plot(im.paco.top,im.pf2.top, col=c("red","blue")[K.PAM$clustering])
title(main=list("PACo - unit branch length trees", cex=0.8))
SPa.top <- silhouette(K.PAM) # measures strength of separation
cat(summary(SPa.top)$avg.width)
SPa.top <- summary(SPa.top)$avg.width
cat("\n")

Sall <- rbind(SPaPf.add, SPa.add, SPaPf.top, SPa.top)
rownames(Sall) <- c("Silhouette PACo-Parafit additive", "Silhouette PACo additive", 
                    "Silhouette PACo-Parafit unit branch length", "Silhouette PACo unit branch length")

## 7.12. write on disk all Silhoutte values ----

write.table(Sall, "results_PACo_analysis/PACo_output_text_files/Silhouette_values_all.txt")


#STEP 8: PLOT RESULTS ----

# 8. Plotting outlier terminals on trees
# 8.1 Set plotting parameters

op <- par(oma=c(3,2,1,1))
par (mfrow=c(1,2),mar = c(1,1,4,4))

# 8.2 Root trees if they are unrooted

NTree <- root(NTree, "Outgroup_Cleomaceae", resolve.root=TRUE)
CPTree <- root(CPTree, "Outgroup_Cleomaceae", resolve.root=TRUE)
NTree  <- ladderize(NTree, right=TRUE)
CPTree  <- ladderize(CPTree, right=TRUE)

# 8.3 Compute median values from squared residual values using additive trees
# Sqrd residual values are coded binary according to a threshold value 

mA <- apply(pac.norm.add, 2, median)
mA[mA > 1/NLinks] <- 1
mA[mA < 1/NLinks] <- 0
mA <- as.data.frame(mA)
labels <- strsplit(colnamesPACo, "_")
labels  <- unique(unlist(labels))
out <- mA$mA
names(out) <- labels
out

# 8.4 Plot potential outliers on tree set 1

plotTree(NTree, setEnv = T, offset=0.5, fsize=0.5, lwd=1)
title(main="Nuclear tree of Gene 1 - PACo potential ouliers", font.main=1, cex.main=0.8)
tiplabels(pie = to.matrix(out, sort(unique(out))), piecol = c("lightgreen", "lightcoral"), 
          cex = 0.5)
legend("bottomleft", c("Congruent", "Outlier"),
       cex=0.9, pch=16, col=c("lightgreen", "lightcoral"))

# 8.5 Plot potential outliers on tree set 2

plotTree(CPTree, setEnv = T, offset=0.5, fsize=0.5, lwd=1)
title(main="Chloroplast tree of Gene 2 - PACo potential ouliers", font.main=1, cex.main=0.8)
tiplabels(pie = to.matrix(out, sort(unique(out))), piecol = c("lightgreen", "lightcoral"), 
          cex = 0.5)


#STEP 9: PLOT COPHYLOGENY-STYLE OPPOSITE TREES ----
#note: script based on suggestions taken from https://yulab-smu.top/treedata-book/chapter2.html

#Create dataframes to store outlier data for each of the two trees.
out.df_n<-data.frame(tribe=tree_n$tip.label)
rownames(out.df_n)<-out.df_n$tribe
out.df_n$tribe2<-sub("_", " ", out.df_n$tribe) #twice to replace all underscores as for the out dataframe
out.df_n$tribe2<-sub("_", " ", out.df_n$tribe2)
out.df_n$outlier<-out[match(out.df_n$tribe2, names(out))]
out.df_n$main_lineage<-metadata$Lineages_Nikolov_et_al[match(out.df$tribe2, metadata$Tribe)]
out.df_n$main_lineage[is.na(out.df_n$main_lineage)|out.df_n$main_lineage=="basal"]<-rep("", sum(is.na(out.df_n$main_lineage)|out.df_n$main_lineage=="basal"))

out.df_cp<-data.frame(tribe=tree_cp$tip.label)
rownames(out.df_cp)<-out.df_cp$tribe
out.df_cp$tribe2<-sub("_", " ", out.df_cp$tribe) #twice to replace all underscores as for the out dataframe
out.df_cp$tribe2<-sub("_", " ", out.df_cp$tribe2)
out.df_cp$outlier<-out[match(out.df_cp$tribe2, names(out))]
out.df_cp$main_lineage<-metadata$Lineages_Nikolov_et_al[match(out.df$tribe2, metadata$Tribe)]
out.df_cp$main_lineage[is.na(out.df_cp$main_lineage)|out.df_cp$main_lineage=="basal"]<-rep("", sum(is.na(out.df_cp$main_lineage)|out.df_cp$main_lineage=="basal"))

#Create objects from ggtree and ggtree data.
cophylogeny_n<-ggtree(tree_n) %<+% out.df_n +
  geom_tiplab(aes(label=tribe2, color=as.character(outlier)), linetype = "dotted", align = T)+
  scale_color_manual(values=c("black", "red"))
cophylogeny_n<-gheatmap(cophylogeny_n, data = out.df_n[,c("main_lineage"), drop=F], width=.05, offset=0, colnames=F) +
  scale_fill_manual(values=c("orange", "lightblue","green","pink","yellow"), na.translate = F) +
  theme(legend.position="none")

  
cophylogeny_cp<-ggtree(tree_cp) %<+% out.df_cp +
  geom_tiplab(aes(label=tribe2, color=as.character(outlier)), linetype = "dotted", align = T)+
  scale_color_manual(values=c("black", "red"))
cophylogeny_cp<-gheatmap(cophylogeny_cp, data = out.df_n[,c("main_lineage"), drop=F], width=.05, offset=0, colnames=F) +
  scale_fill_manual(values=c("orange", "lightblue","green","pink","yellow"), na.translate = F) +
  theme(legend.position="none")

cophylogeny_n_data<-cophylogeny_n$data
cophylogeny_cp_data<-cophylogeny_cp$data

#Reverse x-axis and
#set offset to make the tree on the right-hand side of the first tree
cophylogeny_cp_data$x <- max(cophylogeny_cp_data$x) - cophylogeny_cp_data$x + max(cophylogeny_n_data$x) + 0.5

# pp <- cophylogeny_n + geom_tree(data=cophylogeny_cp_data) +
#   ggnewscale::new_scale_color() +
#   geom_tiplab(data = out.df_cp, aes(label=tribe2, color=as.character(outlier))) +
#   scale_color_manual(values=c("black", "red"))+
#   theme(legend.position="none")

dd <- bind_rows(cophylogeny_n_data, cophylogeny_cp_data) %>% 
  filter(!is.na(label))

cophylogeny_n_plus_lines<-cophylogeny_n + geom_line(aes(x, y, group=label), data=dd, color='grey')

#Combine plots using cowplot package.
ggsave(cowplot::plot_grid(cophylogeny_n_plus_lines,
                   cophylogeny_cp +
                     scale_x_reverse(),
                   rel_widths = c(2,1)),
       file = "12c.PACo_analysis_cophylogeny_plot.svg")

