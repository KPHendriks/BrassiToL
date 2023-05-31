###This script reads in the SplitsTree4 network svg file and can change annotations.

#Read the SplitsTree4 network svg file.
#An svg file can be read and string-searched like any text file.
splits_network_svg  <- readLines("5a.Nuclear_genes_supermatrix_SplitsTree4.svg")
#SplitsTree4 changed "-2" in "_2"; change that back.
splits_network_svg <- gsub(pattern="_2", replace="-2", x=splits_network_svg)

#Load the annotations table from the previous script.
annotations<-read.csv(file = "2e.metadata.csv")

#Create a dataframe for replacing samples with Genus names.
replace<-annotations[,c(1:4,8)] #Make sure to include here the annotation by lineage cf Nikolov et al. (2019)
replace$Library_ID2 <- gsub(pattern=" ", replace="", x=replace$Library_ID)

#Create a dataframe to score how many times each Genus is represented in the network.
#This way we can add integer numbers to genera present >1 in the network.
genera_used<-data.frame(Genus=unique(replace$Genus), times_used=0)

#Loop though the possible samples in the network to replace their sample numbers with Genus names.
for (i in 1:nrow(replace)){
  #for (i in 1:5){
  #Check for each sample in the svg file.
  pattern = replace$Library_ID2[i]
  Genus = replace$Genus[i]
  #Check if sample i is used in the svg file at all.
  if (length(grep(pattern, splits_network_svg, ))>0){
    #If so, list in in the used genera dataframe.
    genera_used$times_used[genera_used$Genus==Genus]<-genera_used$times_used[genera_used$Genus==Genus]+1
  }
  #Now replace the sample name with the Genus in the svg file and add a number if already used.
  if (genera_used$times_used[genera_used$Genus==Genus]-1==0){
    splits_network_svg2  <- gsub(pattern = pattern, replace = Genus, x = splits_network_svg)
  } else {
    splits_network_svg2  <- gsub(pattern = pattern, replace = paste0(Genus,genera_used$times_used[genera_used$Genus==Genus]-1) , x = splits_network_svg)
  }
  splits_network_svg <- splits_network_svg2
}

#Write out a new svg file with the updated names.
writeLines(splits_network_svg, con="5d.Brassicaceae_nuclear_superstrict_supermatrix_renamed_Genus.svg")


#Repeat for Tribes.

#Read the original SplitsTree4 network svg file.
#An svg file can be read and string-searched like any text file.
splits_network_svg  <- readLines("5a.Nuclear_genes_supermatrix_SplitsTree4.svg")
#SplitsTree4 changed "-2" in "_2"; change that back.
splits_network_svg <- gsub(pattern="_2", replace="-2", x=splits_network_svg)

#Create a dataframe to score how many times each Tribe is represented in the network.
#This way we can add integer numbers to genera present >1 in the network.
Tribes_used<-data.frame(Tribe=unique(replace$Tribe), times_used=0)

#Loop though the possible samples in the network to replace their sample numbers with Tribe names.
for (i in 1:nrow(replace)){
  #for (i in 1:5){
  #Check for each sample in the svg file.
  pattern = replace$Library_ID2[i]
  Tribe = replace$Tribe[i]
  #Check if sample i is used in the svg file at all.
  if (length(grep(pattern, splits_network_svg, ))>0){
    #If so, list in in the used genera dataframe.
    Tribes_used$times_used[Tribes_used$Tribe==Tribe]<-Tribes_used$times_used[Tribes_used$Tribe==Tribe]+1
  }
  #Now replace the sample name with the Tribe in the svg file and add a number if already used.
  if (Tribes_used$times_used[Tribes_used$Tribe==Tribe]-1==0){
    splits_network_svg2  <- gsub(pattern = pattern, replace = Tribe, x = splits_network_svg)
  } else {
    splits_network_svg2  <- gsub(pattern = pattern, replace = paste0(Tribe,Tribes_used$times_used[Tribes_used$Tribe==Tribe]-1) , x = splits_network_svg)
  }
  splits_network_svg <- splits_network_svg2
}

#Write out a new svg file with the updated names.
writeLines(splits_network_svg, con="5d.Brassicaceae_nuclear_superstrict_supermatrix_renamed_Tribe.svg")



#Repeat for Lineages cf. Nikolov et al. (2019).

#Read the original SplitsTree4 network svg file.
#An svg file can be read and string-searched like any text file.
splits_network_svg  <- readLines("5a.Nuclear_genes_supermatrix_SplitsTree4.svg")
#SplitsTree4 changed "-2" in "_2"; change that back.
splits_network_svg <- gsub(pattern="_2", replace="-2", x=splits_network_svg)

#Create a dataframe to score how many times each lineage is represented in the network.
#This way we can add integer numbers to genera present >1 in the network.
lineages_used<-data.frame(lineage=unique(replace$Lineages_Nikolov_et_al), times_used=0)

#Loop though the possible samples in the network to replace their sample numbers with lineage names.
for (i in 1:nrow(replace)){
  #for (i in 1:5){
  #Check for each sample in the svg file.
  pattern = replace$Library_ID2[i]
  lineage = replace$Lineages_Nikolov_et_al[i]
  #Check if sample i is used in the svg file at all.
  if (length(grep(pattern, splits_network_svg, ))>0){
    #If so, list in in the used genera dataframe.
    lineages_used$times_used[lineages_used$lineage==lineage]<-lineages_used$times_used[lineages_used$lineage==lineage]+1
  }
  #Now replace the sample name with the lineage in the svg file and add a number if already used.
  splits_network_svg2  <- gsub(pattern = pattern, replace = lineage, x = splits_network_svg)

  splits_network_svg <- splits_network_svg2
}

#Write out a new svg file with the updated names.
writeLines(splits_network_svg, con="5d.Brassicaceae_nuclear_superstrict_supermatrix_renamed_lineage.svg")
