# Global Brassicaceae phylogeny based on filtering of 1,000 gene dataset

These are the scripts associated with the publication of Hendriks et al. (2023) in <i>Current Biology</i>.
With these scripts it is possible to reconstruct the most complete Brassicaceae genus-level family phylogeny (Brassicaceae Tree of Life, or BrassiToL) to date, and analyse the impact of possible paralogy on the species tree.

## Highlights

- An unparalleled Brassicaceae phylogeny based on more than 1,000 genes is presented, covering nearly all 349 genera
- Cytonuclear discordance is omnipresent and a likely sign of rampant hybridisation
- A more recent (late Eocene to late Oligocene) origin of the family follows from fossil calibration
- Results are used to come up with a new family classification

## Summary
The mustard family (Brassicaceae) is a scientifically and economically important family, containing the model plant <i>Arabidopsis thaliana</i> and numerous crop species that feed billions worldwide. Despite its relevance, most phylogenetic trees of the family are incompletely sampled and often contain poorly supported branches. Here, we present the most complete Brassicaceae genus-level family phylogenies to date (Brassicaceae Tree of Life, or BrassiToL) based on nuclear (1,081 genes, 319 of 349 genera; 57 of 58 tribes) and plastome (60 genes, 265 genera; all tribes) data. We found cytonuclear discordance between the two, which is likely a result of rampant hybridisation among closely and more distantly related lineages. To evaluate the impact of such hybridisation on the nuclear phylogeny reconstruction, we performed five different gene sampling routines, which increasingly removed putatively paralog genes. Our cleaned subset of 297 genes revealed high support for the tribes, while support for the main lineages (supertribes) was moderate. Calibration based on the 20 most clock-like nuclear genes suggests a late Eocene to late Oligocene origin of the family. Finally, our results strongly support a recently published new family classification, dividing the family into two subfamilies (one with five supertribes), together representing 58 tribes. This includes five recently described or re-established tribes, including Arabidopsideae, a monogeneric tribe accommodating <i>Arabidopsis</i> without any close relatives. With a worldwide community of thousands of researchers working on Brassicaceae and its diverse members, our new genus-level family phylogeny will be an indispensable tool for studies on biodiversity and plant biology.

## Instructions

The shell and R scripts to start analyses are all within the main directory. The subdirectories hold either input data, functions, reference files, and resulting figures and tables. Directories that will contain (intermediate) results have been pre-created, but are empty to save on the amount of storage on GitHub. The scripts will create the revelant files in the process.

The various scripts have commonly been written for use on an HPC that allows parallel computing. Whenever you start using these scripts, make sure to first install all relevant tools and to update the scripts to direct to the location of these tools on you specific machine. 

[![DOI](https://zenodo.org/badge/647656354.svg)](https://zenodo.org/badge/latestdoi/647656354)
