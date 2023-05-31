#!/bin/bash

parallel -j 47 '~/iqtree-2.1.3-Linux/bin/iqtree2 -s results_assessment_saturation_and_clocklikeliness/input_gene_alignments_cleaned/{}_cleaned.fasta -B 1000 -nm 1000 -m GTR+F+R -nt AUTO --prefix results_assessment_saturation_and_clocklikeliness/results_iqtree_gene_trees/result_{}_iqtree' :::: 4c.gene_alignments_cleaned_list.txt
