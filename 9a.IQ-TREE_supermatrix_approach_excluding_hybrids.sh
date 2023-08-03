#!/bin/bash

#Create directories to store results.
mkdir -p results_IQ-TREE_supermatrix_excluding_hybrids/
mkdir -p results_IQ-TREE_supermatrix_excluding_hybrids/input_tree_files/
mkdir -p results_IQ-TREE_supermatrix_excluding_hybrids/input_multiple_sequence_alignments/
mkdir -p results_IQ-TREE_supermatrix_excluding_hybrids/output_concatenated_sequences/
mkdir -p results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/
mkdir -p results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/

#Get the latest multiple sequence alignments for concatenation into a supermatrix.
cp results_superstrict_excluding_hybrids/results_taper_final_superstrict_excluding_hybrids/*.fasta results_IQ-TREE_supermatrix_excluding_hybrids/input_multiple_sequence_alignments/
#Remove a single temp file that might have ended up in the mafft sequence alignment directory.
rm results_IQ-TREE_supermatrix_excluding_hybrids/input_multiple_sequence_alignments/temp*
#Remove empty files that had no sequences left after cleaning loops from step 3.
find results_IQ-TREE_supermatrix_excluding_hybrids/input_multiple_sequence_alignments/ -size 0 -print -delete


#Create a list of fasta files.
ls results_IQ-TREE_supermatrix_excluding_hybrids/input_multiple_sequence_alignments/*taper*.fasta > list.all.genes

#Use the cat sequences tool (https://github.com/ChrisCreevey/catsequences) to concatenate the multiple sequence alignments into a single supermatrix.
~/catsequences/catsequences list.all.genes

# #Remove unwanted gene alignment files for later use in calculation of concordance factors (see this script further down).
# comm -1 -3 list.genes.to.keep list.all.genes > list.genes.to.remove
# xargs rm < list.genes.to.remove

#Move the results to a directory with a sensible name.
mv allseqs.fas results_IQ-TREE_supermatrix_excluding_hybrids/output_concatenated_sequences/
mv allseqs.partitions.txt results_IQ-TREE_supermatrix_excluding_hybrids/output_concatenated_sequences/
mv list* results_IQ-TREE_supermatrix_excluding_hybrids/output_concatenated_sequences/

#Infer a concatenation-based species tree with 1000 ultrafast bootstrap and an edge-linked partition model.
~/iqtree-2.2.2-Linux/bin/iqtree2 -s results_IQ-TREE_supermatrix_excluding_hybrids/output_concatenated_sequences/allseqs.fas -o "MYZV" -m GTR+F+R --prefix results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/ -bb 1000 -wbtl -nt AUTO #-wbtl makes sure to write bootstrap trees including branch lengths.
#iqtree2 -p $ALIGN_FOLDER --prefix "$TREE_FOLDER/concat_$NAME" -bb 1000 -nt $THREADS
#~/iqtree-2.1.3-Linux/bin/iqtree2 -s ${gene}_taper_preliminary_strict.fasta -B 1000 -nm 1000 -m GTR+F+R -nt AUTO --prefix ${gene}_iqtree_preliminary_strict

#Let's rename some files with more sensible names. Note that IQ-TREE stored hidden files, so we directly make these visible.
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/.iqtree results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/iqtree_report.txt
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/.treefile results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/iqtree_ML_tree.tree
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/.mldist results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/iqtree_ML_distances.txt
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/.splits.nex results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/iqtree_bootstrap_support_values.nex
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/.contree results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/iqtree_bootstrap_concensus_tree.tree
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/.ufboot results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/iqtree_bootstrap_bootstrapped_trees.tree
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/.log results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/iqtree_bootstrap_logfile.log

#IQ-TREE writes out hidden files, so let's make all remaining hidden files visible.
find results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/ -type f -name '.*' -execdir sh -c 'mv -i "$0" "./${0#./.}"' {} \;

#Calculate concordance factors.

#First get the subset of gene trees calculated before for a the coalescent approach for the 'inclusive' routine.
cp results_superstrict_excluding_hybrids/results_iqtree_final_superstrict_excluding_hybrids/*.treefile \
  results_IQ-TREE_supermatrix_excluding_hybrids/input_tree_files/

#Create a list of gene tree files.
ls results_IQ-TREE_supermatrix_excluding_hybrids/input_tree_files/*.treefile > list.all.gene.trees

#Remove tree files as done for gene alignment files above.

#Create a list of gene trees.
ls results_IQ-TREE_supermatrix_excluding_hybrids/input_tree_files/*.treefile > list.all.gene.trees

# grep -v "4757\|6110\|6016\|6620\|5200\|6780\|4889\|5366\|7602\|5404\|6072\|5578\|6738\|6946\|6376\|5921\|5816\|5168\|6532\|5502\|5318\|5177\|5348\|7141\|6420\|5664\|5770\|6148\|6176\|4932\|5280\|6412\|5614\|7174\|7136\|5206\|5064\|6498\|5131\|6450\|4691\|5138\|6459\|5702\|6914\|5599\|6487\|5913\|7067\|6979\|5339\|6048\|5870\|6387\|6492\|6034\|7363\|7111\|7241\|6050\|5942\|5910\|5220\|5163\|6882\|5463\|6652\|5434\|6859\|4724\|5477\|6528\|6114\|5090\|5326\|6639\|4954\|6404\|5990\|5716\|6298\|5460\|4527\|6238\|5562\|6270\|4802\|6221\|5802\|5945\|6026\|6500\|6265\|6641\|5639\|6258\|5357\|6000\|6363\|5899\|5894\|5427\|5335\|6496\|5347\|6164\|6825\|6667\|6003\|6860\|5660\|5821\|5422\|5454\|6955\|5840\|6439\|6909\|5536\|5644\|6494\|5421\|5018\|5842\|5941\|6958\|6139\|5594\|6550\|6405\|6533\|6128\|5703\|6785\|6004\|5841\|5430\|6198\|6226\|5858\|6679\|4942\|5859\|5958\|5299\|6565\|5865\|5815\|6130\|6460\|5699\|6098\|5328\|7336\|5950\|5936\|5744\|5918\|5449\|6119\|6875\|6462\|4890\|6051\|6282\|6320\|4989\|5032\|5981\|6570\|6913\|6557\|6961\|4796\|5188\|5123\|6563\|6036\|6886\|6447\|7628\|4806\|7029\|5843\|6883\|5343\|5772\|6064\|6685\|5791\|5620\|5428\|6430\|6544\|5333\|5354\|6733\|4793\|6636\|6373\|5893\|5949\|5977\|5853\|6401\|6746\|6274\|5940\|5469\|6947\|6284\|6038\|5944\|5919\|5551\|5554\|5355\|4951\|5426\|6449\|6068\|5857\|4893\|5116\|5464\|6507\|6175\|6954\|6782\|5406\|5721\|7128\|6318\|6295\|5257\|5296\|1G01910\|1G02400\|1G02420\|1G02860\|1G02970\|1G03100\|1G03390\|1G04110\|1G04690\|1G04910\|1G04970\|1G05850\|1G05910\|1G06260\|1G06270\|1G06440\|1G06560\|1G06950\|1G07590\|1G08460\|1G09010\|1G09340\|1G09870\|1G09900\|1G10020\|1G10330\|1G10385\|1G10780\|1G11915\|1G12150\|1G12330\|1G12470\|1G13970\|1G15440\|1G15740\|1G16070\|1G16570\|1G19025\|1G19860\|1G20080\|1G20560\|1G22620\|1G24610\|1G25500\|1G25570\|1G27752\|1G27760\|1G29690\|1G29900\|1G30000\|1G30440\|1G33410\|1G47380\|1G47670\|1G48090\|1G48880\|1G49820\|1G49970\|1G50120\|1G50200\|1G50940\|1G51550\|1G51720\|1G52520\|1G52760\|1G53270\|1G55150\|1G55270\|1G55480\|1G55760\|1G55910\|1G58120\|1G60060\|1G60490\|1G60995\|1G61850\|1G63160\|1G63660\|1G63700\|1G65320\|1G65840\|1G66830\|1G67120\|1G67560\|1G68570\|1G68740\|1G68830\|1G70260\|1G70570\|1G70820\|1G72280\|1G72970\|1G73180\|1G73930\|1G74460\|1G74680\|1G75200\|1G76400\|1G78280\|1G79080\|1G80350\|1G80360\|1G80410\|1G80460\|2G01070\|2G02560\|2G04660\|2G13370\|2G13540\|2G15620\|2G15695\|2G16440\|2G16880\|2G17020\|2G17510\|2G17760\|2G18710\|2G20050\|2G20190\|2G21470\|2G21610\|2G22120\|2G23140\|2G25800\|2G26180\|2G26780\|2G26930\|2G27090\|2G27500\|2G27590\|2G28070\|2G28780\|2G29050\|2G30520\|2G31530\|2G32040\|2G32590\|2G32900\|2G33770\|2G35040\|2G35150\|2G35610\|2G35720\|2G36840\|2G37230\|2G37310\|2G37320\|2G37370\|2G38000\|2G38280\|2G38500\|2G38510\|2G38770\|2G39260\|2G39450\|2G39620\|2G39730\|2G39830\|2G39970\|2G40070\|2G40190\|2G40730\|2G40760\|2G40890\|2G41190\|2G42490\|2G42700\|2G42850\|2G43235\|2G43420\|2G43770\|2G43890\|2G45500\|2G46060\|2G46370\|2G46890\|2G47020\|2G47790\|3G01100\|3G01460\|3G01670\|3G02130\|3G02570\|3G03210\|3G03380\|3G03710\|3G04480\|3G04740\|3G06270\|3G06483\|3G06530\|3G06810\|3G06860\|3G06880\|3G07080\|3G07140\|3G07420\|3G08670\|3G08760\|3G08950\|3G08960\|3G09090\|3G10030\|3G10230\|3G10380\|3G11220\|3G11460\|3G11540\|3G11960\|3G11964\|3G12280\|3G12610\|3G13490\|3G15380\|3G15550\|3G16270\|3G16910\|3G17640\|3G17810\|3G17830\|3G17900\|3G18290\|3G18524\|3G19553\|3G19630\|3G19810\|3G19970\|3G19990\|3G20050\|3G20240\|3G20260\|3G20630\|3G20920\|3G21420\|3G21720\|3G22170\|3G23020\|3G23430\|3G24170\|3G25430\|3G25660\|3G26090\|3G26410\|3G26670\|3G26700\|3G27820\|3G28040\|3G29320\|3G46220\|3G47400\|3G47700\|3G48150\|3G48380\|3G48610\|3G48820\|3G50620\|3G50660\|3G51470\|3G51490\|3G51930\|3G52200\|3G52970\|3G53150\|3G53180\|3G54690\|3G54790\|3G54860\|3G55070\|3G56120\|3G56310\|3G56370\|3G56940\|3G57180\|3G57300\|3G57610\|3G60830\|3G61960\|4G00450\|4G00550\|4G00740\|4G00910\|4G01037\|4G01320\|4G01570\|4G01730\|4G01880\|4G02030\|4G02260\|4G02460\|4G02680\|4G02730\|4G02780\|4G02990\|4G03020\|4G03220\|4G03260\|4G04640\|4G05090\|4G08920\|4G11120\|4G15420\|4G16144\|4G16570\|4G16660\|4G17090\|4G17230\|4G17830\|4G18010\|4G18060\|4G18340\|4G18830\|4G18950\|4G19180\|4G19460\|4G20050\|4G20060\|4G20070\|4G20090\|4G21350\|4G21680\|4G21770\|4G22720\|4G24190\|4G24220\|4G24610\|4G24620\|4G24740\|4G24790\|4G24830\|4G26750\|4G27640\|4G28080\|4G28220\|4G29010\|4G29310\|4G29810\|4G30490\|4G30600\|4G30790\|4G31200\|4G31540\|4G31770\|4G31850\|4G32430\|4G32910\|4G33210\|4G33330\|4G33410\|4G33440\|4G33760\|4G33945\|4G34260\|4G34350\|4G34850\|4G35560\|4G35760\|4G35870\|4G35880\|4G36390\|4G36630\|4G36790\|4G37030\|4G38010\|4G39850\|5G01360\|5G01380\|5G02270\|5G02480\|5G02810\|5G03070\|5G03280\|5G03430\|5G04360\|5G04480\|5G04550\|5G04710\|5G04930\|5G05200\|5G05680\|5G05840\|5G08030\|5G08280\|5G08310\|5G08415\|5G08490\|5G08530\|5G08580\|5G08660\|5G10900\|5G11390\|5G11560\|5G11960\|5G12040\|5G12290\|5G12470\|5G12900\|5G13020\|5G13030\|5G13230\|5G13420\|5G13560\|5G13630\|5G13640\|5G13770\|5G14210\|5G14230\|5G14450\|5G14510\|5G14700\|5G14950\|5G15300\|5G15400\|5G15540\|5G15710\|5G15880\|5G16280\|5G16300\|5G16310\|5G17230\|5G17250\|5G17620\|5G17860\|5G18580\|5G18670\|5G18700\|5G18830\|5G19130\|5G19180\|5G19350\|5G19640\|5G19680\|5G19690\|5G19730\|5G20270\|5G20350\|5G20890\|5G20990\|5G21930\|5G21970\|5G21990\|5G22030\|5G22260\|5G22330\|5G22350\|5G22450\|5G22510\|5G22800\|5G22850\|5G23300\|5G23340\|5G24240\|5G24810\|5G25230\|5G26570\|5G30510\|5G36880\|5G38530\|5G39040\|5G39500\|5G40405\|5G40740\|5G43920\|5G44000\|5G44240\|5G44370\|5G46070\|5G46100\|5G46220\|5G46580\|5G47010\|5G47690\|5G47780\|5G48270\|5G48385\|5G48800\|5G49030\|5G49070\|5G49430\|5G49890\|5G51150\|5G51690\|5G51970\|5G52280\|5G52580\|5G52850\|5G53000\|5G53320\|5G53480\|5G53920\|5G53970\|5G54180\|5G54260\|5G54690\|5G54830\|5G55540\|5G55960\|5G56220\|5G56290\|5G56580\|5G56730\|5G57800\|5G58020\|5G58230\|5G58330\|5G58410\|5G58430\|5G58470\|5G58510\|5G58750\|5G58940\|5G59900\|5G60020\|5G60450\|5G60750\|5G60790\|5G60980\|5G61400\|5G61560\|5G61865\|5G62130\|5G62310\|5G62370\|5G63010\|5G63050\|5G63410\|5G63640\|5G63770\|5G63840\|5G63860\|5G64270\|5G64370\|5G65500\|5G66680\|5G67430" list.all.gene.trees > list.gene.trees.to.keep
# 
# #Remove unwanted gene tree files for later use in calculation of concordance factors (see this script further down).
# comm -1 -3 list.gene.trees.to.keep list.all.gene.trees > list.gene.trees.to.remove
# xargs rm < list.gene.trees.to.remove
# 
# #Remove gene tree files we want to exclude.
# mv list* results_IQ-TREE_supermatrix_excluding_hybrids/input_tree_files/

#Concatenate all treefiles into a single input file for IQ-TREE.
cat results_IQ-TREE_supermatrix_excluding_hybrids/input_tree_files/*.treefile > results_IQ-TREE_supermatrix_excluding_hybrids/input_tree_files/gene_trees_combined.tree

#Start IQ-TREE to run calculation of concordance factors from species tree + gene trees + gene alignments.
#We followed suggestions from http://www.iqtree.org/doc/Concordance-Factor.
~/iqtree-2.2.2-Linux/bin/iqtree2 \
  -t results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/iqtree_ML_tree.tree \
  --gcf results_IQ-TREE_supermatrix_excluding_hybrids/input_tree_files/gene_trees_combined.tree \
  -p results_IQ-TREE_supermatrix_excluding_hybrids/input_multiple_sequence_alignments/ \
  --scf 1000 \
  --prefix results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/ \
  -nt 55 ##nb cannot be set to AUTO when calculating concordance factors
# ~/iqtree-2.1.3-Linux/bin/iqtree2 -t results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree/iqtree_ML_tree.tree --gcf results_IQ-TREE_supermatrix_excluding_hybrids/input_tree_files/gene_trees_combined.tree -p results_IQ-TREE_supermatrix_excluding_hybrids/input_multiple_sequence_alignments/ --scf 1000 --prefix results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/ -nt 1

#Like before, IQ-TREE writes out hidden files. Let's rename these with sensible names and make them visible.
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/.cf.tree results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/iqtree_ML_tree_with_concordance_factors.tree
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/.cf.tree.nex results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/iqtree_ML_tree_with_concordance_factors.nex
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/.cf.branch results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/iqtree_ML_tree_with_branch_IDs.tree
mv results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/.cf.stat results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/iqtree_ML_tree_concordance_factors_per_branch.stat

#IQ-TREE writes out hidden files, so let's make all remaining hidden files visible.
find results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/ -type f -name '.*' -execdir sh -c 'mv -i "$0" "./${0#./.}"' {} \;

#Save a copy of the annotated supermatrix gene tree for publication.
cp results_IQ-TREE_supermatrix_excluding_hybrids/output_iqtree_condordance_factors/iqtree_ML_tree_with_concordance_factors.tree 9b.IQ-TREE_supermatrix_supermatrix_excluding_hybrids_with_concordance_factors.tree

#Use a short R script to root the phylogeny.
#Rscript 10b.treePL_root_phylogeny.R

#Run treePL and save output from screen to a log file for later review.
#[This script assumes the use of treePL on the Naturalis HighMem computing cluster, which requires sudo to run.]
#sudo treePL 10c.treePL_config > results_IQ-TREE_supermatrix_excluding_hybrids/output_treePL/10d.treePL_calibrated_phylogeny.log
