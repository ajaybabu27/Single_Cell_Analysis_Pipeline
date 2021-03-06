Tree_Roll
----------

Description:

Probabilistic matching of single cell gene expression vectors to appropriate mouse Anatomical Structures guided by EMAP-GXD DAG Tree. This code is mainly written in python and it also serves as a wrapper for couple of R scripts that are used for statistical calculations. 

Main Command: 

tree_roll_as_matching_clean.py arg1

1) Here arg1 is Theiler Stage (eg. 'TS13')

2) Edit file paths on the top of main body of code to set location of various input files. The default file path for single cell gene expression dataset uses subset of data for quick demonstration (Runs in about 20 minutes (10x bootstraps), Full data with 1000X bootstraps runs for about 4 hours). 

Folder Structure:

experiment_data - contains the single cell boolean gene expression dataset (qPCR data)
reference_data - contains file with DAG relationships (*.obo) and GXD gene expression data filtered for 96 qPCR markers (biomark_96_*)
ouput - contains results from the main script. Contains two main outputs:
	1) Graph file showing enrichment of various single cells to different anatomical strucutures. This is shown as nodes shaded by a color gradient.
	Each node is annotated in the following format: anatomical strucuture name_num of genes annotated_(number of matching single cells)
	2) TSV file showing cells assigned to anatomical structures at each level.
	

Other Scripts:

create_Robject.R - Called by main script to create a background model from single cell gene expression data.
tree_matching.R - Called by main script to perform statistical analysis by computing p-value (by bootstrapping from background model) across different nodes at a given level in the tree.
