# autoSeurat

Scripts for processing cellranger output through Seurat.

Submission script is written for UCSF's Wynton SGE submission system.

Requires cellranger output, including the .mtx cell x gene matrix, the .tsv  barcode list, and the .tsv gene list.
Parameters for seurat processing are specified in a parameters.csv file.

Output is a folder (named by the variable version) containing QC plots, UMAP plots, and .rds Seurat object.

Scripts to process data with Soupx are in the Soupx_scripts folder.
