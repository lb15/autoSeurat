# autoSeurat

Scripts for processing cellranger output through Seurat. Works with Seurat V4.

Submission script is written for UCSF's Wynton HPC submission system.

Requires cellranger output, including the .mtx cell x gene matrix, the .tsv  barcode list, and the .tsv gene list.
Parameters for seurat processing are specified in a parameters.csv file. Right now only tested with "basic_analysis".

Output is a folder (named by the variable version) containing QC plots, UMAP plots, and .rds Seurat object.

Scripts to process data with SoupX are in the Soupx_scripts folder. Works with SoupX v1.5.0


