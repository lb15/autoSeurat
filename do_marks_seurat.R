## read in pre-made Seurat object and perform DE expression analysis

library(Seurat)
library(dplyr)

args=commandArgs(trailingOnly=TRUE)
basename=args[1]
version=args[2]
parameters_file = args[3]

setwd(basename)

### load in parameter arguments
parameters = read.csv(parameters_file,sep=",",header=T)

res=parameters$resolution
doublets = parameters$path_to_doublet_csv
analyses=parameters$analyses
setwd(version)

for(analysis in analyses){
	
	seur=readRDS(paste0(basename,"_",version,"_",analysis,".rds"))

	
	seur.markers <- FindAllMarkers(seur, only.pos = TRUE, assay="RNA", min.pct = 0.25, logfc.threshold = 0.25,group.by=paste0("RNA_snn_res.",res))
	seur.markers.ordered <- arrange(seur.markers, cluster, desc(avg_logFC))
	
	if(analysis == "basic"){
		write.csv(seur.markers.ordered,file=paste0(analysis,"_analysis/",basename,"_",version,"_res_",res,"_markers.csv"))
	}else{
		write.csv(seur.markers.ordered,file=paste0(analysis,"/",basename,"_",version,"_res_",res,"_markers.csv"))  	
	}
}

