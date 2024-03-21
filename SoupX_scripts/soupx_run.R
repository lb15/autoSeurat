library(SoupX)
library(Matrix)

args=commandArgs(trailingOnly=T)
basename=args[1]
version=args[2]
log_file=paste0(basename,"_",version,"_soupx_log.txt")

sink(file=log_file)
print(paste0("SoupX Processing for ", basename, " version ", version))
sink()

setwd(paste0("/wynton/group/reiter/lauren/",basename,"/"))

#sc=load10X("outs/")

toc = Seurat::Read10X("outs/filtered_feature_bc_matrix/")
tod = Seurat::Read10X("outs/raw_feature_bc_matrix/")

if(!is.na(version)){
	sink(file=log_file,append=T)
	print(paste0("Susbetting to version ",version," cell barcodes"))
	sink()
	seur=readRDS(paste0(version,"/basic_analysis/",basename,"_",version,"_basic.rds"))
	toc_red=toc[,match(colnames(seur),colnames(toc))]
	clusters=seur$seurat_clusters
}else{
	toc_red=toc
	clus = read.csv("outs/analysis/clustering/graphclust/clusters.csv")
	clusters=clus$Cluster
	names(clusters) = clus$Barcode
}
sc = SoupChannel(tod, toc_red)
sc=setClusters(sc, clusters)
sink(file=log_file,append=T)
sc=autoEstCont(sc)
sink()

out=adjustCounts(sc)

dir.create("soupx_results")

writeMM(out, file="soupx_results/matrix.mtx")
write(colnames(out), file="soupx_results/barcodes.tsv")
write(rownames(out),file="soupx_results/features.tsv")

saveRDS(sc, file="soupx_results/soupx.rds")
