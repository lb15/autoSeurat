## Seurat V3 analysis
## initialize object

library(Seurat)
library(dplyr)


### load in functions
source("/wynton/group/reiter/lauren/autoSeurat/seurat_functions.R")
source("/wynton/group/reiter/lauren/autoSeurat/seurat_analyses.R")

args=commandArgs(trailingOnly = T)

basename=args[1]
version=args[2]
parameters_file = args[3]
do_marks = args[4]

setwd(basename)

### load in parameter arguments
parameters = read.csv(parameters_file,sep=",",header=T)
parameter_checks(parameters)

horm = parameters$human_mouse
MT_filter = parameters$Mt_filter
nCount_high = parameters$nCount_high
nCount_low = parameters$nCount_low
nFeature_high = parameters$nFeature_high
nFeature_low = parameters$nFeature_low
num_pcs = parameters$PCs
res=parameters$resolution
doublets = parameters$path_to_doublet_csv 
analyses=parameters$analyses

## read in data
data_dir ="filtered_feature_bc_matrix/"
data = Read10X(data.dir=data_dir)
log_file = paste(basename,version,analyses,"log.txt",sep="_")

## create version folder and set directory
if(!dir.exists(version)){
        dir.create(version)
        setwd(version)
        sink(file=log_file)
        print("Created new directory")
        
}else{
        setwd(version)
        sink(file=log_file)
        print("Version directory already exists")
        
        
}

if(!dir.exists("QCplots")){
        dir.create("QCplots")
        print("Created new QCplots directory")
        
}

print("initializing Seurat object")
## initialize seurat object
seur <- CreateSeuratObject(counts = data, project = basename, min.cells = 3, min.features = 200)

if(horm == "mouse"){
        seur[["percent.mt"]] <- PercentageFeatureSet(object = seur, pattern = "^mt-")
        print("Using pattern '^mt-' for mouse")
} else{
        seur[["percent.mt"]] <- PercentageFeatureSet(object = seur, pattern = "^MT-")
        print("Using pattern '^MT-' for human")
}

########################################################################
## to subset or not to subset?
print("Subsetting based on filters")
if(sum(!is.na(c(MT_filter, nCount_high,nFeature_high))) == 0 & sum(nCount_low,nFeature_low)==0){
        print("No Subsetting")
}else if(sum(!is.na(c(MT_filter, nCount_high,nFeature_high))) == 0){
        seur=subset(seur, subset = nCount_RNA > nCount_low & nFeature_RNA > nFeature_low)
        print("Subsetted only on nCount_low and nFeature_low")
}else if(!is.na(MT_filter)){
        seur = subset(seur, subset = percent.mt <MT_filter & nCount_RNA > nCount_low & nFeature_RNA > nFeature_low)
        if(!is.na(nCount_high) & !is.na(nFeature_high)){
                seur=subset(seur, subset = nCount_RNA < nCount_high)
                seur=subset(seur, subset=nFeature_RNA < nFeature_high)
                print("Subsetted on percent.MT, nCount, and nFeature")
        }else if(!is.na(nCount_high)){
                seur=subset(seur, subset = nCount_RNA < nCount_high)
                print("Subsetted on percent.MT, nCount, nFeature_low")
        }else if(!is.na(nFeature_high)){
                seur=subset(seur, subset=nFeature_RNA < nFeature_high)
                print("Subsetted on percent.MT, nFeature, and nCount_low")
        }else{
                print("Subsetted on percent.MT, nCount_low, and nFeature_low")
        }
}else if(!is.na(nCount_high) & !is.na(nFeature_high)){
        seur=subset(seur, subset = nCount_RNA < nCount_high & nCount_RNA > nCount_low)
        seur=subset(seur, subset=nFeature_RNA < nFeature_high & nFeature_RNA > nFeature_low)
        print("Subsetted on nCount and nFeature")
}else if(!is.na(nCount_high)){
        seur=subset(seur, subset = nCount_RNA < nCount_high & nCount_RNA > nCount_low & nFeature_RNA > nFeature_low)
        print("Subsetted on nCount nad nFeature_low")
}else{
        seur=subset(seur, subset=nFeature_RNA < nFeature_high & nFeature_RNA > nFeature_low & nCount_RNA > nCount_low)
        print("Subsetted on nFeature and nCount_low")
}
                
sink()

setwd("QCplots/")
qc_plots_stats(seur, basename, version)
setwd("../")

#######################################################################
## taking care of doublets (everyday!)

if(is.na(doublets)){
        sink(file=log_file,append=T)
        print("Not removing any doublets")
        sink()
}else{
	print(doublets)
	seur <- add_doublets(seur,doublets)
	#plot_doublets(seur, basename, version) cannot plot because no dimension reduction yet run
	sink(file=log_file, append=T)
	print(paste0("Doublets removed:",table(seur$predicted_doublet)))
	sink()

	seur <- remove_doublets(seur)

	dir.create("QCplots_nodubs")
	setwd("QCplots_nodubs")
	qc_plots_stats(seur, basename, version)
	setwd("../")
}


#######################################################################
print("starting normalization")
seur_for_sctrans <- seur
seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seur), 10)
plot1_vg <- VariableFeaturePlot(seur)
plot2_vg <- LabelPoints(plot = plot1_vg, points = top10, repel = TRUE)

png(paste0("QCplots/",basename,"_",version,"_variablegenes.png"),height = 800,width=1100)
CombinePlots(plots = list(plot1_vg, plot2_vg))
dev.off()

######################################################################
## Seurat analyses

sink(log_file,append=T)
print(paste0("Analyses to perform: ", analyses))
sink()

if("all" %in% analyses){
        sink(file=log_file,append=T,split=T)
        print("Running all analyses")
        basic_noregress(seur,basename,version,res,num_pcs,do_marks)
        print("Finished 'basic'")
        basic_regress_pcmt(seur,basename,version,res,num_pcs,do_marks)
        print("Finished 'basic_regress_pcmt'")
        basic_regress_both(seur,basename,version,res,num_pcs,do_marks)
        print("Finished 'basic_regress_both'")
        basic_regress_count(seur,basename,version,res,num_pcs,do_marks)
        print("Finished 'basic_regress_cnt'")
        sctrans(data,basename,version,res,num_pcs,do_marks)
        print("Finished 'sctransform'")
        sctrans_regress(data,basename,version,res,num_pcs,do_marks)
        print("Finished 'sctransform_regress'")
        scran_analysis(seur,basename,version,res,num_pcs,do_marks)
        print("Finished 'scran'")
        sink()
}else for(analysis in list_of_analyses(analyses)){
        sink(file=log_file,append=T,split=T)
        if(analysis == "basic"){
                print("Starting 'Basic' analysis")
                basic_noregress(seur,basename,version,res,num_pcs,do_marks)
                print("Finished 'Basic' analysis")
        }else if(analysis =="basic_regress_pcmt"){
                print("Starting 'Basic_regress_pcmt' analysis")
                basic_regress_pcmt(seur,basename,version,res,num_pcs,do_marks)
                print("Finished 'Basic_regress_pcmt' analysis")
        }else if(analysis == "basic_regress_cnt"){
                print("Starting 'Basic_regress_cnt' analysis")
                basic_regress_count(seur,basename,version,res,num_pcs,do_marks)
                print("Finished 'Basic_regress_cnt' analysis")
        }else if(analysis == "basic_regress_both"){
                print("Starting 'Basic_regress_both' analysis")
                basic_regress_both(seur,basename,version,res,num_pcs,do_marks)
                print("Finished'Basic_regress_both' analysis")
        }else if(analysis == "sctransform"){
                print("Starting 'sctransform' analysis")
                sctrans(seur_for_sctrans,basename,version,res,num_pcs,do_marks)
                print("Finished 'sctransform' analysis")
        }else if(analysis == "sctransform_regress"){
                print("Starting 'sctransform_regress' analysis")
                sctrans_regress(seur_for_sctrans,basename,version,res,num_pcs,do_marks)
                print("Finished 'sctransform_regress' analysis")
        }else if(analysis == "scran"){
                print("Starting 'scran' analysis")
                scran_analysis(seur,basename,version,res,num_pcs,do_marks)
                print("Finished 'scran' analysis")
        }
        sink()
}


