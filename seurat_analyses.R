### seurat analyses

###################################################################
### basic, no linear regression
options(echo=TRUE)
basic_noregress <-function(seur,basename,version,res,num_pcs,do_marks){
        
        dir.create("basic_analysis")
        setwd("basic_analysis/")
        
        all.genes <- rownames(seur)
        seur <- ScaleData(seur, features = all.genes)
        seur <- RunPCA(seur, features = VariableFeatures(object = seur))
        
        png(paste(basename,version,"basic_PCAelbow.png",sep="_"),height = 800,width=1100)
        print(ElbowPlot(seur,ndims=40))
        dev.off()
        
        seur <- RunUMAP(seur, dims = 1:num_pcs)
        seur <- FindNeighbors(seur, dims = 1:num_pcs)
        
        seur <- run_dr(seur, basename, version, res, num_pcs)  
        
        saveRDS(seur, file = paste(basename,version,"basic.rds",sep="_"))
        
        setwd("../")
}
##########################################################
## linear regression with percent.mt

basic_regress_pcmt <- function(seur_lin_reg,basename,version,res,num_pcs,do_marks){
        dir.create("basic_regress_pcmt")
        setwd("basic_regress_pcmt/")
        
        all.genes <- rownames(seur_lin_reg)
        seur_lin_reg <- ScaleData(seur_lin_reg, features = all.genes,vars.to.regress = "percent.mt")
        seur_lin_reg <- RunPCA(seur_lin_reg, features = VariableFeatures(object = seur_lin_reg))
        
        png(paste(basename,version,"basic_regress_PCAelbow.png",sep="_"),height = 800,width=1100)
        print(ElbowPlot(seur_lin_reg,ndims=40))
        dev.off()
        
        seur_lin_reg <- RunUMAP(seur_lin_reg, dims = 1:num_pcs)
        seur_lin_reg <- FindNeighbors(seur_lin_reg, dims = 1:num_pcs)
        
        seur_lin_reg <- run_dr(seur_lin_reg, basename, version, res, num_pcs) 
        
        saveRDS(seur_lin_reg, file = paste(basename,version,"regress_pcmt.rds",sep="_"))
        
        setwd("../")
}
### linear regression with both percent.mt and nCount_RNA

basic_regress_both <- function(seur_lin_both,basename,version,res,num_pcs,do_marks){
        dir.create("basic_regress_pcmt_ncount")
        setwd("basic_regress_pcmt_ncount/")
        
        all.genes <- rownames(seur_lin_both)
        seur_lin_both <- ScaleData(seur_lin_both, features = all.genes,vars.to.regress = c("percent.mt","nCount_RNA"))
        seur_lin_both <- RunPCA(seur_lin_both, features = VariableFeatures(object = seur_lin_both))
        
        png(paste(basename,version,"basic_regress_both_PCAelbow.png",sep="_"),height = 800,width=1100)
        print(ElbowPlot(seur_lin_both,ndims=40))
        dev.off()
        
        seur_lin_both <- RunUMAP(seur_lin_both, dims = 1:num_pcs)
        seur_lin_both <- FindNeighbors(seur_lin_both, dims = 1:num_pcs)
        
        seur_lin_both <- run_dr(seur_lin_both, basename, version, res, num_pcs) 
       
        saveRDS(seur_lin_both, file = paste(basename,version,"regress_both.rds",sep="_"))
        
        
        setwd("../")
}

### linear regression with  nCount_RNA

basic_regress_count <- function(seur_regresscnt,basename,version,res,num_pcs,do_marks){
        dir.create("basic_regress_pcmt_ncount")
        setwd("basic_regress_pcmt_ncount/")
        
        all.genes <- rownames(seur_regresscnt)
        seur_regresscnt <- ScaleData(seur_regresscnt, features = all.genes,vars.to.regress = c("nCount_RNA"))
        seur_regresscnt <- RunPCA(seur_regresscnt, features = VariableFeatures(object = seur_regresscnt))
        
        png(paste(basename,version,"basic_regress_nCount_PCAelbow.png",sep="_"),height = 800,width=1100)
        print(ElbowPlot(seur_regresscnt,ndims=40))
        dev.off()
        
        seur_regresscnt<- RunUMAP(seur_regresscnt, dims = 1:num_pcs)
        seur_regresscnt <- FindNeighbors(seur_regresscnt, dims = 1:num_pcs)
        
        seur_regresscnt <- run_dr(seur_regresscnt, basename, version, res, num_pcs) 
        
        saveRDS(seur_regresscnt, file = paste(basename,version,"regress_nCount.rds",sep="_"))
   
        setwd("../")
}
##########################################################
##sctransform
## need to adjust default assay for get_marks -> not SCT assay
sctrans <- function(seur_sctrans,basename,version,res,num_pcs,do_marks){
        
        dir.create("sctransform")
        setwd("sctransform")
        
        seur_sctrans = SCTransform(seur_sctrans)
        seur_sctrans <- RunPCA(seur_sctrans, verbose = FALSE)
        seur_sctrans <- RunUMAP(seur_sctrans, dims = 1:num_pcs, verbose = FALSE)
        seur_sctrans <- FindNeighbors(seur_sctrans, dims = 1:num_pcs, verbose = FALSE)
        
        seur_sctrans <- run_dr(seur_sctrans, basename, version, res, num_pcs) 
        
        saveRDS(seur_sctrans, file = paste(basename,version,"sctrans.rds",sep="_"))
        
        setwd("../")
}
###############################################################
##sctransform with linear regression

sctrans_regress <- function(seur_sctrans,basename,version,res,num_pcs,do_marks){
        
        dir.create("sctransform_regress")
        setwd("sctransform_regress")
        
        seur_sctrans = SCTransform(seur_sctrans, vars.to.regress = "percent.mt")
        seur_sctrans <- RunPCA(seur_sctrans, verbose = FALSE)
        seur_sctrans <- RunUMAP(seur_sctrans, dims = 1:num_pcs, verbose = FALSE)
        seur_sctrans <- FindNeighbors(seur_sctrans, dims = 1:num_pcs, verbose = FALSE)

        seur_sctrans <- run_dr(seur_sctrans, basename, version, res, num_pcs) 
        
        saveRDS(seur_sctrans, file = paste(basename,version,"sctrans_pcmtregress.rds",sep="_"))
        
        setwd("../")
}
#######################################################################
## scran

scran_analysis <- function(seur_scran,basename,version,res,num_pcs,do_marks){
        library(scran)
        dir.create("scran")
        setwd("scran")
        
        sce <- SingleCellExperiment(assays = list(counts = as.matrix(x = seur_scran[["RNA"]]@counts))) # read data from Seurat
        clusters = quickCluster(sce, min.size=100)
        sce = computeSumFactors(sce, cluster=clusters)
        ## update normalize function to logNormCounts from scater
        #sce = normalize(sce, return_log = FALSE) # without(!) log transform
        sce= logNormCounts(sce, log=F)
        seur_scran@misc[["seurat_norm_data"]] = as.matrix(x = seur[["RNA"]]@data) # backup Seurat's norm data
        
        
        seur_scran[["RNA"]]@data = log(x = assay(sce, "normcounts") + 1)
        seur_scran <- FindVariableFeatures(seur_scran, selection.method = "vst", nfeatures = 2000)
        all.genes <- rownames(seur_scran)
        
        seur_scran <- ScaleData(seur, features=all.genes)
        seur_scran <- RunPCA(seur_scran, features = VariableFeatures(object = seur_scran))
        png(paste(basename,version,"scran_PCAelbow.png",sep="_"),height = 800,width=1100)
        print(ElbowPlot(seur_scran,ndims=40))
        dev.off()
        
        seur_scran <- RunUMAP(seur_scran, dims = 1:num_pcs, verbose = FALSE)
        seur_scran <- FindNeighbors(seur_scran, dims = 1:num_pcs, verbose = FALSE)
        
        seur_scran <- run_dr(seur_scran, basename, version, res, num_pcs) 
        
        saveRDS(seur_scran, file = paste(basename,version,"scran.rds",sep="_"))
        
        setwd("../")
        
        #### with regression for percent_mt
        
        dir.create("scran_pcmt")
        setwd("scran_pcmt")
        
        seur_scran <- ScaleData(seur, features=all.genes,vars.to.regress = "percent.mt")
        seur_scran <- RunPCA(seur_scran, features = VariableFeatures(object = seur_scran))
        
        png(paste(basename,version,"scran_regress_PCAelbow.png",sep="_"),height = 800,width=1100)
        print(ElbowPlot(seur_scran,ndims=40))
        dev.off()
        
        seur_scran <- RunUMAP(seur_scran, dims = 1:num_pcs, verbose = FALSE)
        seur_scran <- FindNeighbors(seur_scran, dims = 1:num_pcs, verbose = FALSE)
        
        seur_scran <- run_dr(seur_scran, basename, version, res, num_pcs) 
        
        saveRDS(seur_scran, file = paste(basename,version,"scran_regress.rds",sep="_"))
        
        setwd("../")
}
