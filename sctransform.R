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
