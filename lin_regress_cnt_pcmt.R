### linear regression with both percent.mt and nCount_RNA

basic_regress_both <- function(seur_lin_both,basename,version,res,num_pcs,do_marks){
        dir.create("basic_regress_pcmt_ncount")
        setwd("basic_regress_pcmt_ncount/")
        
        all.genes <- rownames(seur_lin_both)
        seur_lin_both <- ScaleData(seur_lin_both, features = all.genes,vars.to.regress = c("percent.mt","nCount_RNA"))
        seur_lin_both <- RunPCA(seur_lin_both, features = VariableFeatures(object = seur_lin_both))
        
        png(paste(basename,version,"basic_regress_both_PCAelbow.png",sep="_"),height = 800,width=1100)
        print(ElbowPlot(seur_lin_both))
        dev.off()
        
        seur_lin_both <- RunUMAP(seur_lin_both, dims = 1:num_pcs)
        seur_lin_both <- FindNeighbors(seur_lin_both, dims = 1:num_pcs)
        
        seur_lin_both <- run_dr(seur_lin_both, basename, version, res, num_pcs) 
       
        saveRDS(seur_lin_both, file = paste(basename,version,"regress_both.rds",sep="_"))
        
        
        setwd("../")
}
