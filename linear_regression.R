library(Seurat)

##########################################################
## linear regression with percent.mt

basic_regress_pcmt <- function(seur_lin_reg,basename,version,res,num_pcs,do_marks){
        dir.create("basic_regress_pcmt")
        setwd("basic_regress_pcmt/")
        
        all.genes <- rownames(seur_lin_reg)
        seur_lin_reg <- ScaleData(seur_lin_reg, features = all.genes,vars.to.regress = "percent.mt")
        seur_lin_reg <- RunPCA(seur_lin_reg, features = VariableFeatures(object = seur_lin_reg))
        
        png(paste(basename,version,"basic_regress_PCAelbow.png",sep="_"),height = 800,width=1100)
        print(ElbowPlot(seur_lin_reg))
        dev.off()
        
        seur_lin_reg <- RunUMAP(seur_lin_reg, dims = 1:num_pcs)
        seur_lin_reg <- FindNeighbors(seur_lin_reg, dims = 1:num_pcs)
        
        seur_lin_reg <- run_dr(seur_lin_reg, basename, version, res, num_pcs) 
        
        saveRDS(seur_lin_reg, file = paste(basename,version,"regress_pcmt.rds",sep="_"))
        
        setwd("../")
}
