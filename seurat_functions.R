### functions for Seurat processing
options(echo=TRUE)
qc_plots_stats <- function(seur, basename, version){

        png(paste0(basename,"_",version,"_vlnQC_nodots.png"),height = 800,width=1100)
        print(VlnPlot(object = seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0,group.by = "orig.ident"))
        dev.off()
        
        png(paste0(basename,"_",version,"_vlnQC_withdots.png"),height = 800,width=1100)
        print(VlnPlot(object = seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident"))
        dev.off()
        
        plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        
        png(paste0(basename,"_",version,"_featurescatter.png"),height = 800,width=1100)
        print(CombinePlots(plots = list(plot1, plot2)))
        dev.off()
        
        png(paste0(basename,"_",version,"_percentmt_histogram.png"),height = 800,width=1100)
        hist(seur@meta.data$percent.mt)
        dev.off()
        png(paste0(basename,"_",version,"_nCountRNA_histogram.png"),height = 800,width=1100)
        hist(seur@meta.data$nCount_RNA)
        dev.off()
        png(paste0(basename,"_",version,"_nFeatureRNA_histogram.png"),height = 800,width=1100)
        hist(seur@meta.data$nFeature_RNA)
        dev.off()
        
        mean_percent.mt = mean(seur$percent.mt)
        median_percent.mt = median(seur$percent.mt)
        sd_percent.mt = sd(seur$percent.mt)
        mean_nCount = mean(seur$nCount_RNA)
        median_nCount= median(seur$nCount_RNA)
        sd_nCount = sd(seur$nCount_RNA)
        mean_nFeature = mean(seur$nFeature_RNA)
        median_nFeature= mean(seur$nFeature_RNA)
        sd_nFeature=sd(seur$nFeature_RNA)
        
        all_stats = list(mean_percent.mt, median_percent.mt, sd_percent.mt, mean_nCount, median_nCount, sd_nCount, mean_nFeature, median_nFeature, sd_nFeature)
        
        statsdf = as.data.frame(do.call(rbind, all_stats))
        stats_print = cbind(statsdf, "Stat"=c("Mean_Percent.mt","Median_Percent.mt", "Standard_deviation_Percent.mt","Mean_nCount","Median_nCount", "Standard_deviation_nCount","Mean_nFeature","Median_nFeature", "Standard_deviation_nFeature"))
        
        write.csv(stats_print, file=paste(basename,version,"stats.csv", sep="_"))
}

run_dr <- function(seur, basename, version, res, num_pcs){
        if(is.numeric(res)){
                seur <- FindClusters(seur, resolution = res)
                umap_plotting(seur, basename,version,res)
                qc_plotting(seur,basename,version)
                if(do_marks){
                        get_marks(seur,basename,version,x)
                }
                return(seur)
        }else{
                qc_plotting(seur, basename,version)
                res=c(0.4, 0.8, 1.2, 1.6)
                for(x in res){
                        seur <- FindClusters(seur, resolution = x)
                        umap_plotting(seur, basename, version, x)
                        if(do_marks){
                                get_marks(seur,basename,version,x)
                        }
                }
                return(seur)
        }
}

umap_plotting <- function(seur,basename,version,res){
        png(paste(basename, version, "res",res,"UMAP.png",sep="_"),height=800,width=1100)
        print(DimPlot(seur, reduction = "umap")) 
        dev.off()
}

qc_plotting <- function(seur,basename,version){
        
        genes_plot=c("nFeature_RNA","nCount_RNA","percent.mt")
        
        for(gene in genes_plot){
                png(paste(basename,version, gene,"plot.png",sep="_"),height = 800,width=1100)
                print(FeaturePlot(seur,features = gene,pt.size = 2))
                dev.off()
        }
        
}

get_marks <- function(seur,basename,version,res){
        seur.markers <- FindAllMarkers(seur, only.pos = TRUE, assay="RNA", min.pct = 0.25, logfc.threshold = 0.25)
        seur.markers.ordered <- arrange(seur.markers, cluster, desc(avg_logFC))
        
        write.csv(seur.markers.ordered,file=paste(basename,version,"res",res,"markers.csv",sep="_"))
}

#######################################
##check parameters

list_of_analyses <- function(analyses){
        loa =strsplit(analyses,",")[[1]]
        list_to_test =sapply(loa, function(x) gsub(" ", "", x, fixed = TRUE))
        return(list_to_test)
}

parameter_checks <- function(parameters){
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
        
        check.integer <- function(x){x%%1==0}
        
        assertthat::assert_that(horm %in% c("mouse","human"), 
                                msg="Parameters must contain 'mouse' or 'human' in first column")
        assertthat::assert_that(is.numeric(MT_filter) | is.na(MT_filter),
                                msg = "MT_filter must be numeric or NA")
        assertthat::assert_that(is.numeric(nCount_high) | is.na(nCount_high),
                                msg = "nCount_high must be numeric or NA")
        assertthat::assert_that(is.numeric(nFeature_high) | is.na(nFeature_high),
                                msg = "nFeature_high must be numeric or NA")
        assertthat::is.number(nCount_low)
        assertthat::is.number(nFeature_low)
        assertthat::assert_that(check.integer(num_pcs),
                                msg = "num_pcs must be an integer")
        assertthat::assert_that(is.numeric(res) | res == "range",
                                msg="res must be numeric or 'range'")
        assertthat::assert_that(is.character(doublets) | is.na(doublets),
                                msg="doublets must be path to doublet csv file or NA")
        list_to_test = list_of_analyses(analyses)
        assertthat::assert_that(sum(!list_to_test %in% c("all","basic","basic_regress_pcmt","basic_regress_cnt","basic_regress_both","sctransform","sctransform_regress","scran"))==0,
                                msg = "Analyses must be one or more of the following, separated by a comma: 'all', 'basic','basic_regress_pcmt', 'basic_regress_cnt','basic_regress_both','sctransform','sctransform_regress','scran'")
}
