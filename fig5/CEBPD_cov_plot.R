


library(Signac)
library(Seurat)

ranges.show <- StringToGRanges("chr8-47738000-47738350") 
ranges.show$color <- "gray20" 

Idents(obj) <- obj@meta.data$ID
obj@meta.data$ID <- factor(obj@meta.data$ID, levels=c("IL15_baseline","IL21_baseline",  "IL15_GBM_3","IL21_GBM_3","IL15_GBM_9","IL21_GBM_9"))

cov_plot<-CoveragePlot(
        object = obj,
        assay = "peaks",
        region =r,
        group.by=c("ID"),
        annotation = TRUE,
        peaks = TRUE,
        region.highlight = ranges.show, 
        ranges.title = r,
        links = FALSE,
        idents = c("IL15_GBM_3", "IL21_GBM_3", "IL15_GBM_9", "IL21_GBM_9"),
        extend.upstream = 400,
        extend.downstream = 600
      ) & 
  scale_fill_manual(values = rep(c("#84B761", "#C9503E"), 2)) 
 
cov_plot[[1]]$labels$title <- "CEBPD" 
cov_plot[[1]]$theme$plot.title <- element_text(size = 20, face = "bold", hjust = 0.5) 
cov_plot[[1]]$theme$strip.text.y.left <- element_text(size = 16, face = "bold", angle = 0) 
cov_plot[[1]]$theme$axis.text.x <- element_blank() 
 
#cov_plot[[1]]$layers[[4]]$aes_params$size <- 6 
cov_plot[[1]]$theme$axis.title.x <- element_text(size = 16) 
cov_plot[[1]]$theme$axis.text.x <- element_text(size = 12) 

