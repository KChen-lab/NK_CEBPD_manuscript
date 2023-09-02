
################## Fig 5 G-I CEBPD expression and regulon visualization ###################

library(tidyverse)
library(RColorBrewer)
# atac seurat object for atac profiles 
# rna  seurat object for rna profiles 
DefaultAssay(atac) <- "chromvar"
p <- FeaturePlot(atac,
            features = "CEBPD", label = TRUE, repel = TRUE) +
            scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p <- FeaturePlot(rna, label.size=7,
            features = "CEBPD", label = TRUE, repel = TRUE,  max.cutoff=3)

pdf(file="./CEBPD_geneExpression.umap.pdf",width=6, height=5)
print(p)
dev.off()

p2 <- FeaturePlot(rna,
            features = "CEBPD_regulon1", label = TRUE, repel = TRUE) +
            scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#p2 <- VlnPlot(rna,features=c("CEBPD_regulon1"),pt.size=0, cols = vega_20_scanpy) + ylab("Regulon Score") + ggtitle("CEBPD regulon")
pdf(file="./CEBPD_regulon.umap.pdf",width=6, height=5)
print(p2)
dev.off()
