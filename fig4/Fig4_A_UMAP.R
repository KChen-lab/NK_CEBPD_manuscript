############################################################################
##########   May 29 2022  scATAC-seq  ######################################
############################################################################


vega_20 = c(
    '#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
    '#8C564B', '#E377C2', '#BCBD22', '#17BECF', '#AEC7E8',
    '#FFBB78', '#98DF8A', '#FF9896', '#C5B0D5', '#C49C94',
    '#F7B6D2', '#DBDB8D', '#9EDAE5', '#AD494A', '#8C6D31')
library(Seurat)

##########   Fig 4A UMAP for scATAC-seq: IL12, IL15 on 3 timepoints ########

part_rmIL12 <- readRDS(file="/rsrch1/bcb/kchen_group/data/Rezvani/SC122_scRNA_ATAC/result_May29_2022/atac_IL21_15.RDS")
part_rmIL12@meta.data$product <- part_rmIL12@meta.data$orig.ident
part_rmIL12@meta.data$ID  <- paste0(part_rmIL12@meta.data$product,"_", part_rmIL12@meta.data$time)
Idents(part_rmIL12) <-part_rmIL12@meta.data$seurat_clusters

pdf(file="~/project/nk_multiOme/atac_umap_il12_il15.pdf", width=8, height=6)
p <- DimPlot(part_rmIL12,ncol=3,split.by=c("ID"),cols=vega_20,label=T)
print(p)
dev.off()

Idents(part_rmIL12) <- part_rmIL12@meta.data$product 

Idents(part_rmIL12) <- part_rmIL12@meta.data$orig.ident
pdf(file="~/project/nk_multiOme/rna_umap_il12_il15_merge.pdf", width=6, height=5)
p <- DimPlot(part_rmIL12,cols=vega_20[10:20],label=F)
print(p)
dev.off()
