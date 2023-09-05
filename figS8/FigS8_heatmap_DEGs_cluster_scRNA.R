

################## Fig S8 heatmap for cluster level DEGs (scRNA) #####################


####################  heatmap by cluster  
library(dplyr) 
DEGs <- readRDS(file="/rsrch1/bcb/kchen_group/data/Rezvani/SC122_scRNA_ATAC/result_May29_2022/DEGs_byClst_rna.RDS")
DEGs %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_logFC) -> top10
tp <-  subset(part_rmIL12, downsample = 100)
lst <-c()
for(i in seq(0,13,1)){
  lst <-c(lst, top10$gene[top10$cluster==i])
  if(i==2){
    lst<-c(lst,c("CISH","DUSP2","BAX","ATM","TNFRSF1B","TGFB1","KLRD1","PTPN22","PTN7","CD96"))
  }
  if(i==3){
    lst<-c(lst,c("CEBPD","CLIC3","GZMK","IRF1",'REST',"EOMES","E2F3","BIRC3","NFKBIA","ETS1"))
  }
}


library(ggplot2)
#lst <-c(top10$gene, )
p1 <- DoHeatmap(object =tp, features=lst, label=F,  disp.min=-2,disp.max=2, group.colors=vega_20_scanpy) + scale_color_manual(values=vega_20_scanpy)

tiff("./rna_clst_top10gene_heatmap.tiff", width=7, height =12, res =300, units = "in", compression = "lzw")
print(p1)
dev.off()






